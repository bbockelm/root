#include <ROOT/RDF/Utils.hxx>
#include <ROOT/RRootBulkDS.hxx>
#include <ROOT/TSeq.hxx>
#include <TBufferFile.h>
#include <TClass.h>
#include <TDataType.h>
#include <TError.h>
#include <TROOT.h>         // For the gROOTMutex
#include <TVirtualMutex.h> // For the R__LOCKGUARD
#include <ROOT/RMakeUnique.hxx>
#include <ROOT/TBulkBranchRead.hxx>

#include <algorithm>
#include <vector>

namespace ROOT {

namespace Internal {

namespace RDF {

/// A struct to handle the state of each column in a given slot.
//
//  We assume that the TBulkBufferMgr is accessed from a single thread and optimize
//  access for scanning monotonically, one-at-a-time through the corresponding entries.
class TBulkBufferMgr
{
private:
   // This is the mapping from col # to void* address.
   std::vector<void*> fAddressMap;

   // This is the mapping from col # to buffer.
   std::vector<TBufferFile*> fBufferMap;

   // This is the mapping from col # to data type.
   std::vector<EDataType> fDataTypeMap;

   // These are all the buffers that must be advanced by n bytes for each event.
   // 1 byte for bool/char, 2 bytes for short, 4 bytes for int/float and 8 bytes for double/long64
   std::vector<TBufferFile*> fOneByteBuffers;
   std::vector<TBufferFile*> fTwoByteBuffers;
   std::vector<TBufferFile*> fFourByteBuffers;
   std::vector<TBufferFile*> fEightByteBuffers;
   // A list of n-byte values that are the targets of the void* pointer handed back to
   // the RDF.  Each event advances n bytes as the above buffers.
   std::vector<Char_t*>      fOneByteValues;
   std::vector<Short_t*>     fTwoByteValues;
   std::vector<Int_t*>       fFourByteValues;
   std::vector<Long64_t*>    fEightByteValues;

   // The current absolute entry in the TTree.
   ULong64_t fCurAbsEntry{0};
   // The current entry relative to the cluster beginning.
   ULong64_t fCurRelEntry{0};

   // The current TTree* being processed by the buffer manager.
   TTree *fCurTree{nullptr};

public:
   TBulkBufferMgr(TTree *curTree)
    : fCurTree(curTree)
   {
      const TObjArray *branchList = curTree->GetListOfBranches();
      Int_t branchCount = branchList->GetEntriesFast();
      for (auto idx : ROOT::TSeqU(branchCount)) {
         fAddressMap.push_back(nullptr);
         fBufferMap.push_back(nullptr);
         auto br = dynamic_cast<TBranch*>((*branchList)[idx]);
         if (!br || !br->SupportsBulkRead()) {
            printf("Skipping branch %s as it does not support bulk reads.\n", br->GetName());
            continue;
         }
         EDataType dt;
         TClass *cls = nullptr;
         if (br->GetExpectedType(cls, dt)) {
            printf("Skipping branch %s as we failed to retrieve the expected type info.\n", br->GetName());
            continue;
         }
         fDataTypeMap.push_back(dt);
         if (dt == kChar_t || dt == kUChar_t || dt == kBool_t) {
            TBufferFile *bf = new TBufferFile(TBuffer::kWrite, 32*1024);
            fOneByteBuffers.push_back(bf);
            fBufferMap.back() = fOneByteBuffers.back();
            Char_t *value = new Char_t[1];
            fOneByteValues.push_back(value);
            fAddressMap.back() = fOneByteValues.back();
         } else if (dt == kShort_t || dt == kUShort_t) {
            TBufferFile *bf = new TBufferFile(TBuffer::kWrite, 32*1024);
            fTwoByteBuffers.push_back(bf);
            fBufferMap.back() = fTwoByteBuffers.back();
            Short_t *value = new Short_t[1];
            fTwoByteValues.push_back(value);
            fAddressMap.back() = fTwoByteValues.back();
         } else if (dt == kFloat_t || dt == kInt_t || dt == kUInt_t) {
            TBufferFile *bf = new TBufferFile(TBuffer::kWrite, 32*1024);
            fFourByteBuffers.push_back(bf);
            fBufferMap.back() = fFourByteBuffers.back();
            Int_t *value = new Int_t[1];
            fFourByteValues.push_back(value);
            fAddressMap.back() = fFourByteValues.back();
         } else if (dt == kDouble_t || dt == kLong64_t || dt == kULong64_t) {
            TBufferFile *bf = new TBufferFile(TBuffer::kWrite, 32*1024);
            fEightByteBuffers.push_back(bf);
            fBufferMap.back() = fEightByteBuffers.back();
            Long64_t *value = new Long64_t[1];
            fEightByteValues.push_back(value);
            fAddressMap.back() = fEightByteValues.back();
         } else {
            printf("Skipping branch %s as its data type (%d) is not supported.\n", br->GetName(), dt);
         }
      }
   }

   ~TBulkBufferMgr() {
      for (auto idx : ROOT::TSeqU(fBufferMap.size())) {
         if (fBufferMap[idx]) delete fBufferMap[idx];
      }
      for (auto idx : ROOT::TSeqU(fAddressMap.size())) {
         if (fAddressMap[idx]) {
            if (fDataTypeMap[idx] == kChar_t || fDataTypeMap[idx] == kUChar_t || fDataTypeMap[idx] == kBool_t) {
               delete static_cast<Char_t*>(fAddressMap[idx]);
            } else if (fDataTypeMap[idx] == kShort_t || fDataTypeMap[idx] == kUShort_t) {
               delete static_cast<Short_t*>(fAddressMap[idx]);
            } else if (fDataTypeMap[idx] == kFloat_t || fDataTypeMap[idx] == kInt_t || fDataTypeMap[idx] == kUInt_t) {
               delete static_cast<Int_t*>(fAddressMap[idx]);
            } else if (fDataTypeMap[idx] == kDouble_t || fDataTypeMap[idx] == kLong64_t || fDataTypeMap[idx] == kULong64_t) {
               delete static_cast<Long64_t*>(fAddressMap[idx]);
            } else {
               printf("Unknown type\n");
            }
         }
      }
   }

   void *getColumnTargetPtr(size_t idx) {
      return &fAddressMap[idx];
   }

   bool SetEntry(ULong64_t entry) {
      // TODO: handle random skips.
      if (R__unlikely(fCurAbsEntry != entry)) {
          return false;
      }
      auto idx1b = 0;
      auto idx2b = 0;
      auto idx4b = 0;
      auto idx8b = 0;
      for (UInt_t idx = 0; idx < fBufferMap.size(); ++idx) {
         EDataType dt = fDataTypeMap[idx];
         if (dt == kChar_t || dt == kUChar_t || dt == kBool_t) {
            Char_t *raw_buffer = reinterpret_cast<Char_t*>(fBufferMap[idx]->GetCurrent());
            Char_t tmp = *reinterpret_cast<Char_t*>(&raw_buffer[fCurRelEntry]);
            char *tmp_ptr = reinterpret_cast<char *>(&tmp);
            frombuf(tmp_ptr, fOneByteValues[idx1b]);
            idx1b++;
         } else if (dt == kShort_t || dt == kUShort_t) {
            Short_t *raw_buffer = reinterpret_cast<Short_t*>(fBufferMap[idx]->GetCurrent());
            Short_t tmp = *reinterpret_cast<Short_t*>(&raw_buffer[fCurRelEntry]);
            char *tmp_ptr = reinterpret_cast<char *>(&tmp);
            frombuf(tmp_ptr, fTwoByteValues[idx2b]);
            idx2b++;
         } else if (dt == kFloat_t || dt == kInt_t || dt == kUInt_t) {
            Int_t *raw_buffer = reinterpret_cast<Int_t*>(fBufferMap[idx]->GetCurrent());
            Int_t tmp = *reinterpret_cast<Int_t*>(&raw_buffer[fCurRelEntry]);
            char *tmp_ptr = reinterpret_cast<char *>(&tmp);
            frombuf(tmp_ptr, fFourByteValues[idx4b]);
            idx4b++;
         } else if (dt == kDouble_t || dt == kLong64_t || dt == kULong64_t) {
            Long64_t *raw_buffer = reinterpret_cast<Long64_t*>(fBufferMap[idx]->GetCurrent());
            Long64_t tmp = *reinterpret_cast<Long64_t*>(&raw_buffer[fCurRelEntry]);
            char *tmp_ptr = reinterpret_cast<char *>(&tmp);
            frombuf(tmp_ptr, fEightByteValues[idx8b]);
            idx8b++;
         } else {
            printf("Unknown data type %d.\n", dt);
         }
      }
      fCurRelEntry++;
      fCurAbsEntry++;
      return true;
   }

   // Initialize a cluster range for processing.
   // Returns true if successful
   // On success, sets entry count to the number of available entries.
   bool SetEntryRange(ULong64_t firstEntry, Int_t &entryCount)
   {
      const TObjArray *branchList = fCurTree->GetListOfBranches();
      Int_t branchCount = branchList->GetEntriesFast();
      Int_t count = -1;
      for (auto idx : ROOT::TSeqU(branchCount)) {
         if (fBufferMap[idx] == nullptr) continue;

         auto br = static_cast<TBranch*>((*branchList)[idx]);
         auto result = br->GetBulkRead().GetEntriesSerialized(firstEntry, *fBufferMap[idx]);
         // TODO: this fails if all the baskets in the cluster do not have the same size.
         if (result < 0) return false;
         else if ((count >= 0) && (result != count)) return false;
         else count = result;
      }
      fCurAbsEntry = firstEntry;
      fCurRelEntry = 0;
      entryCount = count;
      return true;
   }
};

} // ns RDF

} // ns Internal

namespace RDF {

namespace Experimental {

// Return a list of type-erased pointers, one per slot.
// The target of the pointer is the location where we will deserialize the column's value for
// the given slot.
std::vector<void *> RRootBulkDS::GetColumnReadersImpl(std::string_view name, const std::type_info &id)
{
   const auto colTypeName = GetTypeName(name);
   const auto &colTypeId = ROOT::Internal::RDF::TypeName2TypeID(colTypeName);
   if (id != colTypeId) {
      std::string err = "The type of column \"";
      err += name;
      err += "\" is ";
      err += colTypeName;
      err += " but a different one has been selected.";
      throw std::runtime_error(err);
   }

   const auto index =
      std::distance(fListOfBranches.begin(), std::find(fListOfBranches.begin(), fListOfBranches.end(), name));

   std::vector<void *> ret(fNSlots);
   for (auto slot : ROOT::TSeqU(fNSlots)) {
      ret[slot] = static_cast<void *>(fBufferMgrs[slot]->getColumnTargetPtr(index));
   }
   return ret;
}

RRootBulkDS::RRootBulkDS(std::string_view treeName, std::string_view fileNameGlob)
   : fTreeName(treeName), fFileNameGlob(fileNameGlob), fModelChain(std::string(treeName).c_str())
{
   fModelChain.Add(fFileNameGlob.c_str());

   const TObjArray &lob = *fModelChain.GetListOfBranches();
   fListOfBranches.resize(lob.GetEntries());

   TIterCategory<TObjArray> iter(&lob);
   std::transform(iter.Begin(), iter.End(), fListOfBranches.begin(), [](TObject *o) { return o->GetName(); });
}

RRootBulkDS::~RRootBulkDS()
{
}

std::string RRootBulkDS::GetTypeName(std::string_view colName) const
{
   if (!HasColumn(colName)) {
      std::string e = "The dataset does not have column ";
      e += colName;
      throw std::runtime_error(e);
   }
   // TODO: we need to factor out the routine for the branch alone...
   // Maybe a cache for the names?
   auto typeName =
      ROOT::Internal::RDF::ColumnName2ColumnTypeName(std::string(colName), /*nsID=*/0, &fModelChain, /*ds=*/nullptr,
                                                     /*isCustomCol=*/false);
   // We may not have yet loaded the library where the dictionary of this type is
   TClass::GetClass(typeName.c_str());
   return typeName;
}

const std::vector<std::string> &RRootBulkDS::GetColumnNames() const
{
   return fListOfBranches;
}

bool RRootBulkDS::HasColumn(std::string_view colName) const
{
   if (!fListOfBranches.empty())
      GetColumnNames();
   return fListOfBranches.end() != std::find(fListOfBranches.begin(), fListOfBranches.end(), colName);
}

void RRootBulkDS::InitSlot(unsigned int slot, ULong64_t firstEntry)
{
   Int_t eventCount;
   if (!fBufferMgrs[slot]->SetEntryRange(firstEntry, eventCount)) {
      throw std::runtime_error("Failed to initialize slot");
   }
   // TODO: compare eventCount against calculated event ranges -- should be the same sizes!
}


void RRootBulkDS::InitialiseSlot(unsigned int slot)
{
   auto chain = new TChain(fTreeName.c_str());
   chain->ResetBit(kMustCleanup);
   chain->Add(fFileNameGlob.c_str());
   chain->GetEntry(0);
   fBufferMgrs[slot].reset(new ROOT::Internal::RDF::TBulkBufferMgr(chain->GetTree()));
   fChains[slot].reset(chain);
}


std::vector<std::pair<ULong64_t, ULong64_t>> RRootBulkDS::GetEntryRanges()
{
   if (fReturnedFirstRange) return {};
   // TODO: Improve this to allow for multi-file chains.

   Long64_t clusterStart = 0;
   fModelChain.GetEntry(0);
   auto tree = fModelChain.GetTree();
   auto clusterIter = tree->GetClusterIterator(clusterStart);
   std::vector<std::pair<ULong64_t, ULong64_t>> entryRanges;

   while ( (clusterStart = clusterIter()) < tree->GetEntries() ) {
      entryRanges.emplace_back(clusterStart, clusterIter.GetNextEntry());
   }

   fReturnedFirstRange = true;
   return entryRanges;
}


bool RRootBulkDS::SetEntry(unsigned int slot, ULong64_t entry)
{
   //printf("RRootBulkDS::SetEntry slot %u, entry %llu\n", slot, entry);
   bool result = fBufferMgrs[slot]->SetEntry(entry);
   //printf("Value: %d\n", **static_cast<int**>(fBufferMgrs[slot]->getColumnTargetPtr(0)));
   return result;
}


void RRootBulkDS::SetNSlots(unsigned int nSlots)
{
   R__ASSERT(0U == fNSlots && "Setting the number of slots even if the number of slots is different from zero.");

   //printf("Setting parallelism to %u\n", nSlots);
   fNSlots = nSlots;
   fChains.resize(fNSlots);
   fBufferMgrs.resize(fNSlots);

   for (auto slot : ROOT::TSeqU(fNSlots)) {
      InitialiseSlot(slot);
   }
}

std::string RRootBulkDS::GetDataSourceType()
{
   return "RootBulk";
}

RDataFrame MakeRootBulkDataFrame(std::string_view treeName, std::string_view fileNameGlob)
{
   ROOT::RDataFrame tdf(std::make_unique<RRootBulkDS>(treeName, fileNameGlob));
   return tdf;
}

} // ns Experimental

} // ns RDF

} // ns ROOT

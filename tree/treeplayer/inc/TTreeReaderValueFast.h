// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-08-02

/*************************************************************************
 * Copyright (C) 1995-2013, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderValueFast
#define ROOT_TTreeReaderValueFast


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReaderValueFast                                                   //
//                                                                        //
// A simple interface for reading data from trees or chains.              //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TBufferFile.h"
#include "TTreeReaderFast.h"

#include <type_traits>

class TBranch;

namespace ROOT {
namespace Internal {

/* All the common code shared by the fast reader templates.
 */
class TTreeReaderValueFastBase {
   public:
      TTreeReaderValueFastBase(const TTreeReaderValueFastBase&) = delete;

      TTreeReaderValueBase::ESetupStatus GetSetupStatus() const { return fSetupStatus; }
      virtual TTreeReaderValueBase::EReadStatus GetReadStatus() const { return fReadStatus; }

      //////////////////////////////////////////////////////////////////////////////
      /// Construct a tree value reader and register it with the reader object.
      TTreeReaderValueFastBase(TTreeReaderFast* reader, const std::string &branchName) :
         fBranchName(branchName),
         fLeafName(branchName), // TODO: only support single-leaf branches for now.
         fTreeReader(reader),
         fBuffer(TBuffer::kWrite, 32*1024),
         fEvtIndex(reader->GetIndexRef())
      {
         if (fTreeReader) fTreeReader->RegisterValueReader(this);
      }

      Int_t GetEvents(Long64_t eventNum) {
          //printf("Getting events starting at %lld.  Current remaining is %d events with base %lld.\n", eventNum, fRemaining, fEventBase);
          if (fEventBase >= 0 && (fRemaining + fEventBase > eventNum)) {
             Int_t adjust = (eventNum - fEventBase);
             if (R__unlikely(Adjust(adjust) < 0)) {
                return -1;
             }
             fRemaining -= adjust;
          }
          fRemaining = fBranch->GetEntriesSerialized(eventNum, fBuffer);
          if (R__unlikely(fRemaining < 0)) {
             fReadStatus = TTreeReaderValueBase::kReadError;
             return -1;
          }
          fEventBase = eventNum;
          //printf("After getting events, the base is %lld with %d remaining.\n", fEventBase, fRemaining);
          fReadStatus = TTreeReaderValueBase::kReadSuccess;
          return fRemaining;
      }

      virtual const char *GetTypeName() {return "{UNDETERMINED}";}


   protected:

      // Adjust the current buffer offset forward N events.
      virtual Int_t Adjust(Int_t eventCount) {
         Int_t bufOffset = fBuffer.Length();
         fBuffer.SetBufferOffset(bufOffset + eventCount*GetSize());
         return 0;
      }
      virtual Int_t GetSize() = 0;

      void MarkTreeReaderUnavailable() {
         fTreeReader = nullptr;
      }

      // Create the linkage between the TTreeReader's current tree and this ReaderValue
      // object.  After CreateProxy() is invoked, if fSetupStatus doesn't indicate an
      // error, then we are pointing toward a valid TLeaf in the current tree
      void CreateProxy();

      virtual ~TTreeReaderValueFastBase();

      // Returns the name of the branch type; will be used when the TBranch version to
      // detect between the compile-time and runtime type names.
      virtual const char *BranchTypeName() = 0;

      std::string  fBranchName;          // Name of the branch we should read from.
      std::string  fLeafName;            // The branch's leaf we should read from.  NOTE: currently only support single-leaf branches.
      TTreeReaderFast *fTreeReader{nullptr}; // Reader we belong to
      TBranch *    fBranch{nullptr};     // Actual branch object we are reading.
      TLeaf *      fLeaf{nullptr};       // Actual leaf we are reading.
      TBufferFile  fBuffer;              // Buffer object holding the current events.
      Int_t        fRemaining{0};        // Number of events remaining in the buffer.
      Int_t       &fEvtIndex;            // Current event index.
      Long64_t     fLastChainOffset{-1}; // Current chain in the TTree we are pointed at.
      Long64_t     fEventBase{-1};       // Event number of the current buffer position.

      TTreeReaderValueBase::ESetupStatus fSetupStatus{TTreeReaderValueBase::kSetupNotSetup}; // setup status of this data access
      TTreeReaderValueBase::EReadStatus  fReadStatus{TTreeReaderValueBase::kReadNothingYet}; // read status of this data access

      friend class ::TTreeReaderFast;
};

}  // Internal
}  // ROOT

template <typename T>
class TTreeReaderValueFast final : public ROOT::Internal::TTreeReaderValueFastBase {

   public:
       TTreeReaderValueFast(TTreeReaderFast* reader, const std::string &branchname) : ROOT::Internal::TTreeReaderValueFastBase(reader, branchname) {}

      T* Get() {
         return Deserialize(reinterpret_cast<char *>(reinterpret_cast<T*>(fBuffer.GetCurrent()) + fEvtIndex));
      }
      T* operator->() { return Get(); }
      T& operator*() { return *Get(); }

   protected:
      T* Deserialize(char *) {return nullptr;}

      virtual const char *GetTypeName() override {return "{INCOMPLETE}";}
      virtual Int_t GetSize() override {return sizeof(T);}
};

template <>
class TTreeReaderValueFast<float> final : public ROOT::Internal::TTreeReaderValueFastBase {

   public:

      TTreeReaderValueFast(TTreeReaderFast& tr, const std::string &branchname) :
            TTreeReaderValueFastBase(&tr, branchname) {}

      // TODO: why isn't template specialization working here?
      float* Get() {
         //printf("Float: Attempting to deserialize buffer %p from index %d.\n", fBuffer.GetCurrent(), fEvtIndex);
         return Deserialize(reinterpret_cast<char *>(reinterpret_cast<float*>(fBuffer.GetCurrent()) + fEvtIndex));
      }
      float* operator->() { return Get(); }
      float& operator*() { return *Get(); }

   protected:
      virtual const char *GetTypeName() override {return "float";}
      virtual const char *BranchTypeName() override {return "float";}
      virtual Int_t GetSize() override {return sizeof(float);}
      float * Deserialize(char *input) {frombuf(input, &fTmp); return &fTmp;}

      float fTmp;
};

#endif // ROOT_TTreeReaderValueFast

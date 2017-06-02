// @(#)root/tree:$Id$
// Author: Axel Naumann, 2010-08-02

/*************************************************************************
 * Copyright (C) 1995-2017, Rene Brun and Fons Rademakers.               *
 * All rights reserved.                                                  *
 *                                                                       *
 * For the licensing terms see $ROOTSYS/LICENSE.                         *
 * For the list of contributors see $ROOTSYS/README/CREDITS.             *
 *************************************************************************/

#ifndef ROOT_TTreeReaderFast
#define ROOT_TTreeReaderFast


////////////////////////////////////////////////////////////////////////////
//                                                                        //
// TTreeReader                                                            //
//                                                                        //
// A simple interface for reading trees or chains.                        //
//                                                                        //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TTreeReader.h"

#include <deque>

namespace ROOT {
namespace Internal {
class TTreeReaderValueFast;
class TTreeReaderValueFastBase;
}
}

class TTreeReaderFast: public TObject {
public:

   // A simple iterator based on TTreeReader::Iterator_t; allows use of the
   // TTreeReaderFast.
   //
   // NOTE that an increment may invalidate previous copies of the iterator.
   class Iterator_t:
      public std::iterator<std::input_iterator_tag, const Long64_t, Long64_t> {
   private:
      Int_t fIdx{0}; ///< Current offset inside this cluster.
      Int_t fCount{0}; ///< Number of entries inside this cluster.
      Long64_t fEntry{-1}; ///< Entry number of the tree referenced by this iterator; -1 is invalid.
      TTreeReaderFast* fReader{nullptr}; ///< The reader we select the entries on.

      /// Whether the iterator points to a valid entry.
      bool IsValid() const { return fEntry >= 0; }

   public:
      /// Default-initialize the iterator as "past the end".
      Iterator_t() {}

      /// Initialize the iterator with the reader it steers and a
      /// tree entry number; -1 is invalid.
      Iterator_t(TTreeReaderFast& reader, Long64_t first, Long64_t count):
         fCount(count), fEntry(first), fReader(&reader) {}

      /// Compare two iterators for equality.
      bool operator==(const Iterator_t& lhs) const {
         // From C++14: value initialized (past-end) it compare equal.
         if (R__unlikely(!IsValid() && !lhs.IsValid())) {return true;}
         return R__unlikely(fEntry == lhs.fEntry && fReader == lhs.fReader);
      }

      /// Compare two iterators for inequality.
      bool operator!=(const Iterator_t& lhs) const {
         return !(*this == lhs);
      }

      /// Increment the iterator (postfix i++).
      Iterator_t operator++(int) {
         Iterator_t ret = *this;
         this->operator++();
         return ret;
      }

      /// Increment the iterator (prefix ++i).
      Iterator_t& operator++() {
         fIdx++;
         if (R__unlikely(fIdx == fCount)) {
             fEntry += fCount;
             fIdx = 0;
             fCount = fReader->GetNextRange(fCount);
             if (R__unlikely(!fCount)) {
                 fEntry = -1;
             }
         }
         return *this;
      }

      /// Set the entry number in the reader and return it.
      Long64_t operator*() {
         return fEntry + fIdx;
      }

      Long64_t operator*() const {
         return **const_cast<Iterator_t*>(this);
      }
   };

   typedef Iterator_t iterator;

   TTreeReaderFast():
      fTree(nullptr),
      fEntryStatus(TTreeReader::kEntryNoTree),
      fLastEntry(-1)
   {}

   TTreeReaderFast(TTree* tree);
   TTreeReaderFast(const char* keyname, TDirectory* dir = NULL );

   ~TTreeReaderFast();

   TTreeReader::EEntryStatus SetEntriesRange(Long64_t first, Long64_t last);

   TTreeReader::EEntryStatus GetEntryStatus() const { return fEntryStatus; }

   TTree* GetTree() const { return fTree; }

   /// Return an iterator to the 0th TTree entry.
   Iterator_t begin() {
      return Iterator_t(*this, 0, this->GetNextRange(0));
   }
   Iterator_t end() const { return Iterator_t(); }

protected:

   // Returns a reference to the current event index in the various value buffers.
   Int_t &GetIndexRef() {return fEvtIndex;}

   void RegisterValueReader(ROOT::Internal::TTreeReaderValueFastBase* reader);
   void DeregisterValueReader(ROOT::Internal::TTreeReaderValueFastBase* reader);

private:

   Int_t GetNextRange(Int_t);
   void Initialize();

   TTree* fTree{nullptr}; ///< tree that's read
   TDirectory* fDirectory{nullptr}; ///< directory (or current file for chains)
   ROOT::Internal::TBranchProxyDirector* fDirector{nullptr}; ///< proxying director, owned
   TTreeReader::EEntryStatus fEntryStatus{TTreeReader::kEntryNoTree}; ///< status of most recent read request
   std::deque<ROOT::Internal::TTreeReaderValueFastBase*> fValues; ///< readers that use our director

   Int_t    fEvtIndex{-1};
   Long64_t fBaseEvent{-1};
   Long64_t fLastEntry{-1};

   friend class ROOT::Internal::TTreeReaderValueFastBase;
   friend class ROOT::Internal::TTreeReaderArrayBase;

   ClassDef(TTreeReaderFast, 0); // A simple interface to read trees
};

#endif // defined TTreeReader

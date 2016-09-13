#include "TMTrackTrigger/TMTrackFinder/interface/L1track2D.h"

// Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.

L1track2D L1track2D::mergeTracks(const L1track2D B) const {

  vector<const Stub*> aStubs=this->getStubs(), bStubs=B.getStubs();
  vector<const Stub*> mStubs;
  unsigned int aL=aStubs.size(), bL=bStubs.size();
  unsigned int aP=0, bP=0;

  // This loop relies on the stubs being ordered in the same way inside the two tracks.
  while ((aP<aL) && (bP<bL))
    {  unsigned int aS=aStubs[aP]->index(), bS=bStubs[bP]->index();
      if (aS == bS)
	{ mStubs.push_back(aStubs[aP]);
	  ++aP; ++bP;
	}
      else if (aS < bS)
	{ mStubs.push_back(aStubs[aP]);
	  ++aP;
	}
      else // bS < aS
	{ mStubs.push_back(bStubs[bP]);
	  ++bP;
	}
    }
  for (; aP<aL; aP++) { mStubs.push_back(aStubs[aP]);} //if any left in A -- at most one of these 2 loops
  for (; bP<bL; bP++) { mStubs.push_back(bStubs[bP]);} //if any left in B -- should be executed

  // N.B. This defines the HT cell location as that of the first track, meaning that the merged tracks depends
  // on which track is first and which is second. This will make it hard to get identical results from hardware 
  // & software.
  return L1track2D(settings_, mStubs, this->getCellLocation(), this->getHelix2D(), this->isRphiTrk());
}

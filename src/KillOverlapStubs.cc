#include "TMTrackTrigger/TMTrackFinder/interface/KillOverlapStubs.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"

#include <vector>

// Considering distinct pairs of stubs on same layer but different modules, find those which correspond to the same track.
std::vector<const Stub*> KillOverlapStubs::pairFinder() const
{
  set<const Stub*> dStubs;

  // consider all distinct pairs of stubs
  for ( auto i1 = vStubs_.begin(); i1 != vStubs_.end(); ++i1 )
    for ( auto i2 = i1+1; i2 != vStubs_.end(); ++i2) {
      const Stub* s1 = *i1; const Stub* s2 = *i2;

      // check if stubs are on neighbouring modules
      if ( ! neighb_modules(s1, s2) ) continue;

      // cuts in r-z and r-phi planes
      double ptOverQ = ( !s1->barrel() && !s1->psModule() )
        ? settings_->invPtToDphi() * (s1->z() - s2->z()) * s1->r() / s1->z() / (s2->phi() - s1->phi()) // for endcap 2S modules
        : settings_->invPtToDphi() * (s1->r() - s2->r()) / (s2->phi() - s1->phi());                    // for all other modules
      double z0 = s2->z() - s2->r() * (s2->z() - s1->z()) / (s2->r() - s1->r());
      if ( fabs(z0) > z0_cut_ ) continue;
      if ( fabs(ptOverQ) < pt_cut_ )              continue;

      // check if qOverPt() of both stubs matches the above
      if ( fabs(1/ptOverQ - s1->qOverPt()) > s1->qOverPtres() ) continue;
      if ( fabs(1/ptOverQ - s2->qOverPt()) > s2->qOverPtres() ) continue;

      // then the pair probably corresponds to the same track
      dStubs.insert(s1);
    }

  // if stub is not a first element of a pair, then add it to vStubs_filtered
  std::vector<const Stub*> vStubs_filtered;
  for (const Stub* s : vStubs_)
    if ( find(dStubs.begin(), dStubs.end(), s ) == dStubs.end() )
      vStubs_filtered.push_back(s);

  return vStubs_filtered;
}

// Check if two stubs are on neighbouring modules
bool KillOverlapStubs::neighb_modules(const Stub* s1, const Stub* s2) const
{
  // check if stubs are on the same layer
  if ( s1->layerId() != s2->layerId() )
    return false;

  // same-layer modules must be either both in barrel, or both in endcap
  if ( s1->barrel() != s2->barrel() )
    return false; // presumably never reached

  // check if stubs are on different modules
  if ( s1->idDet() == s2->idDet() )
    return false;

  // neighbouring modules are about 10 cm apart
  const float delta = 12.0; // 12.0 = about 10cm

  // neighbouring in phi: r*phi differs by about 10 cm
  const float ctrR1   = 0.5*(s1->minR()+s1->maxR());
  const float ctrR2   = 0.5*(s2->minR()+s2->maxR());
  const float ctrPhi1 = 0.5*(s1->minPhi()+s1->maxPhi());
  const float ctrPhi2 = 0.5*(s2->minPhi()+s2->maxPhi());
  const float rPhi1   = ctrR1 * ctrPhi1;
  const float rPhi2   = ctrR2 * ctrPhi1;
  if ( fabs(rPhi1-rPhi2) > delta )
    return false;

  // neighbouring in r-z in barrel: z differs by about 10 cm
  const float ctrZ1   = 0.5*(s1->minZ()+s1->maxZ());
  const float ctrZ2   = 0.5*(s2->minZ()+s2->maxZ());
  if ( s1->barrel() ) {
    if ( fabs( ctrZ1 - ctrZ2 ) > delta )
      return false;
  // neighbouring in r-z in endcap: r differs by about 10 cm
  } else {
    if ( fabs( ctrR1 - ctrR2 ) > delta )
      return false;
  }

  // neighbouring
  return true;
}


std::vector<const Stub*> KillOverlapStubs::truePairFinder() const
{
  set<const Stub*> dStubs;
  for ( auto i1 = vStubs_.begin(); i1 != vStubs_.end(); ++i1 )
    for ( auto i2 = i1+1; i2 != vStubs_.end(); ++i2) {
      // require same layer
      if ( (*i1)->layerId() != (*i2)->layerId() ) continue;
      // ignore stubs from the same module
      if ( (*i1)->idDet() == (*i2)->idDet() )     continue;
      // require at least one TP in common
      if ( commonTP(*i1, *i2) == nullptr )        continue;
      // require TP.pt >= pt_cut_
      if ( commonTP(*i1, *i2)->pt() < pt_cut_ )    continue;
      // stubs in pair (*i1, *i2) are the ones the algorithm needs to find
      dStubs.insert(*i1);
    }

  // if stub is not a first element of a pair, then add it to vStubs_filtered
  std::vector<const Stub*> vStubs_filtered;
  for (const Stub* s : vStubs_)
    if ( find(dStubs.begin(), dStubs.end(), s ) == dStubs.end() )
      vStubs_filtered.push_back(s);

  return vStubs_filtered;
}

// If the two stubs share at least one TP in common, we consider them as belonging to the same track.
const TP* KillOverlapStubs::commonTP(const Stub* s1, const Stub* s2) const
{
  // both stubs must have at least one TP associated with them
  if ( ! (s1->genuine() && s2->genuine()) ) return nullptr;

  // find at least one TP in common between the two stubs, and return it if found
  for ( const TP* tp1 : s1->assocTPs() ) {
    set<const TP*> tp2s = s2->assocTPs();
    if ( tp2s.find(tp1) != tp2s.end() )
      return tp1;
  }

  // stubs have no TPs in common
  return nullptr;
}


#ifndef __SECTOR_H__
#define __SECTOR_H__

#include <vector>
#include <unordered_map>

class Settings;
class Stub;
class TP;

using namespace std;

class Sector {

public:
  
  Sector() {}
  ~Sector() {}

  // Initialization.
  void init(const Settings* settings, unsigned int iPhiSec, unsigned int iEtaSec);

  // Check if stub within the eta and/or phi boundaries of this sector.
  bool inside   ( const Stub* stub ) const {return (this->insideEta(stub) && this->insidePhi(stub));}
  bool insideEta( const Stub* stub ) const;
  bool insidePhi( const Stub* stub ) const;

  // Check if stub is within subsectors in eta that sector may be divided into.
  vector<bool> insideEtaSubSecs( const Stub* stub) const;

  float phiCentre() const { return phiCentre_; } // Return phi of centre of this sector.
  float etaMin()    const { return etaMin_; } // Eta range covered by this sector.
  float etaMax()    const { return etaMax_; } // Eta range covered by this sector.

  // For performance studies, note which stubs on given tracking particle are inside the sector.
  // Returns two booleans for each stub, indicating if they are in phi & eta sectors respectively.
  // You can AND them together to check if stub is in (eta,phi) sector.
  unordered_map<const Stub*, pair<bool, bool>> stubsInside ( const TP& tp) const;

  // Count number of stubs in given tracking particle which are inside this (phi,eta) sector;
  // or inside it if only the eta cuts are applied; or inside it if only the phi cuts are applied.
  // The results are returned as the 3 last arguments of the function.
  void numStubsInside( const TP& tp, 
                       unsigned int& nStubsInsideEtaPhi, unsigned int& nStubsInsideEta, 
		       unsigned int& nStubsInsidePhi) const ;

private: 

  // Check if stub is within eta sector or subsector that is delimated by specified zTrk range.
  bool insideEtaRange( const Stub* stub, float zRangeMin, float zRangeMax) const;

private:

  const Settings* settings_;

  float  beamWindowZ_;
  float  trackerOuterRadius_;
  float  trackerInnerRadius_;
  float  trackerHalfLength_;
  bool   handleStripsPhiSec_;
  bool   handleStripsEtaSec_;

  // Define eta region.
  float  etaMin_; // Range in eta covered by this sector.
  float  etaMax_;
  float  chosenRofZ_; // Use z of track at radius="chosenRofZ" to define eta sectors.
  float  rOuterMax_; // Larger eta boundary point (r,z)
  float  zOuterMax_;
  float  rOuterMin_; // Smaller eta boundary point (r,z)
  float  zOuterMin_;

  // Define phi sector.
  float  phiCentre_; // phi of centre of sector.
  float  sectorHalfWidth_; // sector half-width excluding overlaps.
  float  chosenRofPhi_; // Use phi of track at radius="chosenRofPhi" to define phi sectors.
  bool   useStubPhi_; // Require stub phi to be consistent with track of Pt > HTArraySpec.HoughMinPt that crosses HT phi axis?
  float  minPt_; // Min Pt covered by HT array.
  bool   useStubPhiTrk_;  // Require stub phi0 (or phi65 etc.) as estimated from stub bend, to lie within HT phi axis, allowing tolerance specified below?
  float  assumedPhiTrkRes_; // Tolerance in stub phi0 (or phi65) assumed to be this fraction of phi sector width. (N.B. If > 0.5, then stubs can be shared by more than 2 phi sectors). 
  bool   calcPhiTrkRes_; // If true, tolerance in stub phi0 (or phi65 etc.) will be reduced below AssumedPhiTrkRes if stub bend resolution specified in HTFilling.BendResolution suggests it is safe to do so.

  // Possible subsectors in eta within each sector.
  unsigned int numSubSecsEta_;
  vector<float> zOuterMinSub_;
  vector<float> zOuterMaxSub_;
};
#endif


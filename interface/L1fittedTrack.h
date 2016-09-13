#ifndef __L1fittedTrack_H__
#define __L1fittedTrack_H__

#include "FWCore/Utilities/interface/Exception.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1trackBase.h"
#include "TMTrackTrigger/TMTrackFinder/interface/L1track3D.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Settings.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Utility.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <vector>
#include <utility>

using namespace std;

//=== This represents a fitted L1 track candidate found in 3 dimensions.
//=== It gives access to the fitted helix parameters & chi2 etc.
//=== It also calculates & gives access to associated truth particle (Tracking Particle) if any.
//=== It also gives access to the 3D hough-transform track candidate (L1track3D) on which the fit was run.

class L1fittedTrack : public L1trackBase {

public:

  // Store a new fitted track, specifying the input Hough transform track, the stubs used for the fit,
  // the fitted helix parameters & chi2,
  // and the number of helix parameters being fitted (=5 if d0 is fitted, or =4 if d0 is not fitted).
  // Also specify phi sector and eta region used by track-finding code that this track was in.
  // And if track fit declared this to be a valid track (enough stubs left on track after fit etc.).
  L1fittedTrack(const Settings* settings, const L1track3D& l1track3D, const vector<const Stub*>& stubs,
                float qOverPt, float d0, float phi0, float z0, float tanLambda, 
                float chi2, unsigned int nHelixParam,
                unsigned int iPhiSec, unsigned int iEtaReg, bool accepted = true) :
    L1trackBase(),
    settings_(settings),
    l1track3D_(l1track3D), stubs_(stubs),
    qOverPt_(qOverPt), d0_(d0), phi0_(phi0), z0_(z0), tanLambda_(tanLambda), 
    chi2_(chi2), nHelixParam_(nHelixParam),
    iPhiSec_(iPhiSec), iEtaReg_(iEtaReg), accepted_(accepted)
  {
    nLayers_   = Utility::countLayers(settings, stubs); // Count tracker layers these stubs are in
    matchedTP_ = Utility::matchingTP(settings, stubs, nMatchedLayers_, matchedStubs_); // Find associated truth particle & calculate info about match.
  }

  ~L1fittedTrack() {}

  //--- Get the 3D Hough transform track candididate corresponding to the fitted track,
  //--- Provide direct access to some of the info it contains.

  // Get track candidate from HT (before fit).
  const L1track3D&            getL1track3D()          const  {return l1track3D_;}

  // Get stubs on fitted track (can differ from those on HT track if track fit kicked out stubs with bad residuals)
  const vector<const Stub*>&  getStubs()              const  {return stubs_;}  
  // Get number of stubs on fitted track.
  unsigned int                getNumStubs()           const  {return stubs_.size();}
  // Get number of tracker layers these stubs are in.
  unsigned int                getNumLayers()          const  {return nLayers_;}
  // Get number of stubs deleted from track candidate by fitter (because they had large residuals)
  unsigned int             getNumKilledStubs()        const  {return l1track3D_.getNumStubs() - this->getNumStubs();}

  // Get cell locations of the track candidate corresponding to the fitted track in the r-phi and r-z Hough transforms in units of bin number.
  pair<unsigned int, unsigned int>  getCellLocationRphi() const  {return l1track3D_.getCellLocationRphi();}
  pair<unsigned int, unsigned int>  getCellLocationRz()   const  {return l1track3D_.getCellLocationRz();}

  //--- Get information about its association (if any) to a truth Tracking Particle.
  //--- Can differ from that of corresponding HT track, if track fit kicked out stubs with bad residuals.

  // Get best matching tracking particle (=nullptr if none).
  const TP*                   getMatchedTP()          const  {return matchedTP_;}
  // Get the matched stubs with this Tracking Particle
  const vector<const Stub*>&  getMatchedStubs()       const  {return matchedStubs_;}
  // Get number of matched stubs with this Tracking Particle
  unsigned int                getNumMatchedStubs()    const  {return matchedStubs_.size();}
  // Get number of tracker layers with matched stubs with this Tracking Particle 
  unsigned int                getNumMatchedLayers()   const  {return nMatchedLayers_;}
  // Get purity of stubs on track (i.e. fraction matching best Tracking Particle)
  float                       getPurity()             const   {return getNumMatchedStubs()/float(getNumStubs());}
  // Get number of stubs matched to correct TP that were deleted from track candidate by fitter.
  unsigned int            getNumKilledMatchedStubs()  const  {
    unsigned int nStubCount = l1track3D_.getNumMatchedStubs();
    if (nStubCount > 0) { // Original HT track candidate did match a truth particle
      const TP* tp = l1track3D_.getMatchedTP(); 
      for (const Stub* s : stubs_) {
	set<const TP*> assTPs = s->assocTPs();
        if (assTPs.find(tp) != assTPs.end()) nStubCount--; // We found a stub matched to original truth particle that survived fit.
      }
    }
    return nStubCount;
  }

  //--- Get the fitted track helix parameters.

  float   qOverPt()      const  {return qOverPt_;}
  float   charge()       const  {return (qOverPt_ > 0  ?  1  :  -1);} 
  float   invPt()        const  {return fabs(qOverPt_);}
  float   pt()           const  {return 1./(1.0e-6 + this->invPt());} // includes protection against 1/pt = 0.
  float   d0()           const  {return d0_;}
  float   phi0()         const  {return phi0_;}
  float   z0()           const  {return z0_;}
  float   tanLambda()    const  {return tanLambda_;}
  float   theta()        const  {return atan2(1., tanLambda_);} // Use atan2 to ensure 0 < theta < pi.
  float   eta()          const  {return -log(tan(0.5*this->theta()));}

  // Get the number of helix parameters being fitted (=5 if d0 is fitted or =4 if d0 is not fitted).
  float   nHelixParam()  const  {return nHelixParam_;}
  // Get degrees of freedom.
  unsigned int numDOF()  const  {return 2*this->getNumStubs() - nHelixParam_;}
  // Get the fit chi2 and chi2/DOF
  float   chi2()         const  {return chi2_;}
  float   chi2dof()      const  {return chi2_/this->numDOF();}

  // Comparitor for sorting tracks by q/Pt using std::sort().
  static bool qOverPtSortPredicate(const L1fittedTrack& t1, const L1fittedTrack t2) { return t1.getCellLocationRphi().first < t2.getCellLocationRphi().first; }

  //--- Get phi sector and eta region used by track finding code that this track is in.
  unsigned int iPhiSec() const  {return iPhiSec_;}
  unsigned int iEtaReg() const  {return iEtaReg_;}

  //--- Get whether the track has been rejected or accepted by the fit
  bool accepted()        const  {return accepted_;}

  //--- Function for merging two tracks into a single track, used by by KillDupTracks.h for duplicate track removal.
  L1fittedTrack mergeTracks(const L1fittedTrack B) const {throw cms::Exception("L1fittedTrack ERROR: function mergeTracks(L1fittedTrack) has not yet been implemented!");}

private:

  //--- Configuration parameters
  const Settings*                    settings_; 

  //--- The 3D hough-transform track candidate which was fitted.
  //  const L1track3D&             l1track3D_;
  L1track3D             l1track3D_;

  //--- The stubs on the fitted track (can differ from those on HT track if fit kicked off stubs with bad residuals)
  vector<const Stub*>   stubs_;
  unsigned int          nLayers_;

  //--- The fitted helix parameters and fit chi-squared.
  float qOverPt_;
  float d0_;
  float phi0_;
  float z0_;
  float tanLambda_;
  float chi2_;

  //--- The number of helix parameters being fitted (=5 if d0 is fitted or =4 if d0 is not fitted).
  unsigned int nHelixParam_;

  //--- Phi sector and eta region used track finding code that this track was in.
  unsigned int iPhiSec_;
  unsigned int iEtaReg_; 

  //--- Information about its association (if any) to a truth Tracking Particle.
  const TP*             matchedTP_;
  vector<const Stub*>   matchedStubs_;
  unsigned int          nMatchedLayers_;

  //--- Has the track fit declared this to be a valid track?
  bool accepted_;
};
#endif

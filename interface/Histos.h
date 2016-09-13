#ifndef __HISTOS_H__
#define __HISTOS_H__

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TMTrackTrigger/TMTrackFinder/interface/Settings.h>
#include <TMTrackTrigger/TMTrackFinder/interface/Stub.h>

#include "boost/numeric/ublas/matrix.hpp"
using  boost::numeric::ublas::matrix;

#include <vector>
#include <map>
#include <string>

class InputData;
class TP;
class Sector;
class HTpair;
class L1fittedTrack;
class L1fittedTrk4and5;
class KillOverlapStubs;
class TH1F;
class TH2F;
class TProfile;
class TGraphAsymmErrors;

class Histos {

public:
  // Store cfg parameters.
  Histos(const Settings* settings) : settings_(settings), numPerfRecoTPforAlg_(0) {}

  ~Histos(){}

  // Book all histograms
  void book();

  // Fill histograms with stubs and tracking particles from input data.
  void fillInputData(const InputData& inputData);
  // Fill histograms relating the stub-pairs algorithm
  void fillStubPairs(const vector<const Stub*>& vStubs);
  // Fill histograms that check if choice of (eta,phi) sectors is good.
  void fillEtaPhiSectors(const InputData& inputData, const matrix<Sector>& mSectors);
  // Fill histograms checking filling of r-phi HT array.
  void fillRphiHT(const matrix<HTpair>& mHtPairs);
  // Book histograms about r-z track filters (or other filters applied after r-phi HT array).
  void fillRZfilters(const matrix<HTpair>& mHtPairs);
  // Fill histograms studying track candidates found by Hough Transform.
  void fillTrackCands(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs);
  // Fill histograms studying freak, events with too many stubs..
  void fillStudyBusyEvents(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs);
  // Fill histograms relating to track fitting performance.
  void fillTrackFitting(const InputData& inputData, const vector<std::pair<std::string,L1fittedTrack>>& fittedTracks, float chi2dofCutPlots);

  void endJobAnalysis();

private:

  // Book histograms for specific topics.
  void bookInputData();
  void bookStubPairs();
  void bookEtaPhiSectors();
  void bookRphiHT();
  void bookRZfilters();
  void bookTrackCands();
  void bookStudyBusyEvents();
  void bookTrackFitting();

  // Produce plots of tracking efficiency prior to track fit (run at end of job).
  void plotTrackEfficiency();
  // Produce plots of tracking efficiency after track fit (run at end of job).
  void plotTrackEffAfterFit(string fitName);

  // Understand why not all tracking particles were reconstructed.
  // Returns list of tracking particles that were not reconstructed and an integer indicating why.
  // Only considers TP used for algorithmic efficiency measurement.
  map<const TP*, string> diagnoseTracking(const InputData& inputData, const matrix<Sector>& mSectors, const matrix<HTpair>& mHtPairs) const;

 private:

  const Settings *settings_; // Configuration parameters.

  edm::Service<TFileService> fs_;

  // Histograms of input data.
  TProfile* profNumStubs_;
  TH1F* hisStubsVsEta_;
  TH1F* hisStubsVsR_;
  TProfile* profNumTPs_;
  TH1F* hisNumStubsPerTP_;
  TProfile* hisStubKillFE_;
  TProfile* hisStubIneffiVsInvPt_;
  TProfile* hisStubIneffiVsEta_;
  TProfile* hisStubKillDataCorr_;
  TH1F* hisPtStub_;
  TH1F* hisPtResStub_;
  TH1F* hisBendFilterPower_;
  TH1F* hisDelPhiStub_;
  TH1F* hisDelPhiResStub_;
  TH1F* hisBendStub_;
  TH1F* hisBendResStub_;
  TH1F* hisNumMergedBend_;
  TH2F* hisBendVsLayerOrRing_;
  TH2F* hisBendFEVsLayerOrRing_;
  TH1F* hisPhiStubVsPhiTP_;
  TH1F* hisPhiStubVsPhi0TP_;
  TH1F* hisPhi0StubVsPhi0TP_;
  TH1F* hisPhi0StubVsPhi0TPres_;
  TH1F* hisPhiStubVsPhi65TP_;
  TH1F* hisPhi65StubVsPhi65TP_;
  TH1F* hisPhi65StubVsPhi65TPres_;
  TH1F* hisPitchOverSep_;
  TH1F* hisRhoParameter_;
  TH1F* hisFracStubsSharingClus0_;
  TH1F* hisFracStubsSharingClus1_;

  // -------------------------------------------------------
  // Histograms to do with the pair finding algorithm
  // -------------------------------------------------------

  // for testing the formula for z0 and q/pt from the two stubs
  TH1F* hisPtPairMinusTruth_;
  TH1F* hisZ0PairMinusTruth_;
  TH2F* hisZ0PairVsTruth_;
  TH2F* hisPtPairVsTruth_;
  TH2F* hisZ0PairVsTruth_barrel_;
  TH2F* hisZ0PairVsTruth_endcap_;
  TH2F* hisZ0PairVsTruth_barrel_PS_;
  TH2F* hisZ0PairVsTruth_barrel_2S_;
  TH2F* hisZ0PairVsTruth_endcap_PS_;
  TH2F* hisZ0PairVsTruth_endcap_2S_;
  TH2F* hisZ0PairVsTruth_PS_;
  TH2F* hisZ0PairVsTruth_2S_;
  TH2F* hisPtPairVsTruth_barrel_;
  TH2F* hisPtPairVsTruth_endcap_;
  TH2F* hisPtPairVsTruth_PS_;
  TH2F* hisPtPairVsTruth_2S_;
  TH2F* hisPtPairVsTruth_barrel_PS_;
  TH2F* hisPtPairVsTruth_barrel_2S_;
  TH2F* hisPtPairVsTruth_endcap_PS_;
  TH2F* hisPtPairVsTruth_endcap_2S_;
  TH1F* hisZ0PairMinusTruth_barrel_PS_;
  TH1F* hisZ0PairMinusTruth_barrel_PS_abs_;

  // Histograms for finding module dimensions
  TH2F* his_dim1_;
  TH2F* his_dim2_;

  // stats on my pair-finding algorithm;
  TH1F* his_StubsInFound_Pt      ;
  TH1F* his_AllStubs_Pt          ;
  TH1F* his_AllStubs_AbsPt       ;
  TH1F* his_TrueFoundPairs_Pt    ;
  TH1F* his_TrueFoundStubs_Pt    ;
  TH1F* his_TrueFoundStubs_AbsPt ;
  TH1F* his_WrongFoundPairs_Pt   ;
  TH1F* his_AllTruePairs_Pt      ;
  TH1F* his_AllTrueStubs_Pt      ;
  TH1F* his_AllTrueStubs_AbsPt   ;
  TH1F* his_AllFoundPairs_Pt     ;
  TH1F* his_AllFoundStubs_Pt     ;
  TH1F* his_AllFoundStubs_AbsPt  ;
  TH1F* his_WrongStubs_Pt        ;
  TH1F* his_WrongStubs_AbsPt     ;
  TH1F* his_StubsInFound_Eta     ;
  TH1F* his_AllStubs_Eta         ;
  TH1F* his_TrueFoundPairs_Eta   ;
  TH1F* his_TrueFoundStubs_Eta   ;
  TH1F* his_WrongFoundPairs_Eta  ;
  TH1F* his_AllTruePairs_Eta     ;
  TH1F* his_AllTrueStubs_Eta     ;
  TH1F* his_AllFoundPairs_Eta    ;
  TH1F* his_AllFoundStubs_Eta    ;
  TH1F* his_WrongStubs_Eta       ;
  TH2F* his_StubsInFound_Loc     ;
  TH2F* his_AllStubs_Loc         ;
  TH2F* his_TrueFoundPairs_Loc   ;
  TH2F* his_TrueFoundStubs_Loc   ;
  TH2F* his_WrongFoundPairs_Loc  ;
  TH2F* his_AllTruePairs_Loc     ;
  TH2F* his_AllTrueStubs_Loc     ;
  TH2F* his_AllFoundPairs_Loc    ;
  TH2F* his_AllFoundStubs_Loc    ;
  TH2F* his_WrongStubs_Loc       ;
  TH1F* his_AllTruePairs_dZ      ;
  TH1F* his_AllTruePairs_dZm     ;
  TH1F* his_AllTruePairs_dRPm    ;
  TH1F* his_NeighBarPS_dZm       ;
  TH1F* his_NeighBar2S_dZm       ;
  TH1F* his_NeighEndPS_dZm       ;
  TH1F* his_NeighEnd2S_dZm       ;
  TH1F* his_NeighBarPS_dRPm      ;
  TH1F* his_NeighBar2S_dRPm      ;
  TH1F* his_NeighEndPS_dRPm      ;
  TH1F* his_NeighEnd2S_dRPm      ;
  TH2F* his_RPhiPos_BarPS        ;
  TH2F* his_RPhiPos_Bar2S        ;
  TH2F* his_RPhiPos_EndPS        ;
  TH2F* his_RPhiPos_End2S        ;
  TH1F* his_fracKilled           ;
  TH1F* his_trueFracKilled       ;
  TH1F* his_pt_cut_StubsInFoundPairs      ;
  TH1F* his_pt_cut_StubsInTrueFoundPairs  ;
  TH1F* his_pt_cut_StubsInWrongFoundPairs ;
  TH1F* his_pt_cut_StubsInAllTruePairs    ;
  TH1F* his_pt_cut_AllStubs               ;
  TH1F* his_z0_cut_StubsInFoundPairs      ;
  TH1F* his_z0_cut_StubsInTrueFoundPairs  ;
  TH1F* his_z0_cut_StubsInWrongFoundPairs ;
  TH1F* his_z0_cut_StubsInAllTruePairs    ;
  TH1F* his_z0_cut_AllStubs               ;

  // stats on pairs of clusters
  TH1F* hisNumShares_;
  TH1F* hisStubsPerModule_;


  // Histograms checking that (eta,phi) sector definition is good.
  TH1F* hisFracStubsInSec_;
  TH1F* hisFracStubsInEtaSec_;
  TH1F* hisFracStubsInPhiSec_;
  TH1F* hisNumSecsPerStub_;
  TH1F* hisNumEtaSecsPerStub_;
  TH1F* hisNumPhiSecsPerStub_;
  TH1F* hisNumStubsPerSec_;
  TProfile* profNumStubsPerEtaSec_;
  TH2F* hisLayerIDvsEtaSec_;
  TH2F* hisLayerIDreducedvsEtaSec_;

  // Histograms checking filling of r-phi HT array.
  TH1F* hisIncStubsPerHT_;
  TH1F* hisExcStubsPerHT_;
  TH2F* hisNumStubsInCellVsEta_;
  TH1F* hisStubsOnRphiTracksPerHT_;

  // Histograms about r-z track filters (or other filters applied after r-phi HT array).
  TH1F* hisNumZtrkSeedCombinations_;
  TH1F* hisNumSeedCombinations_;
  TH1F* hisNumGoodSeedCombinations_;
  TH1F* hisCorrelationZTrk_;

  // Histograms studying track candidates found by Hough Transform.
  TProfile* profNumTrackCands_;
  TProfile* profNumTracksVsEta_;
  TH1F*     hisNumTracksVsQoverPt_;
  TProfile* profStubsOnTracks_;
  TProfile* profStubsOnTracksVsEta_;
  TH1F*     hisStubsOnTracksPerSect_;
  TH1F*     hisStubsPerTrack_;
  TH1F*     hisLayersPerTrack_;
  TH1F*     hisPSLayersPerTrack_;
  TProfile* profExcessStubsPerTrackVsPt_;
  TH1F*     hisFracMatchStubsOnTracks_;
  TH1F* hisDeltaPhiRtruePS_;
  TH1F* hisDeltaRorZtruePS_;
  TH1F* hisDeltaPhiRtrue2S_;
  TH1F* hisDeltaRorZtrue2S_;
  TH1F* hisDeltaPhiRfakePS_;
  TH1F* hisDeltaRorZfakePS_;
  TH1F* hisDeltaPhiRfake2S_;
  TH1F* hisDeltaRorZfake2S_;
  TProfile* profNsigmaPhiRvsInvPt_;
  TProfile* profNsigmaPhiRvsFracDist_;
  TH1F* hisDeltaBendTrue_;
  TH1F* hisDeltaBendFake_;
  TProfile* profFracTrueStubsVsLayer_;
  TProfile* profDupTracksVsTPeta_;

  // Histograms of track parameter resolution after HT transform.
  TH1F* hisQoverPtRes_;
  TH1F* hisPhi0Res_;
  TH1F* hisEtaRes_;
  TH1F* hisZ0Res_;

  // Diagnosis of failed tracking.
  TH1F* hisRecoFailureReason_;
  TH1F* hisRecoFailureLayer_;

  // Histograms used to make efficiency plots with track candidates prior to fit.
  TH1F* hisTPinvptForEff_;
  TH1F* hisRecoTPinvptForEff_;
  TH1F* hisTPetaForEff_;
  TH1F* hisRecoTPetaForEff_;
  TH1F* hisTPphiForEff_;
  TH1F* hisRecoTPphiForEff_;
  //
  TH1F* hisPerfRecoTPinvptForEff_;
  //
  TH1F* hisTPinvptForAlgEff_;
  TH1F* hisRecoTPinvptForAlgEff_;
  TH1F* hisTPetaForAlgEff_;
  TH1F* hisRecoTPetaForAlgEff_;
  TH1F* hisTPphiForAlgEff_;
  TH1F* hisRecoTPphiForAlgEff_;
  //
  TH1F* hisPerfRecoTPinvptForAlgEff_;
  //
  TH1F* hisTPd0ForAlgEff_;
  TH1F* hisRecoTPd0ForAlgEff_;
  TH1F* hisTPz0ForAlgEff_;
  TH1F* hisRecoTPz0ForAlgEff_;

  // Histograms for studying freak, large events with too many stubs.
  TH1F*     hisNumBusySecsInPerEvent_;
  TH1F*     hisNumBusySecsOutPerEvent_;
  TProfile* profFracBusyInVsEtaReg_;
  TProfile* profFracBusyOutVsEtaReg_;
  TProfile* profFracStubsKilledVsEtaReg_;
  TProfile* profFracTracksKilledVsEtaReg_;
  TProfile* profFracTracksKilledVsInvPt_;
  TProfile* profFracTPKilledVsEta_;
  TProfile* profFracTPKilledVsInvPt_;
  TH1F*     hisNumTPkilledBusySec_;
  map<string, TH1F*> hisNumInputStubs_;
  map<string, TH1F*> hisQoverPtInputStubs_;
  map<string, TH1F*> hisNumOutputStubs_;
  map<string, TH1F*> hisNumTracks_; 
  map<string, TH1F*> hisNumStubsPerTrack_; 
  map<string, TH1F*> hisTrackQoverPt_; 
  map<string, TH1F*> hisTrackPurity_; 
  map<string, TH1F*> hisNumTPphysics_; 
  map<string, TH1F*> hisNumTPpileup_; 
  map<string, TH1F*> hisSumPtTPphysics_; 
  map<string, TH1F*> hisSumPtTPpileup_; 

  // Histograms for track fitting evaluation, where map index specifies name of track fitting algorithm used.
  map<std::string, TH1F*> hisSeedQinvPt_;
  map<std::string, TH1F*> hisSeedPhi0_;
  map<std::string, TH1F*> hisSeedD0_;
  map<std::string, TH1F*> hisSeedZ0_;
  map<std::string, TH1F*> hisSeedEta_;

  map<std::string, TProfile*> profNumFittedCands_;

  map<std::string, TH1F*> hisFitQinvPtMatched_;
  map<std::string, TH1F*> hisFitPhi0Matched_;
  map<std::string, TH1F*> hisFitD0Matched_;
  map<std::string, TH1F*> hisFitZ0Matched_;
  map<std::string, TH1F*> hisFitEtaMatched_;

  map<std::string, TH1F*> hisFitChi2Matched_;
  map<std::string, TH1F*> hisFitChi2DofMatched_;

  map<std::string, TH1F*> hisFitQinvPtUnmatched_;
  map<std::string, TH1F*> hisFitPhi0Unmatched_;
  map<std::string, TH1F*> hisFitD0Unmatched_;
  map<std::string, TH1F*> hisFitZ0Unmatched_;
  map<std::string, TH1F*> hisFitEtaUnmatched_;

  map<std::string, TH1F*> hisFitChi2Unmatched_;
  map<std::string, TH1F*> hisFitChi2DofUnmatched_;

  map<std::string, TH2F*> hisFitVsTrueQinvPtGoodChi2_;
  map<std::string, TH2F*> hisFitVsTruePhi0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueD0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueZ0GoodChi2_;
  map<std::string, TH2F*> hisFitVsTrueEtaGoodChi2_;

  map<std::string, TH2F*> hisFitVsTrueQinvPtGenCand_;
  map<std::string, TH2F*> hisFitVsTruePhi0GenCand_;
  map<std::string, TH2F*> hisFitVsTrueD0GenCand_;
  map<std::string, TH2F*> hisFitVsTrueZ0GenCand_;
  map<std::string, TH2F*> hisFitVsTrueEtaGenCand_;

  map<std::string, TH1F*> hisFitQinvPtResGoodChi2_;
  map<std::string, TH1F*> hisFitPhi0ResGoodChi2_;
  map<std::string, TH1F*> hisFitD0ResGoodChi2_;
  map<std::string, TH1F*> hisFitZ0ResGoodChi2_;
  map<std::string, TH1F*> hisFitEtaResGoodChi2_;  

  map<std::string, TH1F*> hisTrueVsSeedQinvPtResGoodChi2_;
  map<std::string, TH1F*> hisTrueVsSeedPhi0ResGoodChi2_;
  map<std::string, TH1F*> hisTrueVsSeedD0ResGoodChi2_;
  map<std::string, TH1F*> hisTrueVsSeedZ0ResGoodChi2_;
  map<std::string, TH1F*> hisTrueVsSeedEtaResGoodChi2_;  

  map<std::string, TH2F*> hisFitVsTrueQinvPtFakeCand_;
  map<std::string, TH2F*> hisFitVsTruePhi0FakeCand_;
  map<std::string, TH2F*> hisFitVsTrueD0FakeCand_;
  map<std::string, TH2F*> hisFitVsTrueZ0FakeCand_;
  map<std::string, TH2F*> hisFitVsTrueEtaFakeCand_;

  map<std::string, TProfile*> hisPtResVsTrueEta_;
  map<std::string, TProfile*> hisPhi0ResVsTrueEta_;
  map<std::string, TProfile*> hisEtaResVsTrueEta_;
  map<std::string, TProfile*> hisZ0ResVsTrueEta_;
  map<std::string, TProfile*> hisD0ResVsTrueEta_;

  map<std::string, TProfile*> hisPtResVsTruePt_;
  map<std::string, TProfile*> hisPhi0ResVsTruePt_;
  map<std::string, TProfile*> hisEtaResVsTruePt_;
  map<std::string, TProfile*> hisZ0ResVsTruePt_;
  map<std::string, TProfile*> hisD0ResVsTruePt_;

  map<std::string, TH2F*> hisTrueFittedChiSquaredVsTrueEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredDofVsTrueEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredVsFittedEta_;
  map<std::string, TH2F*> hisTrueFittedChiSquaredDofVsFittedEta_;
  map<std::string, TH2F*> hisFittedChiSquaredFunctionOfStubs_;
  map<std::string, TH2F*> hisFittedChiSquaredDofFunctionOfStubs_;

  map<std::string, TH1F*> hisTrueEtaMatchedGoodChi2_;
  map<std::string, TH1F*> hisTrueEtaMatchedBadChi2_;
  map<std::string, TH1F*> hisStubPurityMatchedGoodChi2_;
  map<std::string, TH1F*> hisStubPurityMatchedBadChi2_;

  map<std::string, TProfile*> profChi2DofVsInvPtPERF_;
  map<std::string, TProfile*> profBigChi2DofVsInvPtPERF_;
  map<std::string, TH1F*>     hisD0TPBigChi2DofPERF_;
  map<std::string, TH1F*>     hisD0TPSmallChi2DofPERF_;

  map<std::string, TH2F*>     hisNumStubsKilledByFit_;
  map<std::string, TProfile*> profTrksKilledByFit_;
  map<std::string, TH2F*>     hisNumStubsVsPurity_;

  map<std::string, TH1F*> hisNumFittingIterations_;
  map<std::string, TH2F*> hisNumFittingIterationsVsPurity_;
  map<std::string, TH2F*> hisNumFittingIterationsVsPurityMatched_;
  map<std::string, TH2F*> hisNumFittingIterationsVsPurityUnmatched_;

  map<std::string, TH2F*> hisFitEfficiencyVsChi2Dof_;
  map<std::string, TH2F*> hisNumStubsVsChi2Dof_;
  map<std::string, TH2F*> hisNumLayersVsChi2Dof_;
  map<std::string, TH2F*> hisAvgNumStubsPerLayerVsChi2Dof_;

  // Histograms used for efficiency plots made with fitted tracks.
  map<std::string, TH1F*> hisFitTPinvptForEff_;
  map<std::string, TH1F*> hisFitTPetaForEff_;
  map<std::string, TH1F*> hisFitTPphiForEff_;
  map<std::string, TH1F*> hisPerfFitTPinvptForEff_;
  map<std::string, TH1F*> hisFitTPinvptForAlgEff_;
  map<std::string, TH1F*> hisFitTPetaForAlgEff_;
  map<std::string, TH1F*> hisFitTPphiForAlgEff_;
  map<std::string, TH1F*> hisPerfFitTPinvptForAlgEff_;
  map<std::string, TH1F*> hisFitTPd0ForAlgEff_;
  map<std::string, TH1F*> hisFitTPz0ForAlgEff_;

  // Histograms of tracking efficiency & fake rate after Hough transform based on tracks prior to track fit.
  TGraphAsymmErrors* graphEffVsInvPt_;
  TGraphAsymmErrors* graphEffVsEta_;
  TGraphAsymmErrors* graphEffVsPhi_;
  //
  TGraphAsymmErrors* graphPerfEffVsInvPt_;
  //
  TGraphAsymmErrors* graphAlgEffVsInvPt_;
  TGraphAsymmErrors* graphAlgEffVsEta_;
  TGraphAsymmErrors* graphAlgEffVsPhi_;
  //
  TGraphAsymmErrors* graphPerfAlgEffVsInvPt_;
  //
  TGraphAsymmErrors* graphAlgEffVsD0_;
  TGraphAsymmErrors* graphAlgEffVsZ0_;

  // Histograms of tracking efficiency & fake rate after Hough transform based on tracks after the track fit.
  map<std::string, TGraphAsymmErrors*> graphEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphEffFitVsEta_;
  map<std::string, TGraphAsymmErrors*> graphEffFitVsPhi_;
  //
  map<std::string, TGraphAsymmErrors*> graphPerfEffFitVsInvPt_;
  //
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsInvPt_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsEta_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsPhi_;
  //
  map<std::string, TGraphAsymmErrors*> graphPerfAlgEffFitVsInvPt_;
  //
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsD0_;
  map<std::string, TGraphAsymmErrors*> graphAlgEffFitVsZ0_;

  // Number of perfectly reconstructed tracks amongst TP used for algorithmic efficiency measurement.
  // Perfectly means that all stubs on track were produced by same TP.
  unsigned int numPerfRecoTPforAlg_;

  // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted.
  map<std::string, unsigned int> numFitAlgEff_;
  map<std::string, unsigned int> numFitPerfAlgEff_;

  // Number of genuine reconstructed and perfectly reconstructed tracks which were fitted post-cut.
  map<std::string, unsigned int> numFitAlgEffPass_;
  map<std::string, unsigned int> numFitPerfAlgEffPass_;

  // For filling histograms related to the stub-pairs algorithm
  void analyse_PairFinding   (const vector<const Stub*>& vStubs);
  void analyse_cuts          (const vector<const Stub*>& vStubs);
  void analyse_Formulae(const vector<const Stub*>& vStubs);
  vector<const Stub*> depair(vector< pair<const Stub*, const Stub*> > vPairs) const;
};

#endif

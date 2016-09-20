#include "TMTrackTrigger/TMTrackFinder/interface/Histos.h"
#include "TMTrackTrigger/TMTrackFinder/interface/KillOverlapStubs.h"
#include "TMTrackTrigger/TMTrackFinder/interface/TP.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TGraphAsymmErrors.h>

#include <algorithm>
#include <array>
#include <unordered_set>
#include <set>


void Histos::bookStubPairs() {
  // for testing the formula for z0 and q/pt from the two stubs
  TFileDirectory formulaeDir = fs_->mkdir("TwoFormulae");
  hisPtPairMinusTruth_        = formulaeDir.make<TH1F>("PtPairMinusTruth","; PT.q/pt - pair.q/pt;"                  , 150, -1.0 , 1.0                  );
  hisZ0PairMinusTruth_        = formulaeDir.make<TH1F>("Z0PairMinusTruth","; PT.z0 - pair.z0;"                      , 150, -20.0, 20.0                 );
  hisZ0PairVsTruth_           = formulaeDir.make<TH2F>("Z0PairVsTruth","; pair z0; truth z0"                        , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisPtPairVsTruth_           = formulaeDir.make<TH2F>("PtPairVsTruth","; pair q/pt; truth q/pt"                    , 150, -1.0 , 1.0 , 50, -1.0 , 1.0 );
  hisZ0PairVsTruth_barrel_    = formulaeDir.make<TH2F>("Z0PairVsTruth_barrel","barrel; pair z0; truth z0"           , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_endcap_    = formulaeDir.make<TH2F>("Z0PairVsTruth_endcap","endcap; pair z0; truth z0"           , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_PS_        = formulaeDir.make<TH2F>("Z0PairVsTruth_PS","PS modules; pair z0; truth z0"           , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_2S_        = formulaeDir.make<TH2F>("Z0PairVsTruth_2S","2S modules; pair z0; truth z0"           , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_barrel_PS_ = formulaeDir.make<TH2F>("Z0PairVsTruth_barrel_PS","barrel, PS; pair z0; truth z0"    , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_barrel_2S_ = formulaeDir.make<TH2F>("Z0PairVsTruth_barrel_2S","barrel, 2S; pair z0; truth z0"    , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_endcap_PS_ = formulaeDir.make<TH2F>("Z0PairVsTruth_endcap_PS","endcap, PS; pair z0; truth z0"    , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisZ0PairVsTruth_endcap_2S_ = formulaeDir.make<TH2F>("Z0PairVsTruth_endcap_2S","endcap, 2S; pair z0; truth z0"    , 150, -20.0, 20.0, 50, -20.0, 20.0);
  hisPtPairVsTruth_barrel_    = formulaeDir.make<TH2F>("PtPairVsTruth_barrel","barrel; pair q/pt; truth q/pt"       , 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_endcap_    = formulaeDir.make<TH2F>("PtPairVsTruth_endcap","endcap; pair q/pt; truth q/pt"       , 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_PS_        = formulaeDir.make<TH2F>("PtPairVsTruth_PS","PS modules; pair q/pt; truth q/pt"       , 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_2S_        = formulaeDir.make<TH2F>("PtPairVsTruth_2S","2S modules; pair q/pt; truth q/pt"       , 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_barrel_PS_ = formulaeDir.make<TH2F>("PtPairVsTruth_barrel_PS","barrel, PS; pair q/pt; truth q/pt", 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_barrel_2S_ = formulaeDir.make<TH2F>("PtPairVsTruth_barrel_2S","barrel, 2S; pair q/pt; truth q/pt", 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_endcap_PS_ = formulaeDir.make<TH2F>("PtPairVsTruth_endcap_PS","endcap, PS; pair q/pt; truth q/pt", 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisPtPairVsTruth_endcap_2S_ = formulaeDir.make<TH2F>("PtPairVsTruth_endcap_2S","endcap, 2S; pair q/pt; truth q/pt", 150, -1.0,  1.0,  50, -1.0,  1.0 );
  hisZ0PairMinusTruth_barrel_PS_ = formulaeDir.make<TH1F>("Z0PairMinusTruth_barrel_PS","barrel, PS; (pair z0) - (truth z0)", 150, -20.0, 20.0);
  hisZ0PairMinusTruth_barrel_PS_abs_ = formulaeDir.make<TH1F>("Z0PairMinusTruth_barrel_PS_abs","barrel, PS; abs((pair z0) - (truth z0))", 100, 0.0, 20.0);

  // stats on my pair-finding algorithm;
  TFileDirectory pairsDir = fs_->mkdir("StubPairs");
  //Pt
  his_StubsInFound_Pt      = pairsDir.make<TH1F>("his_StubsInFound_Pt","number of stubs in found pairs;q/Pt of TP;"                                 ,50,-1.0,1.0);
  his_AllStubs_Pt          = pairsDir.make<TH1F>("his_AllStubs_Pt","number of all stubs;q/Pt of TP;"                                                ,50,-1.0,1.0);
  his_AllStubs_AbsPt       = pairsDir.make<TH1F>("his_AllStubs_AbsPt","number of all stubs;abs(q/Pt) of TP;"                                        ,50, 0.0,1.0);
  his_TrueFoundPairs_Pt    = pairsDir.make<TH1F>("his_TrueFoundPairs_Pt","number of true found pairs;q/Pt of first common TP;"                      ,50,-1.0,1.0);
  his_TrueFoundStubs_Pt    = pairsDir.make<TH1F>("his_TrueFoundStubs_Pt","stubs in true found pairs;q/Pt of first TP;"                              ,50,-1.0,1.0);
  his_TrueFoundStubs_AbsPt = pairsDir.make<TH1F>("his_TrueFoundStubs_AbsPt","stubs in true found pairs;abs(q/Pt) of first TP;"                      ,50,0.0,1.0);
  his_WrongFoundPairs_Pt   = pairsDir.make<TH1F>("his_WrongFoundPairs_Pt","number of wrong found pairs;q/Pt of first TP of first stub;"             ,50,-1.0,1.0);
  his_AllTruePairs_Pt      = pairsDir.make<TH1F>("his_AllTruePairs_Pt","number of all true pairs;q/Pt of first common TP;"                          ,50,-1.0,1.0);
  his_AllTrueStubs_Pt      = pairsDir.make<TH1F>("his_AllTrueStubs_Pt","stubs in all true pairs;q/Pt of first TP;"                                  ,50,-1.0,1.0);
  his_AllTrueStubs_AbsPt   = pairsDir.make<TH1F>("his_AllTrueStubs_AbsPt","stubs in all true pairs;abs(q/Pt) of first TP;"                          ,50,0,1.0);
  his_AllFoundPairs_Pt     = pairsDir.make<TH1F>("his_AllFoundPairs_Pt","number of all found pairs;q/Pt of TP of first stub;"                       ,50,-1.0,1.0);
  his_AllFoundStubs_Pt     = pairsDir.make<TH1F>("his_AllFoundStubs_Pt","stubs in all found pairs;q/Pt of first TP;"                                ,50,-1.0,1.0);
  his_AllFoundStubs_AbsPt  = pairsDir.make<TH1F>("his_AllFoundStubs_AbsPt","stubs in all found pairs;abs(q/Pt) of first TP;"                        ,50,0,1.0);
  his_WrongStubs_Pt        = pairsDir.make<TH1F>("his_WrongStubs_Pt","number of stubs in different-TP found pairs;q/Pt of TP of first stub;"        ,50,-1.0,1.0);
  his_WrongStubs_AbsPt     = pairsDir.make<TH1F>("his_WrongStubs_AbsPt","number of stubs in different-TP found pairs;abs(q/Pt) of TP of first stub;",50,0.0,1.0);
  //Eta
  his_StubsInFound_Eta   = pairsDir.make<TH1F>("his_StubsInFound_Eta","number of stubs in found pairs;eta of TP;"           ,50,0.0,2.5);
  his_AllStubs_Eta       = pairsDir.make<TH1F>("his_AllStubs_Eta","number of all stubs;eta of TP;"                          ,50,0.0,2.5);
  his_TrueFoundPairs_Eta = pairsDir.make<TH1F>("his_TrueFoundPairs_Eta","number of true found pairs;eta of first common TP;",50,0.0,2.5);
  his_TrueFoundStubs_Eta = pairsDir.make<TH1F>("his_TrueFoundStubs_Eta","stubs in true found pairs;eta of first TP;",50,0.0,2.5);
  his_WrongFoundPairs_Eta = pairsDir.make<TH1F>("his_WrongFoundPairs_Eta","number of wrong found pairs;eta of first TP of first stub;",50,0.0,2.5);
  his_AllTruePairs_Eta   = pairsDir.make<TH1F>("his_AllTruePairs_Eta","number of all true pairs;eta of first common TP;"    ,50,0.0,2.5);
  his_AllTrueStubs_Eta   = pairsDir.make<TH1F>("his_AllTrueStubs_Eta","stubs in all true pairs;eta of first TP;"    ,50,0.0,2.5);
  his_AllFoundPairs_Eta  = pairsDir.make<TH1F>("his_AllFoundPairs_Eta","number of all found pairs;eta of TP of first stub;" ,50,0.0,2.5);
  his_AllFoundStubs_Eta  = pairsDir.make<TH1F>("his_AllFoundStubs_Eta","stubs in all found pairs;eta of first TP;" ,50,0.0,2.5);
  his_WrongStubs_Eta     = pairsDir.make<TH1F>("his_WrongStubs_Eta","number of stubs in different-TP found pairs;eta of TP of first stub;",50,0.0,2.5);
  //r-z
  his_StubsInFound_Loc   = pairsDir.make<TH2F>("his_StubsInFound_Loc","number of stubs in found pairs;abs(z);r"                                          ,50,0.0,200.0,50,0.0,120.0);
  his_AllStubs_Loc       = pairsDir.make<TH2F>("his_AllStubs_Loc","number of all stubs;abs(z);r"                                                         ,50,0.0,200.0,50,0.0,120.0);
  his_TrueFoundPairs_Loc = pairsDir.make<TH2F>("his_TrueFoundPairs_Loc","number of true found pairs;abs(z) of first stub in pair;r of first stub in pair",50,0.0,200.0,50,0.0,120.0);
  his_TrueFoundStubs_Loc = pairsDir.make<TH2F>("his_TrueFoundStubs_Loc","stubs in true found pairs;abs(z);r",50,0.0,200.0,50,0.0,120.0);
  his_WrongFoundPairs_Loc = pairsDir.make<TH2F>("his_WrongFoundPairs_Loc","number of wrong found pairs;abs(z) of first stub in pair;r of first stub in pair",50,0.0,200.0,50,0.0,120.0);
  his_AllTruePairs_Loc   = pairsDir.make<TH2F>("his_AllTruePairs_Loc","number of all true pairs;abs(z) of first stub in pair;r of first stub in pair"    ,50,0.0,200.0,50,0.0,120.0);
  his_AllTrueStubs_Loc   = pairsDir.make<TH2F>("his_AllTrueStubs_Loc","stubs in all true pairs;abs(z) of first stub in pair;r of first stub in pair"    ,50,0.0,200.0,50,0.0,120.0);
  his_AllFoundPairs_Loc  = pairsDir.make<TH2F>("his_AllFoundPairs_Loc","number of all found pairs;abs(z) of first stub in pair;r of first stub in pair"  ,50,0.0,200.0,50,0.0,120.0);
  his_AllFoundStubs_Loc  = pairsDir.make<TH2F>("his_AllFoundStubs_Loc","stubs in all found pairs;abs(z);r"  ,50,0.0,200.0,50,0.0,120.0);
  his_WrongStubs_Loc     = pairsDir.make<TH2F>("his_WrongStubs_Loc","number of stubs in different-TP found pairs;abs(z) of first stub in pair;r of first stub in pair",50,0.0,200.0,50,0.0,120.0);

  // determining the cuts
  his_pt_cut_AllStubs               = pairsDir.make<TH1F>("his_pt_cut_AllStubs","number of all stubs;pt_cut;",50,0,4);
  his_pt_cut_StubsInFoundPairs      = pairsDir.make<TH1F>("his_pt_cut_StubsInFoundPairs","number of stubs in found pairs;pt_cut;",50,0,4);
  his_pt_cut_StubsInTrueFoundPairs  = pairsDir.make<TH1F>("his_pt_cut_StubsInTrueFoundPairs","number of stubs in true found pairs;pt_cut;",50,0,4);
  his_pt_cut_StubsInWrongFoundPairs = pairsDir.make<TH1F>("his_pt_cut_StubsInWrongFoundPairs","number of pt>3.00 stubs in wrong found pairs;pt_cut;",50,0,4);
  his_pt_cut_StubsInAllTruePairs    = pairsDir.make<TH1F>("his_pt_cut_StubsInAllTruePairs","number of stubs in all true pairs;pt_cut;",50,0,4);
  his_z0_cut_AllStubs               = pairsDir.make<TH1F>("his_z0_cut_AllStubs","number of all stubs;z0_cut;",50,0,60);
  his_z0_cut_StubsInFoundPairs      = pairsDir.make<TH1F>("his_z0_cut_StubsInFoundPairs","number of stubs in found pairs;z0_cut;",50,0,60);
  his_z0_cut_StubsInTrueFoundPairs  = pairsDir.make<TH1F>("his_z0_cut_StubsInTrueFoundPairs","number of stubs in true found pairs;z0_cut;",50,0,60);
  his_z0_cut_StubsInWrongFoundPairs = pairsDir.make<TH1F>("his_z0_cut_StubsInWrongFoundPairs","number of pt>3.00 stubs in wrong found pairs;z0_cut;",50,0,60);
  his_z0_cut_StubsInAllTruePairs    = pairsDir.make<TH1F>("his_z0_cut_StubsInAllTruePairs","number of stubs in all true pairs;z0_cut;",50,0,60);
}

// Fill histograms relating the stub-pairs algorithm
void Histos::fillStubPairs(const vector<const Stub*>& vStubs)
{
  analyse_Formulae(vStubs);
  analyse_PairFinding(vStubs);
  analyse_cuts(vStubs);
}

// Stats on the input data
void Histos::analyse_Formulae(const vector<const Stub*>& vStubs)
{
  // The formulae for z0 and ptOverQ of a pair of stubs having any q/pt
  KillOverlapStubs killOverlapStubs_ = KillOverlapStubs(vStubs, settings_, -1, 15.0);
  vector< pair<const Stub*, const Stub*> > true_pairs  = killOverlapStubs_.getPairs("truePairFinder");
  for ( auto p : true_pairs ) {
    const Stub* s1   = p.first; const Stub* s2 = p.second;
    double z0_f      = killOverlapStubs_.getTrackParams(s1,s2).first;
    double ptOverQ_f = killOverlapStubs_.getTrackParams(s1,s2).second;
    double z0_t      = killOverlapStubs_.getCommonTP(s1,s2)->z0();
    double qOverPt_t = killOverlapStubs_.getCommonTP(s1,s2)->qOverPt();
    hisZ0PairVsTruth_    -> Fill( z0_f, z0_t              );
    hisPtPairVsTruth_    -> Fill( 1/ptOverQ_f, qOverPt_t  );
    hisPtPairMinusTruth_ -> Fill( qOverPt_t - 1/ptOverQ_f );
    hisZ0PairMinusTruth_ -> Fill( z0_t - z0_f );
    if ( s1->barrel() ) {
      hisZ0PairVsTruth_barrel_ -> Fill( z0_f, z0_t );
      hisPtPairVsTruth_barrel_ -> Fill( 1/ptOverQ_f, qOverPt_t );
      if ( s1->psModule() ) {
        hisZ0PairVsTruth_barrel_PS_ -> Fill( z0_f, z0_t );
        hisPtPairVsTruth_barrel_PS_ -> Fill( 1/ptOverQ_f, qOverPt_t );
        hisZ0PairMinusTruth_barrel_PS_ -> Fill( z0_f - z0_t );
        hisZ0PairMinusTruth_barrel_PS_abs_ -> Fill( fabs(z0_f - z0_t) );
      } else { // 2S
        hisZ0PairVsTruth_barrel_2S_ -> Fill( z0_f, z0_t );
        hisPtPairVsTruth_barrel_2S_ -> Fill( 1/ptOverQ_f, qOverPt_t );
      }
    } else { // endcap
      hisZ0PairVsTruth_endcap_ -> Fill( z0_f, z0_t );
      hisPtPairVsTruth_endcap_ -> Fill( 1/ptOverQ_f, qOverPt_t );
      if ( s1->psModule() ) {
        hisZ0PairVsTruth_endcap_PS_ -> Fill( z0_f, z0_t );
        hisPtPairVsTruth_endcap_PS_ -> Fill( 1/ptOverQ_f, qOverPt_t );
      } else { // 2S
        hisZ0PairVsTruth_endcap_2S_ -> Fill( z0_f, z0_t );
        hisPtPairVsTruth_endcap_2S_ -> Fill( 1/ptOverQ_f, qOverPt_t );
      }
    }
    if ( s1->psModule() ) {
      hisZ0PairVsTruth_PS_ -> Fill( z0_f, z0_t );
      hisPtPairVsTruth_PS_ -> Fill( 1/ptOverQ_f, qOverPt_t );
    } else { // 2S
      hisZ0PairVsTruth_2S_ -> Fill( z0_f, z0_t );
      hisPtPairVsTruth_2S_ -> Fill( 1/ptOverQ_f, qOverPt_t );
    }
  }
}

// Stats on my pair-finding algorithm;
void Histos::analyse_PairFinding(const vector<const Stub*>& vStubs)
{
  // filter out non-genuine stubs
  vector<const Stub*> vStubs_filt;
  for (const Stub* s : vStubs)
    if ( s->genuine() )
      vStubs_filt.push_back(s);

  KillOverlapStubs killOverlapStubs_ = KillOverlapStubs(vStubs_filt, settings_);
  vector< pair<const Stub*, const Stub*> > true_pairs  = killOverlapStubs_.getPairs("truePairFinder");
  vector< pair<const Stub*, const Stub*> > found_pairs = killOverlapStubs_.getPairs("pairFinder");

  for (const Stub* s : vStubs_filt) {
    his_AllStubs_Pt    -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
    his_AllStubs_AbsPt -> Fill( fabs(killOverlapStubs_.getFirstTP(s)->qOverPt()) );
    his_AllStubs_Eta   -> Fill( killOverlapStubs_.getFirstTP(s)->eta()     );
    his_AllStubs_Loc   -> Fill( fabs(s->z()), s->r()  );
  }

  for ( const Stub* s : depair(found_pairs) ) {
    his_StubsInFound_Pt  -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
    his_StubsInFound_Eta -> Fill( killOverlapStubs_.getFirstTP(s)->eta()     );
    his_StubsInFound_Loc -> Fill( fabs(s->z()), s->r()  );
  }

  set<const Stub*> true_stubs;

  for (auto p : true_pairs) {
    his_AllTruePairs_Pt  -> Fill( killOverlapStubs_.getCommonTP(p.first,p.second)->qOverPt() );
    his_AllTruePairs_Eta -> Fill( killOverlapStubs_.getCommonTP(p.first,p.second)->eta()     );
    his_AllTruePairs_Loc -> Fill( fabs(p.first->z()), p.first->r()      );
    true_stubs.insert(p.first);
    true_stubs.insert(p.second);
  }

  for (const Stub* s : true_stubs) {
    his_AllTrueStubs_Pt    -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
    his_AllTrueStubs_AbsPt -> Fill( fabs(killOverlapStubs_.getFirstTP(s)->qOverPt()) );
    his_AllTrueStubs_Eta   -> Fill( killOverlapStubs_.getFirstTP(s)->eta() );
    his_AllTrueStubs_Loc   -> Fill( fabs(s->z()), s->r()  );

    set<const Stub*> wrong_stubs;
    set<const Stub*> true_found_stubs;
    set<const Stub*> found_stubs;

    for (auto p : found_pairs) {
      if ( killOverlapStubs_.getCommonTP(p.first, p.second) ) {
        his_TrueFoundPairs_Pt  -> Fill( killOverlapStubs_.getCommonTP(p.first,p.second)->qOverPt() );
        his_TrueFoundPairs_Eta -> Fill( killOverlapStubs_.getCommonTP(p.first,p.second)->eta()     );
        his_TrueFoundPairs_Loc -> Fill( fabs(p.first->z()), p.first->r()      );
        true_found_stubs.insert(p.first);
        true_found_stubs.insert(p.second);
      } else {
        his_WrongFoundPairs_Pt  -> Fill( killOverlapStubs_.getFirstTP(p.first)->qOverPt() );
        his_WrongFoundPairs_Eta -> Fill( killOverlapStubs_.getFirstTP(p.first)->eta()     );
        his_WrongFoundPairs_Loc -> Fill( fabs(p.first->z()), p.first->r()      );
        wrong_stubs.insert(p.first);
        wrong_stubs.insert(p.second);
      }
      his_AllFoundPairs_Pt  -> Fill( killOverlapStubs_.getFirstTP(p.first) -> qOverPt()    );
      his_AllFoundPairs_Eta -> Fill( killOverlapStubs_.getFirstTP(p.first) -> eta()        );
      his_AllFoundPairs_Loc -> Fill( fabs(p.first->z()), p.first->r() );
      found_stubs.insert(p.first);
      found_stubs.insert(p.second);
    }

    for (const Stub* s : wrong_stubs) {
      his_WrongStubs_Pt    -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
      his_WrongStubs_AbsPt -> Fill( fabs(killOverlapStubs_.getFirstTP(s)->qOverPt()) );
      his_WrongStubs_Eta   -> Fill( killOverlapStubs_.getFirstTP(s)->eta()     );
      his_WrongStubs_Loc   -> Fill( fabs(s->z()), s->r()  );
    }

    for (const Stub* s : true_found_stubs) {
      his_TrueFoundStubs_Pt    -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
      his_TrueFoundStubs_AbsPt -> Fill( fabs(killOverlapStubs_.getFirstTP(s)->qOverPt()) );
      his_TrueFoundStubs_Eta   -> Fill( killOverlapStubs_.getFirstTP(s)->eta() );
      his_TrueFoundStubs_Loc   -> Fill( fabs(s->z()), s->r()  );
    }

    for (const Stub* s : found_stubs) {
      his_AllFoundStubs_Pt    -> Fill( killOverlapStubs_.getFirstTP(s)->qOverPt() );
      his_AllFoundStubs_AbsPt -> Fill( fabs(killOverlapStubs_.getFirstTP(s)->qOverPt()) );
      his_AllFoundStubs_Eta   -> Fill( killOverlapStubs_.getFirstTP(s)->eta() );
      his_AllFoundStubs_Loc   -> Fill( fabs(s->z()), s->r()  );
    }
  }
}


// Determining ideal cuts
void Histos::analyse_cuts(const vector<const Stub*>& vStubs)
{

  // filter out non-genuine stubs
  vector<const Stub*> vStubs_filt;
  for (const Stub* s : vStubs)
    if ( s->genuine() )
      vStubs_filt.push_back(s);

  ////////////////////////////////////////////////////////////////////////////////////////////
  //////   Analyse the pt_cut                             ////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  for (double pt_cut=4.0/100; pt_cut<4.0; pt_cut+=4.0/50) {
    KillOverlapStubs killOverlapStubs_ = KillOverlapStubs(vStubs_filt, settings_, pt_cut, 15.0);
    vector< pair<const Stub*, const Stub*> > true_pairs  = killOverlapStubs_.getPairs("truePairFinder");
    vector< pair<const Stub*, const Stub*> > found_pairs = killOverlapStubs_.getPairs("pairFinder");

    // find the number of all stubs
    size_t numOfStubs = vStubs_filt.size();
    his_pt_cut_AllStubs -> Fill( pt_cut, numOfStubs );

    // find the number of stubs in found pairs
    size_t numOfFoundStubs = depair(found_pairs).size();
    his_pt_cut_StubsInFoundPairs -> Fill( pt_cut, numOfFoundStubs);

    // find the number of stubs in true and wrong found pairs
    size_t numOfTrueFoundStubs = 0;
    size_t numOfWrongHighPtStubs = 0;
    for (auto p : found_pairs)
      if ( killOverlapStubs_.getCommonTP(p.first, p.second) ) {
        ++numOfTrueFoundStubs;
        ++numOfTrueFoundStubs;
      } else {
        if (killOverlapStubs_.getFirstTP(p.first)->pt() >= 3.00 )
          ++numOfWrongHighPtStubs;
        if (killOverlapStubs_.getFirstTP(p.second)->pt() >= 3.00 )
          ++numOfWrongHighPtStubs;
      }
    his_pt_cut_StubsInTrueFoundPairs  -> Fill( pt_cut, numOfTrueFoundStubs );
    his_pt_cut_StubsInWrongFoundPairs -> Fill( pt_cut, numOfWrongHighPtStubs );

    // find the number of stubs in all true pairs
    size_t numOfTrueStubs = depair(true_pairs).size();
    his_pt_cut_StubsInAllTruePairs -> Fill( pt_cut, numOfTrueStubs );
  }


  ////////////////////////////////////////////////////////////////////////////////////////////
  //////   Analyse the z0_cut                             ////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////

  for (double z0_cut=60.0/100; z0_cut<60.0; z0_cut+=60.0/50) {
    KillOverlapStubs killOverlapStubs_ = KillOverlapStubs(vStubs_filt, settings_, 3.0, z0_cut);
    vector< pair<const Stub*, const Stub*> > true_pairs  = killOverlapStubs_.getPairs("truePairFinder");
    vector< pair<const Stub*, const Stub*> > found_pairs = killOverlapStubs_.getPairs("pairFinder");

    // find the number of all stubs
    size_t numOfStubs = vStubs_filt.size();
    his_z0_cut_AllStubs -> Fill( z0_cut, numOfStubs );

    // find the number of stubs in found pairs
    size_t numOfFoundStubs = depair(found_pairs).size();
    his_z0_cut_StubsInFoundPairs -> Fill( z0_cut, numOfFoundStubs);

    // find the number of stubs in true and wrong found pairs
    size_t numOfTrueFoundStubs = 0;
    size_t numOfWrongHighPtStubs = 0;
    for (auto p : found_pairs)
      if ( killOverlapStubs_.getCommonTP(p.first, p.second) ) {
        ++numOfTrueFoundStubs;
        ++numOfTrueFoundStubs;
      } else {
        if (killOverlapStubs_.getFirstTP(p.first)->pt() >= 3.00 )
          ++numOfWrongHighPtStubs;
        if (killOverlapStubs_.getFirstTP(p.second)->pt() >= 3.00 )
          ++numOfWrongHighPtStubs;
      }
    his_z0_cut_StubsInTrueFoundPairs  -> Fill( z0_cut, numOfTrueFoundStubs );
    his_z0_cut_StubsInWrongFoundPairs -> Fill( z0_cut, numOfWrongHighPtStubs );

    // find the number of stubs in all true pairs
    size_t numOfTrueStubs = depair(true_pairs).size();
    his_z0_cut_StubsInAllTruePairs -> Fill( z0_cut, numOfTrueStubs );
  }
}


// convert a vector of pairs into a set, then vector, of elements
vector<const Stub*> Histos::depair(vector< pair<const Stub*, const Stub*> > vPairs) const
{
  set<const Stub*> sElems;
  for (auto p : vPairs) {
    sElems.insert(p.first);
    sElems.insert(p.second);
  }
  vector<const Stub*> vElems;
  for (auto s : sElems)
    vElems.push_back(s);
  return vElems;
}

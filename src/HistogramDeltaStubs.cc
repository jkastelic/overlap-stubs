#include "TMTrackTrigger/TMTrackFinder/interface/Histos.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <TH1F.h>
#include <TH2F.h>

#include <unordered_map>


// Book histograms to do with nearby stubs in a module (which may be produced by delta rays)
void Histos::bookDeltaStubs() {
  TFileDirectory deltaDir = fs_->mkdir("DeltaStubs");

  // on all modules in the detector
  his_delta_U    = deltaDir.make<TH1F>("his_delta_U","same-module pairs of stubs;abs(u1-u2);",50,0,1000);
  his_delta_V    = deltaDir.make<TH1F>("his_delta_V","same-module pairs of stubs;abs(v1-v2);",50,0,1000);
  his_delta_U_hr = deltaDir.make<TH1F>("his_delta_U_hr","same-module pairs of stubs;abs(u1-u2);",100,0,20);
  his_delta_V_hr = deltaDir.make<TH1F>("his_delta_V_hr","same-module pairs of stubs;abs(v1-v2);",100,0,20);

  // barrel PS module stubs only
  his_delta_U_bPS        = deltaDir.make<TH1F>("his_delta_U_bPS","same-module pairs of stubs;abs(u1-u2);",50,0,1000);
  his_delta_V_bPS        = deltaDir.make<TH1F>("his_delta_V_bPS","same-module pairs of stubs;abs(v1-v2);",50,0,1000);
  his_delta_U_hr_bPS     = deltaDir.make<TH1F>("his_delta_U_hr_bPS","same-module pairs of stubs;abs(u1-u2);",100,0,3);
  his_delta_V_hr_bPS     = deltaDir.make<TH1F>("his_delta_V_hr_bPS","same-module pairs of stubs;abs(v1-v2);",100,0,20);
  his_delta_hr_bPS       = deltaDir.make<TH2F>("his_delta_hr_bPS","same-module pairs of stubs;abs(u1-u2);abs(v1-v2)",100,0,3,100,0,16);
  his_delta_bPS          = deltaDir.make<TH2F>("his_delta_bPS","same-module pairs of stubs;abs(u1-u2);abs(v1-v2)",50,0,100,50,0,20);
  his_delta_r_bPS        = deltaDir.make<TH1F>("his_delta_r_bPS","same-module pairs of stubs;sqrt( (u1-u2)^2+(v1-v2)^2 );",50,0,1000);
  his_delta_r_hr_bPS     = deltaDir.make<TH1F>("his_delta_r_hr_bPS","same-module pairs of stubs;sqrt( (u1-u2)^2+(v1-v2)^2 );",100,0,20);
  his_delta_r_hhr_bPS    = deltaDir.make<TH1F>("his_delta_r_hhr_bPS","same-module pairs of stubs;sqrt( (u1-u2)^2+(v1-v2)^2 );",100,0,2);
  his_delta_global_r_bPS = deltaDir.make<TH1F>("his_delta_global_r_bPS","same-module pairs of stubs;global_r;",50,0,20);
  his_delta_global_r_hr_bPS = deltaDir.make<TH1F>("his_delta_global_r_hr_bPS","same-module pairs of stubs;global_r;",50,0,0.2);
}

// Fill histograms to do with nearby stubs in a module (which may be produced by delta rays)
void Histos::fillDeltaStubs(const vector<const Stub*>& vStubs)
{
  // sort stubs into modules
  std::unordered_map<unsigned, vector<const Stub*>> modules;
  for (const Stub* s : vStubs)
    modules[s->idDet()].push_back(s);

  size_t deltaCount_U=0, deltaCount_V=0;

  for (auto m : modules)
    for (const Stub* s1 : m.second) {

      // calculate stub 1 position by averaging cluster positions

      float u1 = ( s1->localU_cluster()[0] + s1->localU_cluster()[1] ) / 2;
      float v1 = ( s1->localV_cluster()[0] + s1->localV_cluster()[1] ) / 2;
      for (const Stub* s2 : m.second) {

        // calculate stub 2 position by averaging cluster positions

        float u2 = ( s2->localU_cluster()[0] + s2->localU_cluster()[1] ) / 2;
        float v2 = ( s2->localV_cluster()[0] + s2->localV_cluster()[1] ) / 2;

        // histogram how close together they are

        his_delta_bPS -> Fill( abs(u1-u2), abs(v1-v2) );
        his_delta_hr_bPS -> Fill( abs(u1-u2), abs(v1-v2) );

        float r  = sqrt( (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2) );
        if ( r != 0 )
          if ( s1->psModule() && s1->barrel()) {
            his_delta_r_bPS -> Fill(r);
            his_delta_r_hr_bPS -> Fill(r);
            his_delta_r_hhr_bPS -> Fill(r);
          }

        // convert to cartesian coordinates; units are presumably centimeters
        double x1 = s1->r() * cos(s1->phi()); double x2 = s2->r() * cos(s2->phi());
        double y1 = s1->r() * sin(s1->phi()); double y2 = s2->r() * sin(s2->phi());
        double z1 = s1->z(); double z2 = s2->z();
        // find distance between two stubs, in centimeters
        float global_r = sqrt( pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2) );
        if ( global_r != 0 )
          if ( s1->psModule() && s1->barrel()) {
            his_delta_global_r_bPS    -> Fill(global_r);
            his_delta_global_r_hr_bPS -> Fill(global_r);
          }

        if ( abs(u1-u2) != 0.0 ) {
          his_delta_U    -> Fill( abs(u1-u2) );
          his_delta_U_hr -> Fill( abs(u1-u2) );
          if ( s1->psModule() && s1->barrel()) {
            his_delta_U_bPS    -> Fill( abs(u1-u2) );
            his_delta_U_hr_bPS -> Fill( abs(u1-u2) );
          }
        }

        if ( abs(v1-v2) != 0.0 ) {
          his_delta_V    -> Fill( abs(v1-v2) );
          his_delta_V_hr -> Fill( abs(v1-v2) );
          if ( s1->psModule() && s1->barrel()) {
            his_delta_V_bPS    -> Fill( abs(v1-v2) );
            his_delta_V_hr_bPS -> Fill( abs(v1-v2) );
          }
        }
      }
    }

  // analyse the delta-stub removal algorithm, on barrel PS modules
  std::vector<const Stub*> vStubs_bPS;
  for (const Stub* s : vStubs)
    if ( s->barrel() )
      vStubs_bPS.push_back(s);
  std::vector<const Stub*> vStubs_filtered = deltaKiller(vStubs_bPS, 2.0);
  std::cout << "deltas!!! " << vStubs_bPS.size() << " " << vStubs_filtered.size() << "\n";
}

std::vector<const Stub*> Histos::deltaKiller(const vector<const Stub*>& vStubs, const double epsilon) const
{
  // sort stubs into modules
  std::unordered_map<unsigned, vector<const Stub*>> modules;
  for (const Stub* s : vStubs)
    modules[s->idDet()].push_back(s);

  // build map
  std::unordered_map<const Stub*, bool> stub_killed;
  for (auto m : modules)
    for (const Stub* s : m.second)
      stub_killed[s] = false;

  // find stubs due to delta rays
  for (auto m : modules)
    for (const Stub* s1 : m.second)
      for (const Stub* s2 : m.second) {
        float u1 = ( s1->localU_cluster()[0] + s1->localU_cluster()[1] ) / 2;
        float v1 = ( s1->localV_cluster()[0] + s1->localV_cluster()[1] ) / 2;
        float u2 = ( s2->localU_cluster()[0] + s2->localU_cluster()[1] ) / 2;
        float v2 = ( s2->localV_cluster()[0] + s2->localV_cluster()[1] ) / 2;
        if ( (s1!=s2) && (abs(u1-u2) < epsilon) )
          stub_killed[s2] = true;
      }

  // filter out killed stubs
  std::vector<const Stub*> vStubs_filtered;
  for (const Stub* s : vStubs)
    if ( ! stub_killed[s] )
      vStubs_filtered.push_back(s);

  return vStubs_filtered;
}

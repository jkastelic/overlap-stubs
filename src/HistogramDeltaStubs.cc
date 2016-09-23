#include "TMTrackTrigger/TMTrackFinder/interface/Histos.h"
#include "TMTrackTrigger/TMTrackFinder/interface/Stub.h"

#include <TH1F.h>
#include <TH2F.h>

#include <unordered_map>


// Book histograms to do with nearby stubs in a module (which may be produced by delta rays)
void Histos::bookDeltaStubs() {
  TFileDirectory deltaDir = fs_->mkdir("DeltaStubs");
  his_delta_U = deltaDir.make<TH1F>("his_delta_U","number of same-module pairs of stubs;abs(u1-u2);",50,0,1000);
  his_delta_V = deltaDir.make<TH1F>("his_delta_V","number of same-module pairs of stubs;abs(v1-v2);",50,0,1000);
  his_delta_U_hr = deltaDir.make<TH1F>("his_delta_U_hr","number of same-module pairs of stubs;abs(u1-u2);",100,0,2);
  his_delta_V_hr = deltaDir.make<TH1F>("his_delta_V_hr","number of same-module pairs of stubs;abs(v1-v2);",100,0,20);
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
        his_delta_U -> Fill( abs(u1-u2) );
        his_delta_V -> Fill( abs(v1-v2) );
        his_delta_U_hr -> Fill( abs(u1-u2) );
        his_delta_V_hr -> Fill( abs(v1-v2) );
      }
    }

}

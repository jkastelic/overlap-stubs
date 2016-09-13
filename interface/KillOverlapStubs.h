#ifndef __KILL_OVERLAP_STUBS_H__
#define __KILL_OVERLAP_STUBS_H__

#include <vector>

class Settings;
class Stub;
class TP;

class KillOverlapStubs {

  public:
    // constructors
    KillOverlapStubs(const std::vector<const Stub*>& vStubs, const Settings* settings, double pt_cut=3.0, double z0_cut=15.0)
      : settings_(settings), vStubs_(vStubs), pt_cut_(pt_cut), z0_cut_(z0_cut) {}

    // get vStubs without stub-duplicates
    const std::vector<const Stub*> getFiltered() const;

  private:
    // functions that do the filtering
    std::vector<const Stub*> pairFinder        () const;
    std::vector<const Stub*> truePairFinder    () const;

    // helper function for the filtering functions
    bool neighb_modules (const Stub* s1, const Stub* s2) const;
    const TP* commonTP  (const Stub* s1, const Stub* s2) const;

    // data members
    const Settings *settings_;
    const std::vector<const Stub*> vStubs_;
    double pt_cut_, z0_cut_;

};

#endif

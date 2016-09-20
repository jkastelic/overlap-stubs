#ifndef __KILL_OVERLAP_STUBS_H__
#define __KILL_OVERLAP_STUBS_H__

#include <vector>
#include <string>

class Settings;
class Stub;
class TP;

class KillOverlapStubs {

  public:
    // TODO: set the cuts in settings or elsewhere
    KillOverlapStubs(const std::vector<const Stub*>& vStubs, const Settings* settings, double pt_cut=3.0, double z0_cut=15.0)
      : settings_(settings), vStubs_(vStubs), pt_cut_(pt_cut), z0_cut_(z0_cut) {}

    const std::vector<const Stub*> getFiltered(std::string method) const;
    std::vector< std::pair<const Stub*, const Stub*> > getPairs(std::string method) const;

    std::pair<double,double> getTrackParams(const Stub* s1, const Stub* s2) const { return trackParams(s1, s2); }
    const TP* getCommonTP (const Stub* s1, const Stub* s2) const { return commonTP(s1, s2); }
    const TP* getFirstTP(const Stub* s) const { return firstTP(s); }

  private:
    // identify pairs of stubs
    std::vector< std::pair<const Stub*, const Stub*> > pairFinder () const;
    std::vector< std::pair<const Stub*, const Stub*> > truePairFinder() const;

    // helper functions for the above
    bool neighb_modules (const Stub* s1, const Stub* s2) const;
    const TP* commonTP (const Stub* s1, const Stub* s2) const;
    std::pair<double,double> trackParams(const Stub* s1, const Stub* s2) const;
    const TP* firstTP(const Stub* s) const;

    // data members
    const Settings *settings_;
    const std::vector<const Stub*> vStubs_;
    double pt_cut_, z0_cut_;

};

#endif

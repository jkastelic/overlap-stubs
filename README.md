### Killing stubs in overlap regions

I implement a `KillOverlapStubs` class inside the `TMTrackFinder` (taken from `/opt/ppd/cms/users/tomalin/jacob2/`) which, given a vector of stubs, returns a filtered vectors, with the so-called duplicates removed. See [header](interface/KillOverlapStubs.h) and [implementation](src/KillOverlapStubs.cc) for details.


#### Defining duplicates

A true duplicate is effectively defined by the `KillOverlapStubs::truePairFinder()` function: a stub is a duplicate of another if the following conditions are satisfied

- it is from the same layer
- is from a different module
- shares at least one TP with `pt>pt_cut`.

`KillOverlapStubs::pairFinder()` tries to find as many of the true duplicate pairs, but will probably also find some false pairs: that is, pairs where the two stubs do not share any TPs.


#### Histogramming the algorithms

In the `Histos` class, the functions `Histos::bookStubPairs()` and `Histos::fillStubPairs()` will, when implemented, book and fill histograms related to the algorithms in `KillOverlapStubs`. For ease of maintenance, these will be found in [src/HistOverlapStubs.cc](src/HistOverlapStubs.cc).

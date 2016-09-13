import FWCore.ParameterSet.Config as cms

#=== Import default values for all parameters.

from TMTrackTrigger.TMTrackFinder.TMTrackProducer_Defaults_cfi import TMTrackProducer

#=== Change below any parameters you don't want to take their default values.

# If this is True, then an r-z HT will run on stubs assigned to any track candidate produced
# by the r-phi HT.

# TMTrackProducer.HTArraySpecRz.EnableRzHT = cms.bool(True)  

# If this is True, then the seed filter (requiring stubs to lie on line in r-z plane) will be run during r-phi HT.
# N.B. It is not recommended to use both the r-z HT and the seed filter, but only one or the other.

#TMTrackProducer.RZfilterOpts.UseSeedFilter  = cms.bool(True)

# If this is True, a filter will be applied to the stub bend information when filling the r-phi HT,
# with the goal of reducing the number of fake tracks found.

# TMTrackProducer.HTFillingRphi.UseBendFilter = cms.bool(False)


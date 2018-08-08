#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

HcalMipTracks = ldmxcfg.Producer( "HcalMipTracks" , "ldmx::HcalMipTrackProducer" )

# name of collection that contains the hcal hits 
HcalMipTracks.parameters["HcalHitCollectionName"] = "hcalDigis"
HcalMipTracks.parameters["HcalHitPassName"] = "recon"

# Name of track collection to be added to event bus
HcalMipTracks.parameters["HcalMipTrackCollectionName"] = "hcalMipTracks"
 
# minimum PE for a single hit to be considered non-noise
HcalMipTracks.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
HcalMipTracks.parameters["MaximumEnergy"] = 4000.0

# minimum number of clusters to be willing to be considered a plausible track
HcalMipTracks.parameters["MinimumNumClusters"] = 5

# minimum fraction of total number of clusters to continue search
HcalMipTracks.parameters["FractionTotalClusters"] = 0.4

# minimum fraction of clusters currently in log to accept track
#  The higher this fraction is, the slower the algorithm will be
#  Lower this parameter to speed up algorithm (but then you might
#  get single tracks splitting into multiple pieces)
HcalMipTracks.parameters["FractionClustersLeft"] = 0.7

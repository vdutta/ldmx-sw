#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

HcalTracks = ldmxcfg.Producer("HcalTracks", "ldmx::HcalTrackProducer")

# name of collection that contains the hcal hits 
#HcalTracks.parameters["HitCollectionName"] = "hcalDigis"

# hcal specifications
#HcalTracks.parameters["NumHcalLayers"] = 81 
#HcalTracks.parameters["NumHcalStrips"] = 34

# minimum PE for a single hit to be considered non-noise
#HcalTracks.parameters["MinimumPE"] = 15.5

# Search Cone around seed specifications
#  Depth - number of layers to search away from seed
#  Angle - number of strips to search across first layer away from seed
#  Minimum Hits - minimum number of hits needed to be found around seed
#HcalTracks.parameters["SearchConeDepth"] = 3
#HcalTracks.parameters["SearchConeAngle"] = 3
#HcalTracks.parameters["MinConeHits"] = 3

# Width of search in each layer (by number of strips) when extending track
#HcalTracks.parameters["TrackWidth"] = 3

# Minimum number of hits in track to be considered valid track
#HcalTracks.parameters["MinTrackLayHits"] = 40


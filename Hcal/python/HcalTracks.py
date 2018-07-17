#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

HcalTracks = ldmxcfg.Producer( "HcalTracks" , "ldmx::HcalTrackProducer" )

# name of collection that contains the hcal hits 
HcalTracks.parameters["HitCollectionName"] = "hcalDigis"
HcalTracks.parameters["HitPassName"] = "recon"

# hcal specifications
HcalTracks.parameters["NumHcalLayers"] = 81 
HcalTracks.parameters["NumHcalStrips"] = 34

# minimum PE for a single hit to be considered non-noise
HcalTracks.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
HcalTracks.parameters["MaximumEnergy"] = 4000.0

# first layer to search for seed in
HcalTracks.parameters["FirstSeedLayer"] = 1

# Search Cone around seed specifications
#  Depth - number of layers to search away from seed
#  Angle - number of strips to search across first layer away from seed
#  Minimum Hits - minimum number of hits needed to be found around seed
HcalTracks.parameters["SearchConeDepth"] = 3
HcalTracks.parameters["SearchConeAngle"] = 3
HcalTracks.parameters["MinConeHits"] = 3

# Width of search in each layer (by number of strips) when extending track
HcalTracks.parameters["TrackWidth"] = 6

# Minimum number of hits in track to be considered valid track
HcalTracks.parameters["MinTrackLayerHits"] = 20

# Maximum number of tracks to allowed to be found (prevents infinite loop)
HcalTracks.parameters["MaxTrackCount"] = 100

# Name of track collection to be added to event bus
HcalTracks.parameters["HcalTrackCollectionName"] = "HcalTracks"
 

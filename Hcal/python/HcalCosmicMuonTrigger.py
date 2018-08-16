#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

hcalCosmicMuonTrigger = ldmxcfg.Producer( "hcalCosmicMuonTrigger" , "ldmx::HcalCosmicMuonTriggerProducer" )

# name of collection that contains the hcal hits 
hcalCosmicMuonTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
hcalCosmicMuonTrigger.parameters["HcalHitPassName"] = "recon"

# minimum PE for a single hit to be considered non-noise
hcalCosmicMuonTrigger.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
hcalCosmicMuonTrigger.parameters["MaximumEnergy"] = 4000.0

# maximum difference (in strips) between a hit and the centerline of a track
hcalCosmicMuonTrigger.parameters["TrackRadius"] = 4.0

# Minimum fraction of total number of hits to be considered valid track
hcalCosmicMuonTrigger.parameters["MinFractionHit"] = 0.8

# Absolute minimum number of hits to attempt to find track
hcalCosmicMuonTrigger.parameters["AbsoluteMinNumberHits"] = 2

# Name of track collection to be added to event bus
hcalCosmicMuonTrigger.parameters["HcalMipTriggerObjectName"] = "hcalCosmicMuonTrigger"
 

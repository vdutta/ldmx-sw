#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

hcalTargetMuonTrigger = ldmxcfg.Producer( "hcalTargetMuonTrigger" , "ldmx::HcalTargetMuonTriggerProducer" )

# name of collection that contains the hcal hits 
hcalTargetMuonTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
hcalTargetMuonTrigger.parameters["HcalHitPassName"] = "recon"

# minimum PE for a single hit to be considered non-noise
hcalTargetMuonTrigger.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
hcalTargetMuonTrigger.parameters["MaximumEnergy"] = 4000.0

# maximum difference (in strips) between a hit and the centerline of a track
hcalTargetMuonTrigger.parameters["TrackRadius"] = 4.0

# Minimum fraction of total number of hits to be considered valid track
hcalTargetMuonTrigger.parameters["MinFractionHit"] = 0.8

# Absolute minimum number of hits to attempt to find track
hcalTargetMuonTrigger.parameters["AbsoluteMinNumberHits"] = 2

# Name of track collection to be added to event bus
hcalTargetMuonTrigger.parameters["HcalTargetMuonTriggerObjectName"] = "hcalTargetMuonTrigger"
 

#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

HcalMipTrigger = ldmxcfg.Producer( "HcalMipTrigger" , "ldmx::HcalMipTriggerProducer" )

# name of collection that contains the hcal hits 
HcalMipTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
HcalMipTrigger.parameters["HcalHitPassName"] = "recon"

# minimum PE for a single hit to be considered non-noise
HcalMipTrigger.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
HcalMipTrigger.parameters["MaximumEnergy"] = 4000.0

# maximum difference (in strips) between a hit and the centerline of a track
HcalMipTrigger.parameters["TrackRadius"] = 4.0

# Minimum fraction of total number of hits to be considered valid track
HcalMipTrigger.parameters["MinFractionHit"] = 0.8

# Absolute minimum number of hits to attempt to find track
HcalMipTrigger.parameters["AbsoluteMinNumberHitsToLook"] = 2

# Absolute minimum number of hits to consider track valid
HcalMipTrigger.parameters["AbsoluteMinNumberHitsToAccept"] = 2

# Name of track collection to be added to event bus
HcalMipTrigger.parameters["HcalMipTriggerObjectName"] = "hcalMipTrigger"
 

#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

HcalMipTrigger = ldmxcfg.Producer( "HcalMipTrigger" , "ldmx::HcalMipTriggerProducer" )

# name of collection that contains the hcal hits 
HcalMipTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
HcalMipTrigger.parameters["HcalHitPassName"] = "recon"

# hcal specifications
HcalMipTrigger.parameters["NumLayersBackHcal"]   = 81 
HcalMipTrigger.parameters["NumLayersTopHcal"]    = 17 
HcalMipTrigger.parameters["NumLayersBottomHcal"] = 17 
HcalMipTrigger.parameters["NumLayersLeftHcal"]   = 17 
HcalMipTrigger.parameters["NumLayersRightHcal"]  = 17 

# minimum PE for a single hit to be considered non-noise
HcalMipTrigger.parameters["MinimumPE"] = 5.5

# maximum Energy of a group of hits to be considered a mip
HcalMipTrigger.parameters["MaximumEnergy"] = 4000.0

# maximum difference (in strips) between a hit and the centerline of a track
HcalMipTrigger.parameters["TrackRadius"] = 4.0

# Minimum fraction of layers hit in track to be considered valid track
HcalMipTrigger.parameters["MinFractionLayersHit"] = 0.8

# Absolute minimum number of hits to attempt to find track
HcalMipTrigger.parameters["AbsoluteMinNumberHits"] = 3

# Name of track collection to be added to event bus
HcalMipTrigger.parameters["HcalMipTriggerObjectName"] = "hcalMipTrigger"
 

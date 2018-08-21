#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

targetMuonTrigger = ldmxcfg.Producer( "targetMuonTrigger" , "ldmx::MuonTrigger" )

# name of collection that contains the hcal hits 
targetMuonTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
targetMuonTrigger.parameters["HcalHitPassName"] = "recon"

# Name of track collection to be added to event bus
targetMuonTrigger.parameters["TriggerObjectName"] = "targetMuonTrigger"

# Cuts on number of Consecutiveutive layers, strips hit for each section
# if any section passes the cuts, then the whole event passes
targetMuonTrigger.parameters["MinConsecutiveLayersHitBackHcal"] = 60
targetMuonTrigger.parameters["MinConsecutiveStripsHitBackHcal"] = 0

targetMuonTrigger.parameters["MinConsecutiveLayersHitTopHcal"] = 10
targetMuonTrigger.parameters["MinConsecutiveStripsHitTopHcal"] = 0

targetMuonTrigger.parameters["MinConsecutiveLayersHitBottomHcal"] = 10
targetMuonTrigger.parameters["MinConsecutiveStripsHitBottomHcal"] = 0

targetMuonTrigger.parameters["MinConsecutiveLayersHitLeftHcal"] = 10
targetMuonTrigger.parameters["MinConsecutiveStripsHitLeftHcal"] = 0

targetMuonTrigger.parameters["MinConsecutiveLayersHitRightHcal"] = 10
targetMuonTrigger.parameters["MinConsecutiveStripsHitRightHcal"] = 0

#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

cosmicMuonTrigger = ldmxcfg.Producer( "cosmicMuonTrigger" , "ldmx::MuonTrigger" )

# name of collection that contains the hcal hits 
cosmicMuonTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
cosmicMuonTrigger.parameters["HcalHitPassName"] = "recon"

# Name of track collection to be added to event bus
cosmicMuonTrigger.parameters["TriggerObjectName"] = "cosmicMuonTrigger"

# Cuts on number of Consecutiveutive layers, strips hit for each section
# if any section passes the cuts, then the whole event passes
cosmicMuonTrigger.parameters["MinConsecutiveLayersHitBackHcal"] = 20
cosmicMuonTrigger.parameters["MinConsecutiveStripsHitBackHcal"] = 10

cosmicMuonTrigger.parameters["MinConsecutiveLayersHitTopHcal"] = 10
cosmicMuonTrigger.parameters["MinConsecutiveStripsHitTopHcal"] = 0

cosmicMuonTrigger.parameters["MinConsecutiveLayersHitBottomHcal"] = 10
cosmicMuonTrigger.parameters["MinConsecutiveStripsHitBottomHcal"] = 0

cosmicMuonTrigger.parameters["MinConsecutiveLayersHitLeftHcal"] = 10
cosmicMuonTrigger.parameters["MinConsecutiveStripsHitLeftHcal"] = 0

cosmicMuonTrigger.parameters["MinConsecutiveLayersHitRightHcal"] = 10
cosmicMuonTrigger.parameters["MinConsecutiveStripsHitRightHcal"] = 0

#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

hcalMuonTrigger = ldmxcfg.Producer( "hcalMuonTrigger" , "ldmx::MuonTrigger" )

# name of collection that contains the hcal hits 
hcalMuonTrigger.parameters["HcalHitCollectionName"] = "hcalDigis"
hcalMuonTrigger.parameters["HcalHitPassName"] = "recon"

# Name of track collection to be added to event bus
hcalMuonTrigger.parameters["HcalMuonTriggerObjectName"] = "hcalMuonTrigger"

# Cuts on number of consecutive layers, strips hit for each section
# if any section passes the cuts, then the whole event passes
hcalMuonTrigger.parameters["MinConsecLayersHitBackHcal"] = 20
hcalMuonTrigger.parameters["MinConsecStripsHitBackHcal"] = 10

hcalMuonTrigger.parameters["MinConsecLayersHitTopHcal"] = 10
hcalMuonTrigger.parameters["MinConsecStripsHitTopHcal"] = 0

hcalMuonTrigger.parameters["MinConsecLayersHitBottomHcal"] = 10
hcalMuonTrigger.parameters["MinConsecStripsHitBottomHcal"] = 0

hcalMuonTrigger.parameters["MinConsecLayersHitLeftHcal"] = 10
hcalMuonTrigger.parameters["MinConsecStripsHitLeftHcal"] = 0

hcalMuonTrigger.parameters["MinConsecLayersHitRightHcal"] = 10
hcalMuonTrigger.parameters["MinConsecStripsHitRightHcal"] = 0

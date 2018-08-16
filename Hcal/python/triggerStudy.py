#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

triggerStudy = ldmxcfg.Analyzer("triggerStudy", "ldmx::TriggerAnalyzer")

# name of hcal  trigger to analyze
triggerStudy.parameters[ "HcalTriggerObjectName" ] = "hcalMipTrigger"

# name of pass that create hcal  trigger
triggerStudy.parameters[ "HcalTriggerPassName" ] = "recon"


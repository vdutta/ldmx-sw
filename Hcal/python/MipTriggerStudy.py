#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

MipTriggerStudy = ldmxcfg.Analyzer("MipTriggerStudy", "ldmx::MipTriggerAnalyzer")

# name of hcal mip trigger to analyze
MipTriggerStudy.parameters[ "HcalMipTriggerObjectName" ] = "hcalMipTrigger"

# name of pass that create hcal mip trigger
MipTriggerStudy.parameters[ "HcalMipTriggerPassName" ] = "recon"


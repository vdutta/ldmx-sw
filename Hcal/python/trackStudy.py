#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

trackStudy = ldmxcfg.Analyzer("trackStudy", "ldmx::HcalTrackAnalyzer")

# name of collection that contains the hcal tracks 
trackStudy.parameters["HcalMipTracksCollectionName"] = "hcalMipTracks"

# name of pass that create the hcal mip tracks
trackStudy.parameters["HcalMipTracksPassName"] = "recon"

# minimum number of non-noise HcalHits for a track to be findable
trackStudy.parameters["MinHcalHits"] = 3

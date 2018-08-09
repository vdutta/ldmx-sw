#!/usr/bin/python

# we need the ldmx configuration package to construct the object
from LDMX.Framework import ldmxcfg

TrackStudy = ldmxcfg.Analyzer("TrackStudy", "ldmx::HcalTrackAnalyzer")

# name of collection that contains the hcal tracks 
TrackStudy.parameters["HcalMipTracksCollectionName"] = "hcalMipTracks"

# name of pass that create the hcal mip tracks
TrackStudy.parameters["HcalMipTracksPassName"] = "recon"


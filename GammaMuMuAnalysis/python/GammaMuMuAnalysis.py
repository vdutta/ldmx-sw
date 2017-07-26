import sys

from LDMX.Framework import ldmxcfg

p=ldmxcfg.Process("recon")

p.libraries.append("libEventProc.so")
p.libraries.append("libGammaMuMuAnalysis.so")

analyzer = ldmxcfg.Analyzer("analyzer","ldmx::kinematicPlots")
analyzer.parameters["simParticleCollection"] = "SimParticles"

p.histogramFile = "histo.root"

p.sequence=[analyzer]

p.inputFiles=[sys.argv[1]]
p.outputFiles=[sys.argv[2]]

p.printMe()

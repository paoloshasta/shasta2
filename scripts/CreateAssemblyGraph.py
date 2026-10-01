#!/usr/bin/python3

import shasta2

shasta2.openPerformanceLog("Python-performance.log")
options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.accessAnchorGraph("")
assembler.createHomopolymerModel(options.homopolymerModelName)
assembler.createAssemblyGraph(options, False)

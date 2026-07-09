#!/usr/bin/python3

import shasta2

options = shasta2.Options()

assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph("J", options)
assemblyGraph.removeZeroLengthSegments()

assemblyGraph.write("Test")




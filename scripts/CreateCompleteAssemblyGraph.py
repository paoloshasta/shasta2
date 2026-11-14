#!/usr/bin/python3

import shasta2

options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()
assembler.accessCompleteAnchorGraph()
assembler.createCompleteAssemblyGraph(options)

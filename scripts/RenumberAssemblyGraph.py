#!/usr/bin/python3

import shasta2

# Get the argument.
import argparse
parser = argparse.ArgumentParser(description = "Renumber assembly grapg segments.")
parser.add_argument("inputStage", type=str, help="Input assembly stage.")
parser.add_argument("outputStage", type=str, help="Output assembly stage.")
arguments = parser.parse_args()

# Access what we need.
options = shasta2.Options()
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.inputStage, options)
assemblyGraph.renumber()
assemblyGraph.write(arguments.outputStage)





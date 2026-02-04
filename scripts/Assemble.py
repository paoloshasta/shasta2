#!/usr/bin/python3

import shasta2

# Get the argument.
import argparse
parser = argparse.ArgumentParser(description = "Run sequence assembly for a given assembly stage.")
parser.add_argument("stage", type=str, help="Assembly stage.")
arguments = parser.parse_args()

shasta2.openPerformanceLog("Python-performance.log")

# Get the options from shasta2.conf.
options = shasta2.Options()

# Create the Assembler and access what we need.
assembler = shasta2.Assembler()
assembler.accessMarkers()
assembler.accessMarkerKmers()
assembler.accessAnchors()
assembler.accessJourneys()

# Load the specified assembly stage and assemble sequence.
assemblyGraph = assembler.getAssemblyGraph(arguments.stage, options)
assemblyGraph.assembleAll()

# Write it out.
assembledName = arguments.stage + "-Assembled"
assemblyGraph.write(assembledName)
assemblyGraph.writeFasta(assembledName)





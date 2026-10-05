#!/usr/bin/python3

import shasta2

# Get the arguments.
import argparse
parser = argparse.ArgumentParser(description = 'Turn the AssemblyGraph into a "pseudo" Assembly graph')
parser.add_argument("inputStage", type=str, help="Input assembly stage. Must be double-stranded, strand-symmetric.")
parser.add_argument("outputStage", type=str, help="Output assembly stage. Double-stranded, strand-symmetric.")
arguments = parser.parse_args()

shasta2.openPerformanceLog("Python-performance.log")

# Get the options from shasta2.conf.
options = shasta2.Options()

# Create the Assembler and access what we need.
assembler = shasta2.Assembler()
assembler.accessAnchors()
assembler.accessJourneys()

assemblyGraph = assembler.getAssemblyGraph(arguments.inputStage, options)
assemblyGraph.makePseudo()
assemblyGraph.compress()
assemblyGraph.removeIsolatedVertices()

# Write it out.
assemblyGraph.write(arguments.outputStage)
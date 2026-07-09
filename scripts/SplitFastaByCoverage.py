#!/usr/bin/python3

import argparse
import csv
import os

# Get the coverage threshold.
parser = argparse.ArgumentParser(description = "Split Assembly-Final.fasta into low and high coverage portions.")
parser.add_argument("coverageThreshold", type=int, help="Coverage threshold.")
arguments = parser.parse_args()
coverageThreshold = arguments.coverageThreshold


# Create files containing lists of the low coverage and high coverage segments.
csvFile = open("Assembly-Final.csv")
csvReader = csv.DictReader(csvFile)
lowCoverageList = open("LowCoverage.txt", "w")
highCoverageList = open("HighCoverage.txt", "w")
for row in csvReader:
    segmentName = row["Segment"]
    coverage = int(row["Average coverage"])
    if coverage <= coverageThreshold:
        lowCoverageList.write("{}\n".format(segmentName))
    else:
        highCoverageList.write("{}\n".format(segmentName))
lowCoverageList.close()
highCoverageList.close()

# Use seqkit to extract the segments in these lists.
os.system("seqkit grep -f LowCoverage.txt Assembly-Final.fasta -o LowCoverage.fasta")
os.system("seqkit grep -f HighCoverage.txt Assembly-Final.fasta -o HighCoverage.fasta")

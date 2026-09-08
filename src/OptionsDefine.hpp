// General options.
SHASTA2_OPTION_DEFINE(
    string, command, "--command", "assemble",
    "Shasta2 command.")

SHASTA2_OPTION_DEFINE(
    string, assemblyDirectory, "--output", "ShastaRun",
    "Assembly directory for output.")

SHASTA2_OPTION_DEFINE(
    string, memoryMode, "--memory-mode", "anonymous",
    "Memory mode. Specify anonymous or filesystem.")

SHASTA2_OPTION_DEFINE(
    string, memoryBacking, "--memory-backing", "4K",
    "Memory backing. Specify disk, 4K, or 2M.")

SHASTA2_BOOL_OPTION_DEFINE(
    keepBinaryData, "--keep-binary-data", false,
    "Keep binary data. Only effective if --memory-mode filesystem.")

SHASTA2_OPTION_DEFINE(
    uint64_t, threadCount, "--threads", 0,
    "Number of threads (0 to use hardware_concurrency).")

SHASTA2_OPTION_DEFINE(
    string, externalAnchorsName, "--external-anchors-name", "",
    "External anchors name. Must be an absolute or relarive path relative base name "
    "for the four external anchors files.")

SHASTA2_OPTION_DEFINE(
    string, externalAnchorGraphName, "--external-anchor-graph-name", "",
    "External anchor graph name. Only in effect if --external-anchors-name is also used.")



// Http server.
SHASTA2_OPTION_DEFINE(
    string, exploreAccess, "--explore-access", "user",
        "Specifies accessibility of Shasta2 http server.\n"
        "Specify user, local, or unrestricted.")

SHASTA2_OPTION_DEFINE(
    uint16_t, port, "--port", 17100,
    "Port number for http server.")



// Reads.
SHASTA2_OPTION_DEFINE(
    uint64_t, minReadLength, "--min-read-length", 0,
    "Read length cutoff.")

SHASTA2_VECTOR_OPTION_DEFINE(
    string, inputFileNames, "--input", inputFileNamesDefault,
    "Input fasta or fastq files (uncompressed).")



// Markers.
SHASTA2_OPTION_DEFINE(
    uint64_t, k, "--k", 60,
    "Marker length.")

SHASTA2_OPTION_DEFINE(
    double, markerDensity, "--marker-density", 0.05,
    "Marker density.")

SHASTA2_OPTION_DEFINE(
    double, maxMarkerErrorRate, "--max-marker-error-rate", 0.5,
    "Maximum marker error rate. Reads with a higher marker error rate are not used.")



// Anchors.
SHASTA2_OPTION_DEFINE(
    uint64_t, minAnchorCoverage, "--min-anchor-coverage", 10,
    "Minimum anchor coverage.")

SHASTA2_OPTION_DEFINE(
    uint64_t, maxAnchorCoverage, "--max-anchor-coverage", 60,
    "Maximum anchor coverage.")

SHASTA2_VECTOR_OPTION_DEFINE(
uint64_t, maxAnchorRepeatLength, "--max-anchor-repeat-length", maxAnchorRepeatLengthDefault,
R"zzz(
Maximum number of copies of repeats of period 1, 2, 3,... allowed in an anchor sequence.
An anchor is not generated if its sequence would contain an exact repeat
consisting of n copies of a unit of length (period) p, if
n > maxAnchorRepeatLength[p-1].
So for example:
- maxAnchorRepeatLength[0] is the maximum allowed length of
  a homopolymer run.
- maxAnchorRepeatLength[1] is the maximum allowed number of
  copies of a repeat with period 2 (e. g. ATATAT).
  Note this is the number of copies, not number of bases.
  So if maxAnchorRepeatLength[1] is 3, the anchor is not
  generated if it contains a 2-base repeat with more than 3 copies
  (a total 6 bases).
)zzz")

SHASTA2_VECTOR_OPTION_DEFINE(
uint64_t, minAnchorDistinctSubkmerCount, "--min-anchor-distinct-subkmer-count", minAnchorDistinctSubkmerCountDefault,
R"zzz(
Minimum number of distinct sub-k-mers of length 1, 2, 3,... required in an anchor sequence.
An anchor is not generated if its sequence would contain less than the specified 
number of distinct sub-k-mers, which indicates low sequence complexity.
For example, --min-anchor-distinct-subkmer-count 4, 12, 24 means
that the k-mer sequence must contain all 4 bases,at least 12 of the 16
possible sub-2-mers, and at least 24 of the 64 possible sub-3-mers.
)zzz")



// Anchor graph.
SHASTA2_OPTION_DEFINE(
    uint64_t, minAnchorGraphEdgeCoverage, "--min-anchor-graph-edge-coverage", 6,
    "Minimum anchor graph edge coverage.")

SHASTA2_OPTION_DEFINE(
    double, minAnchorGraphEdgeCoverageFraction, "--min-anchor-graph-edge-coverage-fraction", 0.,
    "Minimum anchor graph edge coverage fraction.")

SHASTA2_OPTION_DEFINE(
    uint64_t, transitiveReductionMaxEdgeCoverage, "--transitive-reduction-max-edge-coverage", 10,
    "Maximum coverage of an AnchorGraph edge subject to removal during transitive reduction.")

SHASTA2_OPTION_DEFINE(
    uint64_t, transitiveReductionMaxDistance, "--transitive-reduction-max-distance", 10,
    "Maximum distance for transitive reduction of the AnchorGraph.")



// Assembly graph.
SHASTA2_BOOL_OPTION_DEFINE(
    writeIntermediateAssemblyStages, "--write-intermediate-assembly-stages", false,
    "Write intermediate assembly stages.")

SHASTA2_BOOL_OPTION_DEFINE(
    writeAssemblyDetails, "--write-assembly-details", false,
    "Write assembly details in csv format.")

SHASTA2_OPTION_DEFINE(
    uint64_t, bubbleCleanupMaxBubbleLength, "--bubble-cleanup-max-bubble-length", 10000,
    "Maximum bubble length for bubble cleanup.")

SHASTA2_OPTION_DEFINE(
    uint64_t, bubbleCleanupMinCommonCount, "--bubble-cleanup-min-common-count", 6,
    "Minimum number of common oriented reads for bubble cleanup.")

SHASTA2_OPTION_DEFINE(
    uint64_t, findSuperbubblesMaxDistance, "--find-superbubbles-max-distance", 10,
    "Maximum distance (number of BFS edges) when finding superbubbles.")

SHASTA2_OPTION_DEFINE(
    uint64_t, phasingDistance, "--phasing-distance", 12,
    "Maximum bubble-bubble distance (number of edges) for phasing.")

SHASTA2_OPTION_DEFINE(
    double, detangleEpsilon, "--detangle-epsilon", 0.05,
    "Epsilon value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    double, detangleMaxLogP, "--detangle-maxLogP", 30.,
    "MaxLogP value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    double, detangleMinLogPDelta, "--detangle-minLogPDelta", 10.,
    "MinLogPDelta value for likelihood ratio detangling.")

SHASTA2_OPTION_DEFINE(
    uint64_t ,representativeRegionStepCount, "--representative-region-step-count", 10,
    "Number of steps for the representative region of a segment "
    "used in detangling and read following.")

SHASTA2_OPTION_DEFINE(
    uint64_t, pruneLength, "--prune-length", 50000,
    "Maximum leaf length for pruning.")

SHASTA2_OPTION_DEFINE(
    uint64_t, strictPruneLength, "--strict-prune-length", 10000,
    "Maximum leaf length for strict pruning.")


// Read following.
SHASTA2_OPTION_DEFINE(
uint64_t, readFollowingMinCommonCount, "--read-following-min-common-count", 2,
"Minimum number of common oriented reads for read following.")

SHASTA2_OPTION_DEFINE(
	double, readFollowingMinCorrectedJaccard, "--read-following-min-corrected-jaccard", 0.,
	"Minimum corrected Jaccard for read following.")

SHASTA2_OPTION_DEFINE(
	uint32_t, readFollowingPruneLength, "--read-following-prune-length", 100000,
	"Prune length for read following.")

SHASTA2_OPTION_DEFINE(
	uint64_t, readFollowingSegmentLengthThreshold, "--read-following-segment-length-threshold", 500000,
	"Segment length threshold for read following.")

SHASTA2_OPTION_DEFINE(
    uint64_t, readFollowingPathCount, "--read-following-path-count", 100,
    "Number of random paths per vertex for read following.")

SHASTA2_OPTION_DEFINE(
    uint64_t, readFollowingPathCountThreshold1, "--read-following-path-count-threshold1", 30,
    "First path count threshold for read following.")

SHASTA2_OPTION_DEFINE(
    uint64_t, readFollowingPathCountThreshold2, "--read-following-path-count-threshold2", 60,
    "Second path count threshold for read following.")


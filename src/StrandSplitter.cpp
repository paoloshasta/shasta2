// Shasta2.
#include "StrandSplitter.hpp"
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "html.hpp"
using namespace shasta2;

// Standard library.



// This takes as input a self-complementary tangle
// in an AssemblyGraph and attempts to split the strands.
// The last two arguments are only used for debug output.
StrandSplitter::StrandSplitter(
    AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraphBaseClass::vertex_descriptor>& tangleVertices,
    uint64_t tangleId,
    const string& debugOutputBaseName) :
    assemblyGraph(assemblyGraph),
    tangleId(tangleId),
    debugOutputBaseName(debugOutputBaseName),
    tangle(assemblyGraph, tangleVertices)
{
    writeInitialDebugOutput();
    SHASTA2_ASSERT(tangle.isSelfComplementary());
    gatherSegments();
    writeSegmentPairs();
    findReadOccurrences();
    createGraph();
    separateStrands();
}



StrandSplitter::~StrandSplitter()
{
    if(debug) {
        writeHtmlEnd(html);
    }
}



void StrandSplitter::writeInitialDebugOutput()
{
    if(debug) {
        cout << "StrandSplitter begins for tangle " << tangleId << endl;
        html.open(debugOutputBaseName + "-StrandSplitter-Tangle-" + to_string(tangleId) + ".html");
        writeHtmlBegin(html, "Tangle " + to_string(tangleId));
        html << "<h1>Self-complementary tangle " << tangleId << "</h1>";
        tangle.writeHtml(html);
    }

}


void StrandSplitter::gatherSegments()
{
    // Fill in allTangleSegments.
    const auto inserter = back_inserter(allTangleSegments);
    std::ranges::copy(tangle.tangleEdges, inserter);
    std::ranges::copy(tangle.entrances, inserter);
    std::ranges::copy(tangle.exits, inserter);
    sort(allTangleSegments.begin(), allTangleSegments.end(), assemblyGraph.orderById);

    // Fill in segmentPairs.
    for(const Segment segment: allTangleSegments) {
        if(assemblyGraph[segment].lengthWeightedAverageCoverage() > maxCoverage) {
            continue;
        }
        const Segment segmentRc = assemblyGraph[segment].eRc;
        SHASTA2_ASSERT(segmentRc != segment);
        if(id(segment) < id(segmentRc)) {
            segmentPairs.push_back({segment, segmentRc});
        }
    }

    // Fill in the segments vector.
    for(const auto&[segment, segmentRc]: segmentPairs) {
        segments.push_back(segment);
        segments.push_back(segmentRc);
    }

}



uint64_t StrandSplitter::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandSplitter::writeSegmentPairs()
{
    if(debug) {
        html << "<h2>Pairs of reverse complemented segments</h2>"
            "This includes entrances, exits, and internal tangle segments "
            " with coverage up to " << maxCoverage <<
            ". These are considered reliable single copy segments "
            "that can be used for strand separation."
            "<br><table>"
            "<tr><th>Pair index<th>Index0<th>Index1<th>Segment0<th>Segment1";
        for(uint64_t segmentPairIndex=0; segmentPairIndex<segmentPairs.size(); segmentPairIndex++) {
            const auto&[segment, segmentRc] = segmentPairs[segmentPairIndex];
            html << "<tr>"
                "<td class=centered>" << segmentPairIndex <<
                "<td class=centered>" << 2*segmentPairIndex <<
                "<td class=centered>" << 2*segmentPairIndex+1 <<
                "<td class=centered>" << id(segment) <<
                "<td class=centered>" << id(segmentRc);
        }
        html << "</table>";
    }

}



void StrandSplitter::findReadOccurrences()
{
    for(uint64_t segmentPairIndex=0; segmentPairIndex<segmentPairs.size(); segmentPairIndex++) {
        const auto&[segment, ignore] = segmentPairs[segmentPairIndex];
        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();
                readOccurrenceMap[readId].emplace_back(ReadOccurrence(segmentPairIndex, strand));
            }
        }
    }

    // Deduplicate and count the Occurrences for each ReadId.
    vector<uint64_t> count;
    for(auto&[readId, occurrences]: readOccurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }

}



void StrandSplitter::Graph::addToEdge(
    uint64_t segmentIndex0,
    uint64_t segmentIndex1,
    uint64_t frequency
    )
{
    auto[e, edgeExists] = boost::edge(segmentIndex0, segmentIndex1, *this);
    if(edgeExists) {
        (*this)[e].frequency += frequency;
    } else {
        boost::add_edge(segmentIndex0, segmentIndex1, Edge(frequency), *this);
    }
}



void StrandSplitter::createGraph()
{
    for(const Segment segment: segments) {
        boost::add_vertex(segment, graph);
    }
    for(const auto&[ignore, occurrences]: readOccurrenceMap) {
        for(uint64_t i0=0; i0<occurrences.size(); i0++) {
            const ReadOccurrence& occurrence0 = occurrences[i0];
            const uint64_t segmentPairIndex0 = occurrence0.segmentPairIndex;
            const uint64_t strand0 = occurrence0.strand;
            for(uint64_t i1=i0+1; i1<occurrences.size(); i1++) {
                const ReadOccurrence& occurrence1 = occurrences[i1];
                const uint64_t segmentPairIndex1 = occurrence1.segmentPairIndex;
                const uint64_t strand1 = occurrence1.strand;
                const uint64_t frequency = occurrence0.frequency * occurrence1.frequency;
                if(strand0 == strand1) {
                    uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                    uint64_t segmentIndex1 = 2 * segmentPairIndex1;
                    graph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                    ++segmentIndex0;
                    ++segmentIndex1;
                    graph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                } else {
                    uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                    uint64_t segmentIndex1 = 2 * segmentPairIndex1 + 1;
                    graph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                    ++segmentIndex0;
                    --segmentIndex1;
                    graph.addToEdge(segmentIndex0, segmentIndex1, frequency);
                }
            }
        }
    }

    graph.findEdgePairs();

}



void StrandSplitter::Graph::findEdgePairs()
{
    Graph& graph = *this;

    std::set<edge_descriptor> edgesFound;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(edgesFound.contains(e)) {
            continue;
        }
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const vertex_descriptor v0Rc = v0 ^ 1;
        const vertex_descriptor v1Rc = v1 ^ 1;
        auto[eRc, edgeExists] = boost::edge(v0Rc, v1Rc, graph);
        SHASTA2_ASSERT(edgeExists);
        SHASTA2_ASSERT(graph[eRc].frequency == graph[e].frequency);

        edgesFound.insert(e);
        edgesFound.insert(eRc);

        edgePairs.emplace_back(EdgePair({e, eRc, graph[e].frequency}));
    }
    sort(edgePairs.begin(), edgePairs.end());

}



void StrandSplitter::separateStrands()
{
    // Do strand separation by adding edges in order of decreasing frequency.
    DisjointSets disjointSets(segments.size());
    uint64_t crossStrandEdgeCount = 0;
    for(const auto& edgePair: graph.edgePairs) {
        const Graph::edge_descriptor eA = edgePair.e;
        const Graph::edge_descriptor eB = edgePair.eRc;
        const uint64_t v0A = source(eA, graph);
        const uint64_t v1A = target(eA, graph);
        const uint64_t v0B = source(eB, graph);
        const uint64_t v1B = target(eB, graph);
        const uint64_t v0ARc = v0A ^ 1;
        const uint64_t v1ARc = v1A ^ 1;
        const uint64_t v0BRc = v0B ^ 1;
        const uint64_t v1BRc = v1B ^ 1;
        const bool strandViolationA = (disjointSets.findSet(v1A) == disjointSets.findSet(v0ARc));
        const bool strandViolationB = (disjointSets.findSet(v1B) == disjointSets.findSet(v0BRc));
        const bool strandViolationARc = (disjointSets.findSet(v0A) == disjointSets.findSet(v1ARc));
        const bool strandViolationBRc = (disjointSets.findSet(v0B) == disjointSets.findSet(v1BRc));
        const bool strandViolation = strandViolationA;
        SHASTA2_ASSERT(strandViolationB == strandViolation);
        SHASTA2_ASSERT(strandViolationARc == strandViolation);
        SHASTA2_ASSERT(strandViolationBRc == strandViolation);
        if(strandViolation) {
            crossStrandEdgeCount += 2;
            graph[eA].isCrossStrandEdge = true;
            graph[eB].isCrossStrandEdge = true;
        } else {
            disjointSets.unionSet(v0A, v1A);
            disjointSets.unionSet(v0B, v1B);
        }
    }

    vector< vector<uint64_t> > components;
    disjointSets.gatherComponents(1, components);

    // Store the component of each vertex.
    for(uint64_t componentId=0; componentId<components.size(); componentId++) {
        const vector<uint64_t>& component = components[componentId];
        for(const uint64_t v: component) {
            graph[v].component = componentId;
        }
    }


    if(debug) {
        ofstream dot(debugOutputBaseName + "-StrandSplitter-Tangle-" + to_string(tangleId) + ".dot");
        dot << "graph splitSelfComplementaryTangle {\n";
        BGL_FORALL_VERTICES(segmentIndex, graph, Graph) {
            const string color = randomHslColor(graph[segmentIndex].component, 0.75, 0.5);
            dot << id(graph[segmentIndex].segment) <<
                " [style=filled fillcolor=\"" << color << "\"]"
                "\n";
        }
        BGL_FORALL_EDGES(e, graph, Graph) {
            const uint64_t segmentIndex0 = source(e, graph);
            const uint64_t segmentIndex1 = target(e, graph);
            dot << id(graph[segmentIndex0].segment) << "--";
            dot << id(graph[segmentIndex1].segment) <<
                "[tooltip=\"" << graph[e].frequency << "\""
                " penwidth=\"" << std::log10(double(graph[e].frequency)) << "\"";
            if(graph[e].isCrossStrandEdge) {
                dot << " color=red";
            }
            dot << "];\n";
        }
        dot << "}\n";


        for(uint64_t componentId=0; componentId<components.size(); componentId++) {
            html << "<h2>Component " << componentId << " segments</h2>";
            const vector<uint64_t>& component = components[componentId];
            for(uint64_t i=0; i<component.size(); i++) {
                if(i != 0) {
                    html << ",<wbr>";
                }
                html << id(graph[component[i]].segment);
            }
        }
    }

}


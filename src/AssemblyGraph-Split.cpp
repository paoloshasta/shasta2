// Shasta2.
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "html.hpp"
#include "performanceLog.hpp"
#include "StrandSplitter.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Standard library.
#include "fstream.hpp"



void AssemblyGraph::splitSelfComplementaryTangles(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t lengthThreshold = 100000;

    performanceLog << timestamp << "AssemblyGraph::splitSelfComplementaryTangles begins: " <<
        debugOutputBaseName << endl;

    // Create the superbubbles as tangles consisting of short Segments.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(lengthThreshold, tangles, tangleRc);

    // Extract the self-complementary ones.
    vector< vector<vertex_descriptor> > selfComplementaryTangles;
    vector<uint64_t> selfComplementaryTangleRc;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleRc[tangleId] == tangleId) {
            selfComplementaryTangleRc.push_back(selfComplementaryTangles.size());
            selfComplementaryTangles.emplace_back(tangles[tangleId]);
        }
    }

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    writeTangles(selfComplementaryTangles, selfComplementaryTangleRc,
        debugOutputBaseName + "-SelfComplementary-Tangles-Bandage.csv");

    for(uint64_t tangleId=0; tangleId<selfComplementaryTangles.size(); tangleId++) {
#if 0
        splitSelfComplementaryTangle(
            tangleId,
            selfComplementaryTangles[tangleId],
            debugOutputBaseName);
#else
        StrandSplitter strandSplitter(
            *this,
            selfComplementaryTangles[tangleId],
            tangleId,
            debugOutputBaseName);
#endif
    }

    performanceLog << timestamp << "AssemblyGraph::splitSelfComplementaryTangles ends: " <<
        debugOutputBaseName << endl;
}



void AssemblyGraph::splitSelfComplementaryTangle(
    uint64_t tangleId,
    const vector<vertex_descriptor>& tangleVertices,
    const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double maxCoverage = 16.;

    // Create the Tangle.
    AssemblyGraph& assemblyGraph = *this;
    const Tangle tangle(assemblyGraph, tangleVertices);

    // Initial debug output.
    const bool debug = true;
    ofstream html;
    if(debug) {
        cout << "AssemblyGraph::splitSelfComplementaryTangle begins for tangle " << tangleId << endl;
        html.open(debugOutputBaseName + "-Tangle-" + to_string(tangleId) + ".html");
        writeHtmlBegin(html, "Tangle " + to_string(tangleId));
        html << "<h1>Self-complementary tangle " << tangleId << "</h1>";
        tangle.writeHtml(html);
    }
    SHASTA2_ASSERT(tangle.isSelfComplementary());



    // Gather all Tangle segments (internal, entrances, exits).
    vector<Segment> allTangleEdges = tangle.tangleEdges;
    std::ranges::copy(tangle.entrances, back_inserter(allTangleEdges));
    std::ranges::copy(tangle.exits, back_inserter(allTangleEdges));

    // Gather pairs of reverse complemented segments.
    vector< pair<Segment, Segment> > segmentPairs;
    for(const Segment segment: allTangleEdges) {
        if(assemblyGraph[segment].lengthWeightedAverageCoverage() > maxCoverage) {
            continue;
        }
        const Segment segmentRc = assemblyGraph[segment].eRc;
        SHASTA2_ASSERT(segmentRc != segment);
        if(id(segment) < id(segmentRc)) {
            segmentPairs.push_back({segment, segmentRc});
        }
    }
    if(debug) {
        html << "<h2>Pairs of reverse complemented segments</h2>"
            "<table>"
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

    // Gather segments in the same order.
    vector<Segment> segments;
    for(const auto&[segment, segmentRc]: segmentPairs) {
        segments.push_back(segment);
        segments.push_back(segmentRc);
    }



    // Gather occurrences of oriented reads in the first Segment of each pair.
    class Occurrence {
    public:
        uint64_t segmentPairIndex = invalid<uint64_t>;
        Strand strand = invalid<Strand>;
        uint64_t frequency = 0;
        Occurrence() {}
        Occurrence(uint64_t segmentPairIndex, Strand strand) :
            segmentPairIndex(segmentPairIndex), strand(strand) {}
        bool operator==(const Occurrence& that) const
        {
            return tie(segmentPairIndex, strand) == tie(that.segmentPairIndex, that.strand);
        }
        bool operator<(const Occurrence& that) const
        {
            return tie(segmentPairIndex, strand) < tie(that.segmentPairIndex, that.strand);
        }
    };
    std::map<ReadId, vector<Occurrence> > occurrenceMap;
    for(uint64_t segmentPairIndex=0; segmentPairIndex<segmentPairs.size(); segmentPairIndex++) {
        const auto&[segment, ignore] = segmentPairs[segmentPairIndex];
        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();
                occurrenceMap[readId].emplace_back(Occurrence(segmentPairIndex, strand));
            }
        }
    }

    // Deduplicate and count the Occurrences for each ReadId.
    vector<uint64_t> count;
    for(auto&[readId, occurrences]: occurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }



    // Now for each ReadId look at pairs of Segments and relative orientations.
    // For each segmentPairIndex i, there are two segmentIndexes 2*i and 2*i+1.
    if(debug) {
        ofstream csv(debugOutputBaseName + "-splitSelfComplementaryTangle.csv");
        for(const auto&[ignore, occurrences]: occurrenceMap) {
            for(uint64_t i0=0; i0<occurrences.size(); i0++) {
                const Occurrence& occurrence0 = occurrences[i0];
                const uint64_t segmentPairIndex0 = occurrence0.segmentPairIndex;
                const uint64_t strand0 = occurrence0.strand;
                for(uint64_t i1=i0+1; i1<occurrences.size(); i1++) {
                    const Occurrence& occurrence1 = occurrences[i1];
                    const uint64_t segmentPairIndex1 = occurrence1.segmentPairIndex;
                    const uint64_t strand1 = occurrence1.strand;
                    const uint64_t frequency = occurrence0.frequency * occurrence1.frequency;
                    if(strand0 == strand1) {
                        uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                        uint64_t segmentIndex1 = 2 * segmentPairIndex1;
                        csv << segmentIndex0 << ",";
                        csv << segmentIndex1 << ",";
                        csv << id(segments[segmentIndex0]) << ",";
                        csv << id(segments[segmentIndex1]) << ",";
                        csv << frequency << endl;
                        ++segmentIndex0;
                        ++segmentIndex1;
                        csv << segmentIndex0 << ",";
                        csv << segmentIndex1 << ",";
                        csv << id(segments[segmentIndex0]) << ",";
                        csv << id(segments[segmentIndex1]) << ",";
                        csv << frequency << endl;
                    } else {
                        uint64_t segmentIndex0 = 2 * segmentPairIndex0;
                        uint64_t segmentIndex1 = 2 * segmentPairIndex1 + 1;
                        csv << segmentIndex0 << ",";
                        csv << segmentIndex1 << ",";
                        csv << id(segments[segmentIndex0]) << ",";
                        csv << id(segments[segmentIndex1]) << ",";
                        csv << frequency << endl;
                        ++segmentIndex0;
                        --segmentIndex1;
                        csv << segmentIndex0 << ",";
                        csv << segmentIndex1 << ",";
                        csv << id(segments[segmentIndex0]) << ",";
                        csv << id(segments[segmentIndex1]) << ",";
                        csv << frequency << endl;
                    }
                }
            }
        }
    }



    // Now create a graph with a vertex for each segment.
    class Vertex {
    public:
        Segment segment;
        uint64_t component = invalid<uint64_t>;
        Vertex(Segment segment = assemblyGraphNullEdge) : segment(segment) {}
    };
    class Edge {
    public:
        uint64_t frequency;
        bool isCrossStrandEdge = false;
        Edge(uint64_t frequency) : frequency(frequency) {}
    };
    using GraphBaseClass = boost::adjacency_list<
        boost::setS,
        boost::vecS,
        boost::undirectedS,
        Vertex,
        Edge>;
    class Graph: public GraphBaseClass {
    public:
        void addToEdge(
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
    };
    Graph graph;
    for(const Segment segment: segments) {
        boost::add_vertex(segment, graph);
    }
    for(const auto&[ignore, occurrences]: occurrenceMap) {
        for(uint64_t i0=0; i0<occurrences.size(); i0++) {
            const Occurrence& occurrence0 = occurrences[i0];
            const uint64_t segmentPairIndex0 = occurrence0.segmentPairIndex;
            const uint64_t strand0 = occurrence0.strand;
            for(uint64_t i1=i0+1; i1<occurrences.size(); i1++) {
                const Occurrence& occurrence1 = occurrences[i1];
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



    // Gather pairs of reverse complemented edges.
    // Sort them by decreasing frequency.
    class EdgePair {
    public:
        Graph::edge_descriptor e;
        Graph::edge_descriptor eRc;
        uint64_t frequency;
        bool operator<(const EdgePair& that) const
        {
            return frequency > that.frequency;
        }
    };
    vector<EdgePair> edgePairs;
    std::set<Graph::edge_descriptor> edgesFound;
    BGL_FORALL_EDGES(e, graph, Graph) {
        if(edgesFound.contains(e)) {
            continue;
        }
        const Graph::vertex_descriptor v0 = source(e, graph);
        const Graph::vertex_descriptor v1 = target(e, graph);
        const Graph::vertex_descriptor v0Rc = v0 ^ 1;
        const Graph::vertex_descriptor v1Rc = v1 ^ 1;
        auto[eRc, edgeExists] = boost::edge(v0Rc, v1Rc, graph);
        SHASTA2_ASSERT(edgeExists);
        SHASTA2_ASSERT(graph[eRc].frequency == graph[e].frequency);

        edgesFound.insert(e);
        edgesFound.insert(eRc);

        edgePairs.emplace_back(EdgePair({e, eRc, graph[e].frequency}));
    }
    sort(edgePairs.begin(), edgePairs.end());



    // Do strand separation by adding edges in order of decreasing frequency.
    DisjointSets disjointSets(segments.size());
    uint64_t crossStrandEdgeCount = 0;
    for(const EdgePair& edgePair: edgePairs) {
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
    if(debug) {
        cout << "Found " << crossStrandEdgeCount << " cross-strand edges out of " <<
            num_edges(graph) << " total." << endl;
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
        ofstream dot("splitSelfComplementaryTangle.dot");
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


    if(debug) {
        writeHtmlEnd(html);
        cout << "AssemblyGraph::splitSelfComplementaryTangle ends for tangle " << tangleId << endl;
    }
}

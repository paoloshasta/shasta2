// Shasta2.
#include "ReadFollowing5.hpp"
#include "deduplicate.hpp"
#include "html.hpp"
#include "Options.hpp"
#include "SegmentStepSupport.hpp"
#include "Tangle.hpp"
using namespace shasta2;
using namespace ReadFollowing5;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "fstream.hpp"



Graph::Graph(
    const AssemblyGraph& assemblyGraph,
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html) :
    assemblyGraph(assemblyGraph),
    tangle(tangle),
    tangleId(tangleId),
    isSelfComplementaryTangle(tangle.isSelfComplementary())
{
    Graph& graph = *this;

    if(html) {
        html << "<h2>Read following</h2>";
    }

    createVertices();
    createEdges();

    if(html) {
        html << "The read following graph for this tangle has " <<
            num_vertices(graph) << " vertices and " <<
            num_edges(graph) << " edges." << endl;
        const string fileName = "ReadFollowingGraph-Tangle-" + to_string(tangleId) + ".dot";
        writeGraphviz(fileName);
        html << "<br>See " << fileName;
    }

    findShortestPaths();
    if(html) {
        writeShortestPaths(html);
    }

}



void Graph::createVertices()
{
    Graph& graph = *this;

    for(const Segment segment: tangle.tangleEdges) {
        add_vertex(Vertex(assemblyGraph, segment, false, false), graph);
    }
    for(const Segment segment: tangle.entrances) {
        entranceVertices.push_back(add_vertex(Vertex(assemblyGraph, segment, true, false), graph));
    }
    for(const Segment segment: tangle.exits) {
        exitVertices.push_back(add_vertex(Vertex(assemblyGraph, segment, false, true), graph));
    }



    // If this is a self-complementary tangle, fill in the vRc fields in the vertices.
    if(tangle.isSelfComplementary()) {

        // Fill in the vRc fields of the entrances and the exits.
        // Each exit is the reverse complements of an entrance.
        std::map<Segment, vertex_descriptor> entranceMap;
        for(const vertex_descriptor entrance: entranceVertices) {
            entranceMap.insert({graph[entrance].segment, entrance});
        }
        for(const vertex_descriptor exit: exitVertices) {
            const Segment exitSegment = graph[exit].segment;
            const Segment entranceSegment = assemblyGraph[exitSegment].eRc;
            const vertex_descriptor entrance = entranceMap.at(entranceSegment);
            graph[exit].vRc = entrance;
            graph[entrance].vRc = exit;
        }

        // Fill in the vRc fields of the internal vertices.
        std::map<Segment, vertex_descriptor> internalMap;
        BGL_FORALL_VERTICES(v, graph, Graph) {
            const Vertex& vertex = graph[v];
            if(not (vertex.isEntrance or vertex.isExit)) {
                internalMap.insert({vertex.segment, v});
            }
        }
        BGL_FORALL_VERTICES(v, graph, Graph) {
            const Vertex& vertex = graph[v];
            if(not (vertex.isEntrance or vertex.isExit)) {
                const Segment segment = vertex.segment;
                const Segment segmentRc = assemblyGraph[segment].eRc;
                const vertex_descriptor vRc = internalMap.at(segmentRc);
                graph[v].vRc = vRc;
            }
        }
        BGL_FORALL_VERTICES(v, graph, Graph) {
            const Vertex& vertex = graph[v];
            if(not (vertex.isEntrance or vertex.isExit)) {
                const vertex_descriptor vRc = vertex.vRc;
                SHASTA2_ASSERT(graph[vRc].vRc == v);
            }
        }
    }
}



Vertex::Vertex(
    const AssemblyGraph& assemblyGraph,
    Segment segment,
    bool isEntrance,
    bool isExit) :
    segment(segment),
    isEntrance(isEntrance),
    isExit(isExit)
{
    const uint32_t representativeRegionStepCount = uint32_t(assemblyGraph.options.representativeRegionStepCount);
    vector<SegmentStepSupport> support;

    // Fill in initial support.
    // Don't include ReadIds that appear on both strands.
    if(not isEntrance) {
        SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, support);
        initialSupport.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            initialSupport.push_back(orientedReadId);
        }
    }

    // Fill in final support.
    // Don't include ReadIds that appear on both strands.
    if(not isExit) {
        SegmentStepSupport::getFinalLast(assemblyGraph, segment, representativeRegionStepCount, support);
        finalSupport.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            finalSupport.push_back(orientedReadId);
        }
    }
}



void Graph::createEdges()
{
    Graph& graph = *this;
    const bool debug = false;

    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);
    const double a = 3.;                // dB
    const double b = 15.;               // dB
    const double logPThreshold = 10.;   // dB

    // No html output from analyzeSegmentPair.
    ostream html(0);



    // For each OrientedReadId, gather the vertices that have that OrientedReadId
    // in their initialSupport/finalSupport.
    std::map<OrientedReadId, vector<vertex_descriptor> > initialSupportMap;
    std::map<OrientedReadId, vector<vertex_descriptor> > finalSupportMap;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        for(const OrientedReadId orientedReadId: vertex.initialSupport) {
            initialSupportMap[orientedReadId].push_back(v);
        }
        for(const OrientedReadId orientedReadId: vertex.finalSupport) {
            finalSupportMap[orientedReadId].push_back(v);
        }
    }

    // Gather vertex pairs (v0, v1) that have common OrientedReadIds
    // between the final support of v0 and the intial support of v1.
    vector< pair<vertex_descriptor, vertex_descriptor> > vertexPairs;
    for(const auto&[orientedReadId, finalSupportVertices]: finalSupportMap) {
        const auto& initialSupportVertices = initialSupportMap[orientedReadId];
        for(const vertex_descriptor v0: finalSupportVertices) {
            for(const vertex_descriptor v1: initialSupportVertices) {
                vertexPairs.push_back(make_pair(v0, v1));
            }
        }
    }

    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(vertexPairs, count, minCommonCount);
    SHASTA2_ASSERT(vertexPairs.size() == count.size());



    // Each vertex pair can generate an edge, but we have to do it in a strand-symmetric way.
    vector<OrientedReadId> commonOrientedReadIds;
    for(const auto&[v0, v1]: vertexPairs) {
        if(v0 == v1) {
            continue;
        }
        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];
        const Segment segment0 = vertex0.segment;
        const Segment segment1 = vertex1.segment;

        // If this is a self-complementary tangle, generate edges in pairs.
        // See if this a vertex pair that we want to use (to ensure strand symmetry).
        if(isSelfComplementaryTangle) {
            const uint64_t id0 = assemblyGraph[segment0].id;
            const uint64_t id1 = assemblyGraph[segment1].id;
            const Segment segment0Rc = assemblyGraph[segment0].eRc;
            const Segment segment1Rc = assemblyGraph[segment1].eRc;
            SHASTA2_ASSERT(segment0Rc != assemblyGraphNullEdge);
            SHASTA2_ASSERT(segment1Rc != assemblyGraphNullEdge);
            SHASTA2_ASSERT(segment0Rc != segment0);
            SHASTA2_ASSERT(segment0Rc != segment1);
            SHASTA2_ASSERT(segment1Rc != segment0);
            SHASTA2_ASSERT(segment1Rc != segment1);
            const uint64_t id0Rc = assemblyGraph[segment0Rc].id;
            const uint64_t id1Rc = assemblyGraph[segment1Rc].id;
            if(min(id0, id1) >= min(id0Rc, id1Rc)) {
                continue;
            }
        }

        const SegmentPairInformation segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
            html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

        if(debug) {
            cout << "Edge candidate " << assemblyGraph.id(segment0) << " " <<
                assemblyGraph.id(segment1) << ", offset " << segmentPairInformation.segmentOffset <<
                ", common count " << segmentPairInformation.commonCount <<
                ", missing " << segmentPairInformation.missing() << endl;
        }

        // If the offset is negative, discard it.
        if(segmentPairInformation.segmentOffset < 0) {
            if(debug) {
                cout << "Discarded due to negative offset." << endl;
            }
            continue;
        }

        // If it does not satisfy our requirements, discard it.
        if(segmentPairInformation.commonCount < minCommonCount) {
            if(debug) {
                cout << "Discarded due low common count." << endl;
            }
            continue;
        }
        const double logP = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing());
        if(logP < logPThreshold) {
            if(debug) {
                cout << "Discarded due low logP." << endl;
            }
            continue;
        }
        if(not assemblyGraph.canConnect(segment0, segment1)) {
            if(debug) {
                cout << "Discarded because cannot connect." << endl;
            }
            continue;
        }


        if(debug) {
            cout << "Keeping " << assemblyGraph.id(segment0) << " " << assemblyGraph.id(segment1) << " " << logP << endl;
        }
        const Edge edge(segmentPairInformation.commonCount, logP);
        add_edge(v0, v1, edge, graph);
        if(isSelfComplementaryTangle) {
            const vertex_descriptor v0Rc = graph[v0].vRc;
            const vertex_descriptor v1Rc = graph[v1].vRc;
            add_edge(v1Rc, v0Rc, edge, graph);
            if(debug) {
                cout << "Keeping " << assemblyGraph.id(graph[v1Rc].segment) << " " <<
                    assemblyGraph.id(graph[v0Rc].segment) << " " << logP << endl;
            }
        }
    }
}



void Graph::writeGraphviz(const string& fileName) const
{
    ofstream dot(fileName);
    writeGraphviz(dot);
}



void Graph::writeGraphviz(ostream& dot) const
{
    const Graph& graph = *this;
    dot << "digraph ReadFollowingGraph {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        dot << assemblyGraph.id(graph[v].segment);
        if(graph[v].isEntrance) {
            dot << "[style=filled fillcolor=pink]";
        }
        if(graph[v].isExit) {
            dot << "[style=filled fillcolor=green]";
        }
        dot << ";\n";
    }

    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        dot << assemblyGraph.id(graph[v0].segment) << "->" <<
            assemblyGraph.id(graph[v1].segment) <<
            "[ label=\"" << int64_t(std::round(graph[e].logP)) << "\"]"
            "\n";
    }

    dot << "}\n";
}


void Graph::findShortestPaths()
{
    const Graph& graph = *this;
    using boost::make_assoc_property_map;

    shortestPaths.resize(entranceVertices.size(), vector<ShortestPath>(exitVertices.size()));

    // Create a vertex index map, needed below.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Loop over entrances.
    for(uint64_t iEntrance=0; iEntrance<entranceVertices.size(); iEntrance++) {
        const vertex_descriptor vEntrance = entranceVertices[iEntrance];

        // Create a shortest path tree rooted at this entrance.
        std::map<vertex_descriptor, double> distanceMap;
        std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
        dijkstra_shortest_paths(graph, vEntrance,
           weight_map(boost::get(&Edge::weight, graph)).
           vertex_index_map(make_assoc_property_map(vertexIndexMap)).
           distance_map(make_assoc_property_map(distanceMap)).
           predecessor_map(make_assoc_property_map(predecessorMap)));

        // Write out the distances of the exits from this entrance.
        for(uint64_t iExit=0; iExit<exitVertices.size(); iExit++) {
            const vertex_descriptor vExit = exitVertices[iExit];

            ShortestPath& shortestPath = shortestPaths[iEntrance][iExit];
            shortestPath.entrance = graph[vEntrance].segment;
            shortestPath.exit = graph[vExit].segment;

            // If this exit is not reachable from this entrance,
            // leave the distance set to invalid<double> and the path empty.
            const double distance = distanceMap.at(vExit);
            if(distance == std::numeric_limits<double>::max()) {
                continue;
            }
            shortestPath.distance = 10. * log10(distance);

            // Walk back the shortest path tree to find the path Segments.
            vertex_descriptor v = vExit;
            while(true) {
                const vertex_descriptor vPrevious = predecessorMap[v];
                if(vPrevious == vEntrance) {
                    break;
                }
                SHASTA2_ASSERT(vPrevious != v);
                v = vPrevious;
                shortestPath.segments.push_back(graph[v].segment);
            }
            std::ranges::reverse(shortestPath.segments);

        }
    }
}



void Graph::writeShortestPaths(ostream& html)
{
    // Write the distance matrix.
    html <<
        "<h3>Distance matrix</h3>"
        "Distance on the read following graph "
        "between between each pair of entrances and exits. "
        "A small value indicates that the entrance and exit are likely to "
        " follow each other in genomic sequence. "
        "Blank values indicate that the exit is unreachable from that entrance."
        "<br><br><table>";

    // Write a header line listing the exits.
    html << "<tr><th>";
    for(const Segment& segment: tangle.exits) {
        html << "<th>" << assemblyGraph.id(segment);
    }

    // Write one row for each entrance.
    for(uint64_t iEntrance=0; iEntrance<entranceVertices.size(); iEntrance++) {
        html << "<tr><th>" << assemblyGraph.id(tangle.entrances[iEntrance]);
        for(uint64_t iExit=0; iExit<exitVertices.size(); iExit++) {
            html << "<td class=centered>";
            const double distance = shortestPaths[iEntrance][iExit].distance;
            if(distance != invalid<double>) {
                html << int64_t(std::round(shortestPaths[iEntrance][iExit].distance));
            }
        }
    }
    html << "</table>";



    // Write the details of the shortest paths.
    html <<
        "<h3>Shortest paths</h3>"
        "<table>"
        "<tr><th>Entrance<th>Exit<th>Distance<th>Path";
    for(uint64_t iEntrance=0; iEntrance<entranceVertices.size(); iEntrance++) {
        const Segment entrance = tangle.entrances[iEntrance];
        for(uint64_t iExit=0; iExit<exitVertices.size(); iExit++) {
            const Segment exit = tangle.exits[iExit];
            const ShortestPath& shortestPath = shortestPaths[iEntrance][iExit];
            SHASTA2_ASSERT(shortestPath.entrance == entrance);
            SHASTA2_ASSERT(shortestPath.exit == exit);
            if(shortestPath.isUnreachable()) {
                continue;
            }

            // Entrance and exit.
            html <<
                "<tr>"
                "<td class=centered>" << assemblyGraph.id(entrance) <<
                "<td class=centered>" << assemblyGraph.id(exit);

            // Distance.
            html << "<td class=centered>";
            html << int64_t(std::round(shortestPath.distance));

            // Path.
            html << "<td class=left>";
            for(uint64_t i=0; i<shortestPath.segments.size(); i++) {
                if(i != 0) {
                    html << ",";
                }
                html << assemblyGraph.id(shortestPath.segments[i]);
            }
         }
    }

    html << "</table>";
}



// Return the number of exits reachable from an entrance.
uint64_t Graph::countReachableExits(uint64_t iEntrance) const
{
    uint64_t n = 0;
    for(uint64_t iExit=0; iExit<tangle.exits.size(); iExit++) {
        if(shortestPaths[iEntrance][iExit].isReachable()) {
            ++n;
        }
    }
    return n;
}



// Return the number of entrances from which an exit is reachable.
uint64_t Graph::countReachableEntrances(uint64_t iExit) const
{
    uint64_t n = 0;
    for(uint64_t iEntrance=0; iEntrance<tangle.entrances.size(); iEntrance++) {
        if(shortestPaths[iEntrance][iExit].isReachable()) {
            ++n;
        }
    }
    return n;
}


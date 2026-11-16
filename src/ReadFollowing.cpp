// Shasta2.
#include "ReadFollowing.hpp"
#include "CondensedGraph.hpp"
#include "deduplicate.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "html.hpp"
#include "Options.hpp"
#include "SegmentStepSupport.hpp"
#include "Tangle.hpp"
#include "transitiveReduction.hpp"
using namespace shasta2;
using namespace ReadFollowing;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "fstream.hpp"



Graph::Graph(
    AssemblyGraph& assemblyGraph,
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html,
    const string& debugOutputBaseName) :
    assemblyGraph(assemblyGraph),
    tangleId(tangleId),
    tangle(tangle),
    html(html),
    debugOutputBaseName(debugOutputBaseName),
    isSelfComplementaryTangle(tangle.isSelfComplementary())
{

    if(html) {
        html << "<h2>Simple read following</h2>";
    }

    createVertices();
    createReadSimilarityEdges();
    writeGraphToHtml("Initial read following graph", "A");

    findShortestPaths();
    writeShortestPaths();

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


    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexMap.insert({graph[v].segment, v});
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



void Graph::createReadSimilarityEdges()
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


    // Fot better display, write vertices in topological order
    // of the CondensedGraph.
    std::map<vertex_descriptor, CondensedGraph::vertex_descriptor> condensedVertexMap;
    const CondensedGraph condensedGraph = createCondensedGraph<Graph, CondensedGraph>(graph, condensedVertexMap);
    vector<CondensedGraph::vertex_descriptor> topologicalOrder;
    boost::topological_sort(condensedGraph, back_inserter(topologicalOrder));
    std::ranges::reverse(topologicalOrder);



    dot << "digraph ReadFollowingGraph {\n";

    BGL_FORALL_VERTICES(cv, condensedGraph, CondensedGraph) {
        for(const vertex_descriptor v: condensedGraph[cv].vertices) {
            dot << assemblyGraph.id(graph[v].segment);
            if(graph[v].isEntrance) {
                dot << "[style=filled fillcolor=pink]";
            }
            if(graph[v].isExit) {
                dot << "[style=filled fillcolor=green]";
            }
            dot << ";\n";
        }
    }



    BGL_FORALL_EDGES(e, graph, Graph) {
        const Edge& edge = graph[e];
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // For clarity, write the edge as one or two arrows.

        if(edge.haveSimilarReadCompositions) {
            dot <<
                assemblyGraph.id(graph[v0].segment) << "->" <<
                assemblyGraph.id(graph[v1].segment) <<
                "["
                "label=\"" << int64_t(std::round(graph[e].logP)) << "\""
                " color=orange"
                "]\n";
        }

        if(edge.areAdjacent) {
            dot <<
                assemblyGraph.id(graph[v0].segment) << "->" <<
                assemblyGraph.id(graph[v1].segment) << "\n";
        }

    }

    dot << "}\n";
}



void Graph::writeGraphToHtml(
    const string& title,
    const string& fileNameSuffix) const
{
    const Graph& graph = *this;

    if(not html) {
        return;
    }

    html <<
        "<h3>" << title << "</h3>"
        "The read following graph has " <<
        num_vertices(graph) << " vertices and " <<
        num_edges(graph) << " edges." << endl;

    const string fileName =
        debugOutputBaseName +
        "-ReadFollowingGraph-" +
        fileNameSuffix +
        "-Tangle-" +
        to_string(tangleId) +
        ".dot";

    writeGraphviz(fileName);
    html << "<br>See " << fileName << "<br>";

    try {
        const double timeout = 30.;
        graphvizToHtml(fileName, "dot", timeout, "-Nshape=rectangle", html, true);
    } catch(std::exception& e) {
        html << "<br>The read following graph could not be displayed: " << e.what();
    }
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



void Graph::writeShortestPaths()
{
    if(not html) {
        return;

    }

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




bool Graph::updateAssemblyGraph1()
{
    // For now, require a one-on-one matching of the entrances
    // and exits.
    for(uint64_t iEntrance=0; iEntrance<tangle.entrances.size(); iEntrance++) {
        if(countReachableExits(iEntrance) != 1) {
            if(html) {
                html << "<br>Simple read following failed because of ambiguous or missing reachability.";
            }
            return false;
        }
    }
    for(uint64_t iExit=0; iExit<tangle.exits.size(); iExit++) {
        if(countReachableEntrances(iExit) != 1) {
            if(html) {
                html << "<br>Simple read following failed because of ambiguous or missing reachability.";
            }
            return false;
        }
    }

    // If getting here, read following can proceed.
    if(html) {
        html << "<br>Simple read following can proceed for this tangle.";
    }



    // We need to disconnect (at the beginning or end) the entrances and exits.
    // To do this we use disconnectAtBeginning and disconnectAtEnd, which change
    // edge_descriptors but leave the ids unchanged.
    // Some Segments can at the same time be an entrance and/or exit
    // of this tangle or its reverse complement. To avoid working with invalidated
    // edge_descriptors, we work with Segment ids instead.
    std::map<uint64_t, Segment> segmentMap;
    tangle.createSegmentMap(segmentMap);
    vector<uint64_t> entranceIds;
    for(const Segment e: tangle.entrances) {
        entranceIds.push_back(assemblyGraph.id(e));
    }
    vector<uint64_t> exitIds;
    for(const Segment e: tangle.exits) {
        exitIds.push_back(assemblyGraph.id(e));
    }

    // Create disconnected versions of the entrances and exits,
    // while keeping the segmentMap up to date.
    for(const uint64_t entranceId: entranceIds) {
        const Segment eNew = assemblyGraph.disconnectAtEnd(segmentMap.at(entranceId));
        segmentMap[assemblyGraph.id(eNew)] = eNew;
        const Segment eNewRc = assemblyGraph[eNew].eRc;
        segmentMap[assemblyGraph.id(eNewRc)] = eNewRc;
    }
    if(not isSelfComplementaryTangle) {
        for(const uint64_t exitId: exitIds) {
            const Segment eNew = assemblyGraph.disconnectAtBeginning(segmentMap.at(exitId));
            segmentMap[assemblyGraph.id(eNew)] = eNew;
            const Segment eNewRc = assemblyGraph[eNew].eRc;
            segmentMap[assemblyGraph.id(eNewRc)] = eNewRc;
        }
    }
    tangle.checkSegmentMap(segmentMap);



    // Now we can connect entrances with exits using the shortest paths stored in the Graph.
    std::set<uint64_t> connectedEntrances;
    for(uint64_t iEntrance=0; iEntrance<tangle.entrances.size(); iEntrance++) {
        const uint64_t entranceId = entranceIds[iEntrance];
        if(isSelfComplementaryTangle and (connectedEntrances.contains(entranceId))) {
            continue;
        }
        const Segment entrance = segmentMap.at(entranceId);
        for(uint64_t iExit=0; iExit<tangle.exits.size(); iExit++) {
            const uint64_t exitId = exitIds[iExit];
            const Segment exit = segmentMap.at(exitId);
            const auto& shortestPath = shortestPaths[iEntrance][iExit];
            if(shortestPath.isUnreachable()) {
                continue;
            }

            // Get the Segments internal to the path.
            // We will make disconnected copies of these segments.
            const vector<Segment>& oldAssemblyPath = shortestPath.segments;



            // Create a complete path consisting of the entrance, the
            // newly created copies of the internal segments, and the exit.
            // We will still have to fill in the portions in-between these segments.
            vector<Segment> assemblyPath;
            assemblyPath.push_back(entrance);
            for(const Segment oldSegment: oldAssemblyPath) {
                const Segment newSegment = assemblyGraph.createDisconnectedCopy(oldSegment);
                assemblyPath.push_back(newSegment);
                const vertex_descriptor v0 = source(newSegment, assemblyGraph);
                const vertex_descriptor v1 = target(newSegment, assemblyGraph);
                assemblyGraph.createReverseComplementVertex(v0);
                assemblyGraph.createReverseComplementVertex(v1);
                assemblyGraph.createReverseComplementEdge(newSegment);
            }
            assemblyPath.push_back(exit);

            // Now we have to fill in the portions in-between these segments.
            for(uint64_t i1=1; i1<assemblyPath.size(); i1++) {
                const uint64_t i0 = i1 - 1;
                const Segment segment0 = assemblyPath[i0];
                const Segment segment1 = assemblyPath[i1];
                const Segment newSegment = assemblyGraph.connect(segment0, segment1);
                assemblyGraph.createReverseComplementEdge(newSegment);
            }



            if(isSelfComplementaryTangle) {
                connectedEntrances.insert(entranceId);
                connectedEntrances.insert(assemblyGraph.id(assemblyGraph[exit].eRc));
            }
        }
    }



    // Remove all the Segments internal to this tangle
    // and their reverse complements.
    for(const Segment e: tangle.tangleEdges) {
        if(not isSelfComplementaryTangle) {
            boost::remove_edge(assemblyGraph[e].eRc, assemblyGraph);
        }
        boost::remove_edge(e, assemblyGraph);
    }

    // If getting here, read following was successful on this tangle.
    return true;

}



bool Graph::updateAssemblyGraph2()
{
    Graph& graph = *this;

    if(html) {
        html << "<h2>Complex read following</h2>" << endl;
    }

    // For now skip the self-complementary case.
    if(isSelfComplementaryTangle) {
        if(html) {
            html << "<br>Skipping complex read following of a self-complementaty tangle." << endl;
        }
        return false;
    }

    partialTransitiveReduction();
    writeGraphToHtml("Read following graph after transitive reduction", "B");

    // We only want unambiguous vertex chains, so we disconnect all branch vertices;
    disconnectBranchVertices();
    writeGraphToHtml("Read following graph after disconnecting branch vertices", "C");

    createAdjacencyEdges();
    writeGraphToHtml("Read following graph after adding adjacency edges", "D");

    const uint64_t adjacencyEdgesRemoved = cleanupEdges();
    if(adjacencyEdgesRemoved == 0) {
        if(html) {
            html << "<br>Complex read following failed because it would not simplify the AssemblyGraph.";
        }
        return false;
    }
    writeGraphToHtml("Read following graph after edge cleanup", "E");

    if(not removeInaccessibleVertices()) {
        if(html) {
            html << "<br>Complex read following failed due to accessibility." << endl;
        }
        return false;
    }
    writeGraphToHtml("Read following graph after removing inaccessible vertices", "F");



    // If getting here, read following was successful.



    // For each vertex still surviving that,
    // make a disconnected copy of the corresponding segment.
    html << flush;
    std::map<Segment, Segment> newSegmentMap;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        const Segment segment = vertex.segment;

        Segment newSegment;
        if(vertex.isEntrance) {
            if(vertex.isExit) {
                SHASTA2_ASSERT(0);  // FOR NOW
                newSegment = assemblyGraph.createDisconnectedCopy(segment); // Does not remove the original.
                if(not isSelfComplementaryTangle) {
                    // Entrance and exit.
                    const AssemblyGraph::vertex_descriptor v0New = source(newSegment, assemblyGraph);
                    const AssemblyGraph::vertex_descriptor v1New = target(newSegment, assemblyGraph);
                    assemblyGraph.createReverseComplementVertex(v0New);
                    assemblyGraph.createReverseComplementVertex(v1New);
                    assemblyGraph.createReverseComplementEdge(newSegment);
                }
            } else {
                // Entrance, but not an exit.
                // This also creates the reverse complement and removes the two original ones.
                newSegment = assemblyGraph.disconnectAtEnd(segment);
#if 0
                if(not isSelfComplementaryTangle) {
                    const AssemblyGraph::vertex_descriptor v1New = target(newSegment, assemblyGraph);
                    assemblyGraph.createReverseComplementVertex(v1New);
                    const Segment newSegmentRc = assemblyGraph.createReverseComplementEdge(newSegment);
                    cout << "BBB created " << assemblyGraph.id(newSegment) << " " << assemblyGraph.id(newSegmentRc) << endl;
                }
#endif
            }
        } else if(vertex.isExit) {
            // Exit but not entrance.
            // This also creates the reverse complement and removes the two original ones.
            SHASTA2_ASSERT(not vertex.isEntrance);
            newSegment = assemblyGraph.disconnectAtBeginning(segment);
        } else {
            // Not entrance and not exit.
            newSegment = assemblyGraph.createDisconnectedCopy(segment); // Does not remove the original.
            if(not isSelfComplementaryTangle) {
                const AssemblyGraph::vertex_descriptor v0New = source(newSegment, assemblyGraph);
                const AssemblyGraph::vertex_descriptor v1New = target(newSegment, assemblyGraph);
                assemblyGraph.createReverseComplementVertex(v0New);
                assemblyGraph.createReverseComplementVertex(v1New);
                assemblyGraph.createReverseComplementEdge(newSegment);
            }
        }

        newSegmentMap.insert({segment, newSegment});
    }



    // For each edge, connect the new segments.
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;
        const Segment newSegment0 = newSegmentMap.at(segment0);
        const Segment newSegment1 = newSegmentMap.at(segment1);
        const Segment newSegment = assemblyGraph.connect(newSegment0, newSegment1);
        if(not isSelfComplementaryTangle) {
            assemblyGraph.createReverseComplementEdge(newSegment);
        }
    }


    // Now we can remove all the internal Segments of the Tangle.
    // This leaves some isolated vertices and some zero length segments
    // that will be taken care of later.
    for(const Segment segment: tangle.tangleEdges) {
        const Segment segmentRc = assemblyGraph[segment].eRc;
        boost::remove_edge(segment, assemblyGraph);
        if(not isSelfComplementaryTangle) {
            boost::remove_edge(segmentRc, assemblyGraph);
        }
    }

    return true;
}



// If a vertex has any out-edges with haveSimilarReadCompositions set
// (there can be at most one), remove all other out-edges.
// If a vertex has any in-edges with haveSimilarReadCompositions set
// (there can be at most one), remove all other in-edges.
uint64_t Graph::cleanupEdges()
{
    Graph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {

        bool removeOut = false;
        BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
            if(graph[e].haveSimilarReadCompositions) {
                removeOut = true;
                break;
            }
        }
        if(removeOut) {
            BGL_FORALL_OUTEDGES(v, e, graph, Graph) {
                if(not graph[e].haveSimilarReadCompositions) {
                    edgesToBeRemoved.push_back(e);
                }
            }
        }

        bool removeIn = false;
        BGL_FORALL_INEDGES(v, e, graph, Graph) {
            if(graph[e].haveSimilarReadCompositions) {
                removeIn = true;
                break;
            }
        }
        if(removeIn) {
            BGL_FORALL_INEDGES(v, e, graph, Graph) {
                if(not graph[e].haveSimilarReadCompositions) {
                    edgesToBeRemoved.push_back(e);
                }
            }
        }
    }
    deduplicate(edgesToBeRemoved);

    uint64_t adjacencyEdgesRemoved = 0;
    for(const edge_descriptor e: edgesToBeRemoved) {
        if(html) {
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            html << "<br>Removing edge " << assemblyGraph.id(graph[v0].segment) <<
                " " << assemblyGraph.id(graph[v1].segment);
        }
        if(graph[e].areAdjacent) {
            ++adjacencyEdgesRemoved;
        }
        boost::remove_edge(e, graph);
    }

    return adjacencyEdgesRemoved;
}



void Graph::createAdjacencyEdges()
{
    Graph& graph = *this;

    for(const AssemblyGraph::vertex_descriptor v: tangle.tangleVertices) {
        BGL_FORALL_INEDGES(v, segment0, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor u0 = vertexMap.at(segment0);
            BGL_FORALL_OUTEDGES(v, segment1, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor u1 = vertexMap.at(segment1);

                const auto[e, edgeExists] = boost::edge(u0, u1, *this);
                if(edgeExists) {
                    graph[e].areAdjacent = true;
                } else {
                    boost::add_edge(u0, u1, *this);
                }
            }
        }
    }
}



void Graph::partialTransitiveReduction()
{
    Graph& graph = *this;

    // Create the CondensedGraph.
    std::map<vertex_descriptor, CondensedGraph::vertex_descriptor> condensedVertexMap;
    CondensedGraph condensedGraph = createCondensedGraph<Graph, CondensedGraph>(graph, condensedVertexMap);

    // Transitive reduction of the CondensedGraph.
    // The CondensedGraph is guaranteed to be acyclic, so this will not fail.
    transitiveReductionAny(condensedGraph);



    // Now remove transitive reduction edges from the Graph.
    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        // If v0 or v1 are in non-trivial strong components,
        // leave this edge alone.
        const auto cv0 = condensedVertexMap.at(v0);
        if(condensedGraph[cv0].vertices.size() > 1) {
            continue;
        }
        const auto cv1 = condensedVertexMap.at(v1);
        if(condensedGraph[cv1].vertices.size() > 1) {
            continue;
        }

        const auto[_, edgeExists] = boost::edge(cv0, cv1, condensedGraph);
        if(not edgeExists) {
            edgesToBeRemoved.push_back(e);
        }
    }
    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }

}



void Graph::disconnectBranchVertices()
{
    Graph& graph = *this;

    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((out_degree(v, graph) > 1) or (in_degree(v, graph) > 1 )) {
            boost::clear_vertex(v, graph);
        }
    }
}



// Remove inaccessible vertices.
// This returns false if any entrances or exits
// are inaccessible.
bool Graph::removeInaccessibleVertices()
{
    Graph& graph = *this;

    // Find vertices reachable from the entrances.
    vector<vertex_descriptor> verticesReachableFromEntrances;
    for(const vertex_descriptor entranceVertex: entranceVertices) {
        std::set<vertex_descriptor> verticesReachableFromThisEntrance;
        findReachableVertices(graph, entranceVertex, 0, verticesReachableFromThisEntrance);
        std::ranges::copy(verticesReachableFromThisEntrance, back_inserter(verticesReachableFromEntrances));
    }
    deduplicate(verticesReachableFromEntrances);

    // Find vertices reachable from the exits.
    vector<vertex_descriptor> verticesReachableFromExits;
    for(const vertex_descriptor exitVertex: exitVertices) {
        std::set<vertex_descriptor> verticesReachableFromThisExit;
        findReachableVertices(graph, exitVertex, 1, verticesReachableFromThisExit);
        std::ranges::copy(verticesReachableFromThisExit, back_inserter(verticesReachableFromExits));
    }
    deduplicate(verticesReachableFromExits);

    // Check that the exits are reachable from the entrances.
    for(const vertex_descriptor exitVertex: exitVertices) {
        if(not std::ranges::binary_search(verticesReachableFromEntrances, exitVertex)) {
            return false;
        }
    }
    for(const vertex_descriptor entranceVertex: entranceVertices) {
        if(not std::ranges::binary_search(verticesReachableFromExits, entranceVertex)) {
            return false;
        }
    }

    vector<vertex_descriptor> reachableVertices;
    std::ranges::set_intersection(verticesReachableFromEntrances, verticesReachableFromExits,
        back_inserter(reachableVertices));

    // Now we can remove the unreachable vertices.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not std::ranges::binary_search(reachableVertices, v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }


    return true;
}

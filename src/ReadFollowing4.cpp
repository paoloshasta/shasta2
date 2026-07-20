// Shasta.
#include "ReadFollowing4.hpp"
#include "Journeys.hpp"
#include "memoryInformation.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "SegmentStepSupport.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta2;
using namespace ReadFollowing4;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"

// Standard library.
#include "fstream.hpp"
#include <queue>
#include <ranges>

// Explicit instantiation.
#include "MultithreadedObject.tpp"
template class MultithreadedObject<ReadFollower>;



ReadFollower::ReadFollower(const AssemblyGraph& assemblyGraph) :
    MultithreadedObject<ReadFollower>(*this),
    assemblyGraph(assemblyGraph)
{
    const bool debug = false;
    writeMemoryStatistics("ReadFollower::ReadFollower begin");

    fillSupportMaps();
    findSegmentPairs();
    createVertices();
    createEdges();

    if(debug) {
        cout << "The initial read following search graph has " << num_vertices(searchGraph) <<
            " vertices and " << num_edges(searchGraph) << " edges." << endl;
    }
    searchGraph.writeGraphviz(assemblyGraph, "Initial");
    searchGraph.check(assemblyGraph);

    // Prune.
    searchGraph.prune();
    searchGraph.writeGraphviz(assemblyGraph, "Pruned");
    searchGraph.check(assemblyGraph);

    if(true) {
        cout << "After pruning, the read following search graph has " << num_vertices(searchGraph) <<
            " vertices and " << num_edges(searchGraph) << " edges." << endl;
        cout << "The read following graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges." << endl;
    }
    graph.writeGraphviz(assemblyGraph, "A");
    graph.check(assemblyGraph);

    // Before we can compute shortest paths we have to create the vertex index map.
    searchGraph.createVertexIndexMap();

    // Use the SearchGraphs to find shortest paths between long segments
    // and store them in the ConnectGraph.
    findShortestPathsMultithreaded(assemblyGraph.options.threadCount);
    graph.check(assemblyGraph);

    if(debug) {
        cout << "After finding shortest paths, the read following graph has " << num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges." << endl;
    }
    graph.writeGraphviz(assemblyGraph, "B");

    graph.removeWeakEdges();
    graph.writeGraphviz(assemblyGraph, "C");
    graph.check(assemblyGraph);

    graph.transitiveReduction();
    graph.writeGraphviz(assemblyGraph, "D");
    graph.check(assemblyGraph);

    writeMemoryStatistics("ReadFollower::ReadFollower end");
}



void ReadFollower::fillSupportMaps()
{
    const uint32_t representativeRegionStepCount = uint32_t(assemblyGraph.options.representativeRegionStepCount);

    vector<SegmentStepSupport> support;
    vector<OrientedReadId> orientedReadIds;

    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {

        // Initial support. Remove duplicate ReadIds.
        SegmentStepSupport::getInitialFirst(assemblyGraph, segment, representativeRegionStepCount, support);
        orientedReadIds.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            orientedReadIds.push_back(orientedReadId);
        }
        initialSupportMap.insert(make_pair(segment, orientedReadIds));



        // Final support. Remove duplicate ReadIds.
        SegmentStepSupport::getFinalLast(assemblyGraph, segment, representativeRegionStepCount, support);
        orientedReadIds.clear();
        for(uint64_t i=0; i<support.size(); i++) {
            const OrientedReadId orientedReadId = support[i].orientedReadId;
            const ReadId readId = orientedReadId.getReadId();
            if((i != 0) and (readId == support[i-1].orientedReadId.getReadId())) {
                continue;
            }
            if((i != support.size()-1) and (readId == support[i+1].orientedReadId.getReadId())) {
                continue;
            }
            orientedReadIds.push_back(orientedReadId);
        }
        finalSupportMap.insert(make_pair(segment, orientedReadIds));
    }


    // Check that the final support for a Segment is
    // the reverse complement of the initial support
    // for the reverse complement of the segment.
    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        SHASTA2_ASSERT(edge.eRc != assemblyGraphNullEdge);
        const Segment eRc = edge.eRc;
        const AssemblyGraphEdge& edgeRc = assemblyGraph[eRc];
        SHASTA2_ASSERT(edgeRc.eRc == e);
        const vector<OrientedReadId>& support = finalSupportMap.at(e);
        const vector<OrientedReadId>& supportRc = initialSupportMap.at(eRc);
        const uint64_t n = support.size();
        SHASTA2_ASSERT(supportRc.size() == n);
        /*
        cout << "AAA " << edge.id << " " << edgeRc.id << endl;
        for(uint64_t i=0; i<n; i++) {
            const OrientedReadId orientedReadId = support[i];
            const OrientedReadId orientedReadIdRc = supportRc[i];
            cout << orientedReadId << " " << orientedReadIdRc << endl;
        }
        */
        for(uint64_t i=0; i<n; i++) {
            const OrientedReadId orientedReadId = support[i];
            const OrientedReadId orientedReadIdRc = supportRc[i];
            SHASTA2_ASSERT(orientedReadId.getReadId() == orientedReadIdRc.getReadId());
            SHASTA2_ASSERT(orientedReadId.getStrand() == 1 - orientedReadIdRc.getStrand());
        }
    }

}



void ReadFollower::findSegmentPairs()
{
    const uint64_t orientedReadCount = assemblyGraph.journeys.size();
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;

    // For each OrientedReadId, gather the Segments that the OrientedReadId
    // appears in, in the initial support.
    vector< vector<Segment> > initialSupportSegments(orientedReadCount);
    for(const auto&[segment, orientedReadIds]: initialSupportMap) {
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            initialSupportSegments[orientedReadId.getValue()].push_back(segment);
        }

    }

    // For each OrientedReadId, gather the Segments that the OrientedReadId
    // appears in, in the final support.
    vector< vector<Segment> > finalSupportSegments(orientedReadCount);
    for(const auto&[segment, orientedReadIds]: finalSupportMap) {
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            finalSupportSegments[orientedReadId.getValue()].push_back(segment);
        }

    }



    // Store all segment pairs (segment0, segment1) such that the final support
    // of segment0 shares one or more OrientedReadIds with the initial support of segment1.
    for(uint64_t i=0; i<orientedReadCount; i++) {
        const vector<Segment>& initialSegments = initialSupportSegments[i];
        const vector<Segment>& finalSegments = finalSupportSegments[i];
        for(const Segment segment0: finalSegments) {
            for(const Segment segment1: initialSegments) {
                if(segment0 != segment1) {
                    segmentPairs.push_back(make_pair(segment0, segment1));
                }
            }
        }
    }



    // Only keep the ones that appear at least minCommonCount times.
    vector<uint64_t> count;
    deduplicateAndCountWithThreshold(segmentPairs, count, minCommonCount);
    SHASTA2_ASSERT(segmentPairs.size() == count.size());

    cout << "Found " << segmentPairs.size() << " segment pairs." << endl;
}



void ReadFollower::createVertices()
{
    const uint64_t lengthThreshold = assemblyGraph.options.readFollowingSegmentLengthThreshold;

    uint64_t totalSegmentCount = 0;
    uint longSegmentCount = 0;

    // Loop over all AssemblyGraph Segments.
    BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
        ++totalSegmentCount;

        // Get the length and check if it qualifies as long.
        const uint64_t length = assemblyGraph[segment].length();
        const bool isLong = (length >= lengthThreshold);

        // Add a SearchGraph vertex.
        searchGraph.createVertex(segment, length, isLong);

        // If isLong, also add a vertex to the ConnectGraph.
        if(isLong) {
            graph.createVertex(segment, length);
            ++longSegmentCount;
        }
    }

    cout << "Read following found " << totalSegmentCount <<
        " segment of which " << longSegmentCount <<
        " are at least " << lengthThreshold << " bases long." << endl;
}



void SearchGraph::createVertex(Segment segment, uint64_t length, bool isLong)
{
    SHASTA2_ASSERT(not vertexMap.contains(segment));
    SearchGraph& graph = *this;
    const vertex_descriptor v = add_vertex(SearchGraphVertex(segment, length, isLong), graph);
    vertexMap.insert(make_pair(segment, v));
}



SearchGraphVertex::SearchGraphVertex(
    Segment segment,
    uint64_t length,
    bool isLong) :
    segment(segment),
    length(length),
    isLong(isLong)
{
}



void ConnectGraph::createVertex(Segment segment, uint64_t length)
{
    SHASTA2_ASSERT(not vertexMap.contains(segment));
    ConnectGraph& graph = *this;
    const vertex_descriptor v = add_vertex(ConnectGraphVertex(segment, length), graph);
    vertexMap.insert(make_pair(segment, v));
}



ConnectGraphVertex::ConnectGraphVertex(
    Segment segment,
    uint64_t length) :
    segment(segment),
    length(length)
{
}



void ReadFollower::createEdges()
{
    uint64_t threadCount = assemblyGraph.options.threadCount;
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }

    setupLoadBalancing(segmentPairs.size(), 100);
    runThreads(&ReadFollower::createEdgesThreadFunction, threadCount);

}



void ReadFollower::createEdgesThreadFunction([[maybe_unused]] uint64_t threadId)
{
    const uint64_t minCommonCount = assemblyGraph.options.readFollowingMinCommonCount;
    const uint32_t representativeRegionStepCount =  uint32_t(assemblyGraph.options.representativeRegionStepCount);

    // No html output from analyzeSegmentPair.
    ostream html(0);



    class SegmentPair {
    public:
        Segment segment0;
        Segment segment1;
        SegmentPairInformation segmentPairInformation;
        double logP;
        double logPForward;
        double logPBackward;
        SegmentPair(
            Segment segment0,
            Segment segment1,
            const SegmentPairInformation& segmentPairInformation) :
                segment0(segment0),
                segment1(segment1),
                segmentPairInformation(segmentPairInformation)
        {
            logP         = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing());
            logPForward  = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing0);
            logPBackward = a * double(segmentPairInformation.commonCount) - b * double(segmentPairInformation.missing1);
        }
    };

    // SegmentPairs that will generate edes in the SearchGraph.
    vector<SegmentPair> edgesToBeAdded;

    // SegmentPairs that will generate edges in the ConnectGraph.
    vector<SegmentPair> edgesToBeAddedLong;


    // To make sure the SearchGraph and the Graph are strand-symmatric,
    // we only keep a SegmentPair if the minimum segment id
    // is less than the minimum segment id of the reverse complemented pair.
    // Then, when generating edges, we generate a pair of reverse
    // complemented edges for each segmentPair.

    // Loop over batches of segment pairs assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over segment pairs in this batch.
        for(uint64_t i=begin; i<end; i++) {
            const auto& [segment0, segment1] = segmentPairs[i];
            SHASTA2_ASSERT(segment0 != segment1);

            // See if this a SegmentPair that we want to use (to ensure strand symmetry).
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

            const SegmentPairInformation segmentPairInformation = SegmentStepSupport::analyzeSegmentPair(
                html, assemblyGraph, segment0, segment1, representativeRegionStepCount);

            // If the offset is negative, discard it.
            if(segmentPairInformation.segmentOffset < 0) {
                continue;
            }

            // Check for long segments.
            const bool isLong0 = isLong(segment0);
            const bool isLong1 = isLong(segment1);

            // Tentatively create the SegmentPair without storing it.
            const SegmentPair segmentPair(segment0, segment1, segmentPairInformation);

            // If it does not satisfy our requirements, discard it.
            if(segmentPair.segmentPairInformation.commonCount < minCommonCount) {
                continue;
            }
            if(
                (segmentPair.logP < logPThreshold) and
                (segmentPair.logPForward < logPThreshold) and
                (segmentPair.logPBackward < logPThreshold)) {
                continue;
            }
            if(not assemblyGraph.canConnect(segment0, segment1)) {
                continue;
            }

            // See if we can add it to the SearchGraphs
            if(segmentPair.logP > logPThreshold) {
                edgesToBeAdded.emplace_back(segmentPair);
            }

            // See if we can add it to the long graph.
            if(isLong0 and isLong1) {
                if( (segmentPair.logP         > logPThreshold) or
                    (segmentPair.logPForward  > logPThreshold) or
                    (segmentPair.logPBackward > logPThreshold)) {
                    edgesToBeAddedLong.emplace_back(segmentPair);
                }
            }

        }
    }



    // Now grab the mutex and add the edges we found.
    std::lock_guard<std::mutex> lock(mutex);
    for(const SegmentPair& segmentPair: edgesToBeAdded) {

        SearchGraphEdge edge(
            segmentPair.segmentPairInformation.commonCount,
            segmentPair.segmentPairInformation.missing0,
            segmentPair.segmentPairInformation.missing1);

        // Add it to the SearchGraph.
        add_edge(
            searchGraph.vertexMap.at(segmentPair.segment0),
            searchGraph.vertexMap.at(segmentPair.segment1),
            edge,
            searchGraph);

        // Also add the reverse complemented edge.
        swap(edge.missingCount0, edge.missingCount1);
        const Segment segment0Rc = assemblyGraph[segmentPair.segment0].eRc;
        const Segment segment1Rc = assemblyGraph[segmentPair.segment1].eRc;
        add_edge(
            searchGraph.vertexMap.at(segment1Rc),
            searchGraph.vertexMap.at(segment0Rc),
            edge,
            searchGraph);

    }



    for(const SegmentPair& segmentPair: edgesToBeAddedLong) {

        ConnectGraphEdge edge(
            segmentPair.segmentPairInformation.commonCount,
            segmentPair.segmentPairInformation.missing0,
            segmentPair.segmentPairInformation.missing1);

        add_edge(
            graph.vertexMap.at(segmentPair.segment0),
            graph.vertexMap.at(segmentPair.segment1),
            edge,
            graph);

        // Also add the reverse complemented edge.
        swap(edge.directConnectInformation.missingCount0, edge.directConnectInformation.missingCount1);
        swap(edge.directConnectInformation.logPForward, edge.directConnectInformation.logPBackward);
        const Segment segment0Rc = assemblyGraph[segmentPair.segment0].eRc;
        const Segment segment1Rc = assemblyGraph[segmentPair.segment1].eRc;
        add_edge(
            graph.vertexMap.at(segment1Rc),
            graph.vertexMap.at(segment0Rc),
            edge,
            graph);
    }
}



// Prune removes all vertices that are not accessible from long
// vertices in both directions.
void SearchGraph::prune()
{
    SearchGraph& graph = *this;

    // Loop over both directions.
    array<std::set<vertex_descriptor>, 2> reachedVertices;
    vector<vertex_descriptor> neighbors;
    for(uint64_t direction=0; direction<2; direction++) {

        // Initialize the BFS.
        std::queue<vertex_descriptor> q;
        BGL_FORALL_VERTICES(v, graph, SearchGraph) {
            if(graph[v].isLong) {
                q.push(v);
                reachedVertices[direction].insert(v);
            }
        }

        // BFS loop in this direction.
        while(not q.empty()) {
            const vertex_descriptor v0 = q.front();
            q.pop();

            // Find the neighbors in this direction.
            neighbors.clear();
            if(direction == 0) {
                // Forward.
                BGL_FORALL_OUTEDGES(v0, e, graph, SearchGraph) {
                    const vertex_descriptor v1 = target(e, graph);
                    neighbors.push_back(v1);
                }
            } else {
                // Backward.
                BGL_FORALL_INEDGES(v0, e, graph, SearchGraph) {
                    const vertex_descriptor v1 = source(e, graph);
                    neighbors.push_back(v1);
                }
            }

            // Loop over the neighbors.
            for(const vertex_descriptor v1: neighbors) {
                if(not reachedVertices[direction].contains(v1)) {
                    reachedVertices[direction].insert(v1);
                    q.push(v1);
                }
            }
        }
    }


    // Remove vertices that are not reachable in both directions.
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        if(not (reachedVertices[0].contains(v) and reachedVertices[1].contains(v))) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        vertexMap.erase(graph[v].segment);
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

    // Sanity check: all leafs must be long vertices.
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        const bool isLeaf = (in_degree(v, graph) == 0) or (out_degree(v, graph) == 0);
        if(isLeaf) {
            SHASTA2_ASSERT(graph[v].isLong);
        }
    }


}



void SearchGraph::writeGraphviz(
    const AssemblyGraph& assemblyGraph,
    const string& name) const
{
    const SearchGraph& graph = *this;

    ofstream dot("ReadFollowing-SearchGraph-" + name + ".dot");
    dot << "digraph ReadFollowingSearchGraph {\n";

    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        const SearchGraphVertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id;

        // Begin attributes.
        dot << " [";

        // Label.
        dot <<
            "label=\"" << assemblyGraphEdge.id << "\\n" <<
            vertex.length << "\\n" <<
            "\"";


        // Color.
        string color;
        if(vertex.isLong) {
            color = "cyan";
        }
        if(not color.empty()) {
            dot << " style=filled fillcolor=" << color;
        }

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    BGL_FORALL_EDGES(e, graph, SearchGraph) {
        const SearchGraphEdge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;

        dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

        // Begin attributes.
        dot << "[";

        // Tooltip.
        dot << " tooltip=\"" <<
            edge.commonCount << "/" <<
            edge.missingCount0 << "/" <<
            edge.missingCount1 << "/" <<
            std::fixed << std::setprecision(1) <<
            edge.logP << "\"";

        // Thickness is determined to logP.
        // Color is always black.
        double logPClipped = max(1., edge.logP);
        logPClipped = min(100., logPClipped);
        const double thickness = 0.05 * logPClipped;
        dot << std::fixed << std::setprecision(2) << " penwidth=" << thickness;

        // End attributes.
        dot << "]";

        // End the line for this edge.
        dot << ";\n";

    }

    dot << "}\n";
}



void ConnectGraph::writeGraphviz(
    const AssemblyGraph& assemblyGraph,
    const string& name) const
{
    const ConnectGraph& graph = *this;

    ofstream dot("ReadFollowing-ConnectGraph-" + name + ".dot");
    dot <<
        "digraph ReadFollowingConnectGraph {\n"
        "tooltip=\" \";\n";

    BGL_FORALL_VERTICES(v, graph, ConnectGraph) {
        const ConnectGraphVertex& vertex = graph[v];
        const Segment segment = vertex.segment;
        const AssemblyGraphEdge& assemblyGraphEdge = assemblyGraph[segment];
        dot << assemblyGraphEdge.id;

        // Begin attributes.
        dot << " [";

        // Label.
        dot << "label=\"" << assemblyGraphEdge.id << "\\n" << vertex.length;
        if(not assemblyGraphEdge.annotation.empty()) {
            dot << "\\n" << assemblyGraphEdge.annotation;
        }
        dot << "\"";

        // End attributes.
        dot << "]";

        // End the line for this vertex.
        dot << ";\n";
    }



    // Each edge is displayed as one or two display edges:
    // - If edge.hasDirectConnection(), an edge representing the
    //   DirectConnectInformation is drawn. This edge uses a solid line.
    //   The edge color depends on edge.directConnectionType():
    //   * If Bidirectional, the edge is black.
    //   * If forward, the edge is blue.
    //   * If backward, the edge is green.
    //   * If ambiguous, the edge is red.
    //   The arrows can be empty or filled:
    //   * The source arrow is filled if directConnectInformation.logPForward  >= logPThreshold.
    //   * The target arrow is filled if directConnectInformation.logPBackward >= logPThreshold.
    // - If one or both of the edge assemblyPaths are non-empty, an edge representing
    //   these assembly paths is drawn. This edge uses a dashed line.
    //   The edge color depends on the two assembly paths:
    //   * If both assembly paths are non-trivial, the edge is black.
    //   * If only the forward  assembly path is non-trivial, the edge is blue.
    //   * If only the backward assembly path is non-trivial, the edge is green.
    //   The arrows can be empty or filled:
    //   * The source arrow is filled if the forward  assembly path is non-trivial.
    //   * The target arrow is filled if the backward assembly path is non-trivial.

    BGL_FORALL_EDGES(e, graph, ConnectGraph) {
        const ConnectGraphEdge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Segment segment0 = graph[v0].segment;
        const Segment segment1 = graph[v1].segment;


        // Draw the edge representing the DirectConnectInformation.
        if(edge.hasDirectConnection()) {

            dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

            // Begin attributes.
            dot << "[";

            // Tooltip.
            dot << " tooltip=\"" <<
                edge.directConnectInformation.commonCount << "/" <<
                edge.directConnectInformation.missingCount0 << "/" <<
                edge.directConnectInformation.missingCount1 << "/" <<
                std::fixed << std::setprecision(1) <<
                edge.directConnectInformation.logP << "/" <<
                edge.directConnectInformation.logPForward << "/" <<
                edge.directConnectInformation.logPBackward << "\"";

            // Thickness is determined to maxLogP.
            double logPClipped = max(1., edge.directConnectInformation.maxLogP());
            logPClipped = min(100., logPClipped);
            const double thickness = 0.05 * logPClipped;
            dot << std::fixed << std::setprecision(2) << " penwidth=" << thickness;

            // Color depends on the edge type.
            string color;
            switch(edge.directConnectionType()) {
            case ConnectGraphEdge::DirectConnectionType::None:
                color = "Red";
                break;
            case ConnectGraphEdge::DirectConnectionType::Bidirectional:
                color = "Black";
                break;
            case ConnectGraphEdge::DirectConnectionType::Forward:
                color = "Blue";
                break;
            case ConnectGraphEdge::DirectConnectionType::Backward:
                color = "Green";
                break;
            case ConnectGraphEdge::DirectConnectionType::Ambiguous:
                color = "Orange";
                break;
            default:
                SHASTA2_ASSERT(0);
            }
            dot << " color=" << color;

            // The source arrow is filled if directConnectInformation.logPForward  >= logPThreshold.
            // The target arrow is filled if directConnectInformation.logPBackward >= logPThreshold.
            dot << " dir=both";
            if(edge.directConnectInformation.logPForward  >= logPThreshold) {
                dot << " arrowtail=inv";
            } else {
                dot << " arrowtail=oinv";
            }
            if(edge.directConnectInformation.logPBackward  >= logPThreshold) {
                dot << " arrowhead=normal";
            } else {
                dot << " arrowhead=onormal";
            }

            // End attributes.
            dot << "]";

            // End the line for this edge.
            dot << ";\n";
        }



        // Draw the edge representing the assembly paths.
        if((edge.assemblyPaths[0].size() > 2) or (edge.assemblyPaths[1].size()  > 2)) {

            dot << assemblyGraph[segment0].id << "->" << assemblyGraph[segment1].id;

            // Begin attributes.
            dot << "[";

            dot << "style=dashed";

            // Tooltip.
            dot << " tooltip=\"" <<
                edge.assemblyPaths[0].size() << "/" <<
                edge.assemblyPaths[1].size() << "\"";

            // The arrowtail is filled if the forward path exists.
            // The arrowhead is filled if the backward path exists.
            dot << " dir=both";
            if(edge.assemblyPaths[0].size() <= 2) {
                dot << " arrowtail=oinv";
            } else {
                dot << " arrowtail=inv";
            }
            if(edge.assemblyPaths[1].size() <= 2) {
                dot << " arrowhead=onormal";
            } else {
                dot << " arrowhead=normal";
            }

            // Color.
            string color = "Black";
            if(edge.assemblyPaths[0].size() <= 2) {
                color = "Green";
            }
            if(edge.assemblyPaths[1].size() <= 2) {
                color = "Blue";
            }
            dot << " color=" << color;


            // End attributes.
            dot << "]";

            // End the line for this edge.
            dot << ";\n";
        }

    }

    dot << "}\n";
}



// This finds a shortest path starting at segment0 and ending at a long Segment,
// with path length defined by Edge::weight.
void ReadFollower::findShortestPath(
    Segment segment0,
    uint64_t direction,     // 0 = forward, 1 = backward
    vector<Segment>& path
    ) const
{
    if(direction == 0) {
        findShortestPathForward(segment0, path);
    } else {
        findShortestPathBackward(segment0, path);
    }
}



void ReadFollower::findShortestPathForward(
    Segment segment0,
    vector<Segment>& path
    ) const
{
    searchGraph.findShortestPath(segment0, path);
}



// This does a forward search in the searchGraph starting at the reverse complement
// of segment0, then returns the reverse complement of the path found in this way.
void ReadFollower::findShortestPathBackward(
    Segment segment0,
    vector<Segment>& path
    ) const
{
    const Segment segment0Rc = assemblyGraph[segment0].eRc;

    vector<Segment> pathRc;
    findShortestPathForward(segment0Rc, pathRc);

    path.clear();
    for(Segment segmentRc: pathRc | std::views::reverse) {
        path.push_back(assemblyGraph[segmentRc].eRc);
    }
}



void ReadFollower::findAndWriteShortestPath(Segment segment0, uint64_t direction) const
{

    vector<Segment> path;

    if(direction == 0) {
        findShortestPath(segment0, direction, path);
    } else {
        findShortestPathBackward(segment0, path);
    }

    cout << "Found a path of length " << path.size() << ":" << endl;
    for(const Segment segment: path) {
        cout << assemblyGraph[segment].id << ",";
    }
    cout << endl;

}



void SearchGraph::findShortestPath(Segment segment0, vector<Segment>& path) const
{
    using namespace boost;
    const SearchGraph& graph = *this;

    // Find the vertex corresponding to this segment.
    const vertex_descriptor v0 = vertexMap.at(segment0);

    // If no outgoing edges, return a path consisting of just segment0.
    if(out_degree(v0, graph) == 0) {
        path.clear();
        path.push_back(segment0);
        return;
    }

    // An exception class used to stop the shortest path computation
    // when a long vertex is encountered.
    class LongVertexReached {
    public:
        vertex_descriptor v;
        LongVertexReached(vertex_descriptor v) : v(v) {}
    };

    // The DijkstraVisitor class throws LongVertexReached when a long vertex is encountered.
    class DijkstraVisitor : public boost::dijkstra_visitor<> {
    public:
        vertex_descriptor v0;
        DijkstraVisitor(vertex_descriptor v0) : v0(v0) {}
        void examine_vertex(vertex_descriptor v, const SearchGraph& graph)
        {
            if((v != v0) and (graph[v].isLong)) {
                throw LongVertexReached(v);
            }
        }
    };
    DijkstraVisitor dijkstraVisitor(v0);

    // The predecessorMap is filled in by the call to dag_shortest_paths
    // and can be used to reconstruct the path from v0 to
    // the first long edge encountered.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;

    // Compute the shortest path using Edge::weight.
    vertex_descriptor v1 = null_vertex();
    try {
        dijkstra_shortest_paths(graph, v0,
           weight_map(boost::get(&SearchGraphEdge::weight, graph)).
           vertex_index_map(make_assoc_property_map(vertexIndexMap)).
           predecessor_map(make_assoc_property_map(predecessorMap)).
           visitor(dijkstraVisitor)
           );
    } catch(LongVertexReached& longVertexReached) {
        v1 = longVertexReached.v;
    } catch(std::exception& e) {
        SHASTA2_ASSERT(0);
    }

    // Use the predecessor map to construct the path.
    path.clear();
    vertex_descriptor v = v1;
    while(true) {
        path.push_back(graph[v].segment);
        if(v == v0) {
            break;
        }
        v = predecessorMap[v];
    }
    std::ranges::reverse(path);
}



void SearchGraph::createVertexIndexMap()
{
    SearchGraph& graph = *this;

    vertexIndexMap.clear();
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, SearchGraph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }
}



SearchGraphEdge::SearchGraphEdge(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    commonCount(commonCount),
    missingCount0(missingCount0),
    missingCount1(missingCount1)
{
    logP         = a * double(commonCount) - b * double(missingCount0 + missingCount1);
    weight = pow(10., 0.1 * logP);
}



ConnectGraphEdge::ConnectGraphEdge(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    directConnectInformation(commonCount, missingCount0, missingCount1)
{}



ConnectGraphEdge::DirectConnectInformation::DirectConnectInformation(
    uint64_t commonCount,
    uint64_t missingCount0,
    uint64_t missingCount1) :
    commonCount(commonCount),
    missingCount0(missingCount0),
    missingCount1(missingCount1)
{
    logP         = a * double(commonCount) - b * double(missingCount0 + missingCount1);
    logPForward  = a * double(commonCount) - b * double(missingCount0);
    logPBackward = a * double(commonCount) - b * double(missingCount1);
}



bool ReadFollower::isLong(Segment segment) const
{
    return graph.vertexMap.contains(segment);
}



double ConnectGraphEdge::DirectConnectInformation::maxLogP() const
{
    return max(logP, max(logPForward, logPBackward));
}



ConnectGraphEdge::DirectConnectionType ConnectGraphEdge::directConnectionType() const
{
    if(not hasDirectConnection()) {
        return DirectConnectionType::None;
    }

    if(directConnectInformation.logP >= logPThreshold) {
        return DirectConnectionType::Bidirectional;
    }

    if(directConnectInformation.logPForward >= logPThreshold) {
        if(directConnectInformation.logPBackward < logPThreshold) {
            return DirectConnectionType::Forward;
        } else {
            return DirectConnectionType::Ambiguous;
        }
    }

    if(directConnectInformation.logPBackward >= logPThreshold) {
        if(directConnectInformation.logPForward < logPThreshold) {
            return DirectConnectionType::Backward;
        } else {
            return DirectConnectionType::Ambiguous;
        }
    }

    return DirectConnectionType::Ambiguous;
}



// Use the SearchGraphs to find shortest paths between long segments
// and store them in the ConnectGraph.
void ReadFollower::findShortestPathsMultithreaded(uint64_t threadCount)
{
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
    FindShortestPathsData& data = findShortestPathsData;

    // Store a vector with all the ConnectGraph vertices to be processed.
    data.graphVertices.clear();
    BGL_FORALL_VERTICES(v, graph, ConnectGraph) {
        data.graphVertices.push_back(v);
    }
    performanceLog << timestamp << "ReadFollower::findShortestPathsMultithreaded begins with " <<
        data.graphVertices.size() << " long segments." << endl;

    // Process them in parallel.
    const uint64_t batchSize = 1;
    setupLoadBalancing(data.graphVertices.size(), batchSize);
    runThreads(&ReadFollower::findShortestPathsThreadFunction, threadCount);

    performanceLog << timestamp << "ReadFollower::findShortestPathsMultithreaded ends." << endl;
}



void ReadFollower::findShortestPathsThreadFunction(uint64_t)
{
    FindShortestPathsData& data = findShortestPathsData;
    vector<Segment> path;

    // Loop over batches assigned to this thread.
    uint64_t begin, end;
    while(getNextBatch(begin, end)) {

        // Loop over vertices in this batch.
        for(uint64_t i=begin; i!=end; i++) {
            if((i % 1000) == 0) {
                std::lock_guard<std::mutex> lock(mutex);
                performanceLog << timestamp << "Read following: " << i << "/" << data.graphVertices.size() << endl;
            }
            const ConnectGraph::vertex_descriptor v0 = data.graphVertices[i];
            const Segment segment0 = graph[v0].segment;

            // Loop for shortest paths in both directions.
            for(uint64_t direction=0; direction<2; direction++) {
                findShortestPath(segment0, direction, path);

                // Discard a trivial path.
                if(path.size() < 2) {
                    continue;
                }

                const Segment s0 = path.front();
                const Segment s1 = path.back();
                if(direction == 0) {
                    SHASTA2_ASSERT(s0 == segment0);
                } else {
                    SHASTA2_ASSERT(s1 == segment0);
                }


                // For all operations on the ConnectGraph we need to acquire the mutex.
                std::lock_guard<std::mutex> lock(mutex);

                // Now we can update the ConnectGraph.
                const ConnectGraph::vertex_descriptor u0 = graph.vertexMap.at(s0);
                const ConnectGraph::vertex_descriptor u1 = graph.vertexMap.at(s1);

                // Look for a ConnectGraph edge between the u0 and u1.
                ConnectGraph::edge_descriptor e;
                bool edgeExists;
                tie(e, edgeExists) = boost::edge(u0, u1, graph);
                if(not edgeExists) {
                    tie(e, edgeExists) = boost::add_edge(u0, u1, graph);
                }
                SHASTA2_ASSERT(edgeExists);

                // Store the path.
                graph[e].assemblyPaths[direction] = path;
            }
        }
    }
}



// This removes edges without a direct connection
// and that don't have paths in both directions.
void ConnectGraph::removeWeakEdges()
{
    ConnectGraph& graph = *this;

    vector<edge_descriptor> edgesToBeRemoved;
    BGL_FORALL_EDGES(e, graph, ConnectGraph) {
        const ConnectGraphEdge& edge = graph[e];

        // If it has a direct connection, keep it.
        if(edge.hasDirectConnection()) {
            continue;
        }

        // Otherwise, only keep it if both assembly paths are present.
        if(edge.assemblyPaths[0].empty() or edge.assemblyPaths[1].empty()) {
            edgesToBeRemoved.push_back(e);
        }
    }


    for(const edge_descriptor e: edgesToBeRemoved) {
        boost::remove_edge(e, graph);
    }
}



void ConnectGraph::transitiveReduction()
{
    ConnectGraph& graph = *this;

    // Gather all edges. They will be processed in this order.
    // It may be better to use a better ordering.
    vector<edge_descriptor> allEdges;
    BGL_FORALL_EDGES(e, graph, ConnectGraph) {
        allEdges.push_back(e);
    }

    // Transitive reduction.
    for(const edge_descriptor e: allEdges) {
        if(transitiveReductionCanRemove(e)) {
            boost::remove_edge(e, graph);
        }
    }

}



bool ConnectGraph::transitiveReductionCanRemove(edge_descriptor e) const
{
    const ConnectGraph& graph = *this;
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);

    // Do a forward BFS starting at v0.
    // If we encounter v1, return true.
    std::queue<vertex_descriptor> q;
    q.push(v0);

    // A map to store vertices already encountered and their distance from v0.
    std::map<vertex_descriptor, uint64_t> m;
    m.insert(make_pair(v0, 0));



    // Main BFS loop.
    while(not q.empty()) {

        // Dequeue a vertex.
        const vertex_descriptor vA = q.front();
        q.pop();
        const auto itA = m.find(vA);
        SHASTA2_ASSERT(itA != m.end());
        const uint64_t distanceA = itA->second;
        const uint64_t distanceB = distanceA + 1;

        // Loop over its out-edges, without using edge e.
        BGL_FORALL_OUTEDGES(vA, eAB, graph, ConnectGraph) {
            if(eAB == e) {
                continue;
            }

            // If we reached v1, return true;
            const vertex_descriptor vB = target(eAB, graph);
            if(vB == v1) {
                return true;
            }

            // If we already encountered vB, don't do anything.
            if(m.contains(vB)) {
                continue;
            }

            if(distanceB < transitiveReductionMaxDistance) {
                q.push(vB);
                m.insert(make_pair(vB, distanceB));
            }
        }
    }

    // If getting here we did not encounter v1 in the BFS loop.
    return false;
}



// Get the assembly path to be used for this edge.
// This includes the source and target Segment.
vector<Segment> ConnectGraph::getAssemblyPath(edge_descriptor e) const
{
    const ConnectGraph& graph = *this;
    const ConnectGraphEdge& edge = graph[e];
    const vertex_descriptor v0 = source(e, graph);
    const vertex_descriptor v1 = target(e, graph);
    const Segment segment0 = graph[v0].segment;
    const Segment segment1 = graph[v1].segment;

    // If both assemblyPaths are non-trivial, return the longer one.
    const uint64_t length0 = edge.assemblyPaths[0].size();
    const uint64_t length1 = edge.assemblyPaths[1].size();
    if((length0 > 2) and (length1 > 2)) {
        if(length0 > length1) {
            return edge.assemblyPaths[0];
        } else {
            return edge.assemblyPaths[1];
        }
    }

    // If the DirectConnectInformation exists and is bidirectional,
    // return an assembly path consisting of just the source and target segments.
    if(edge.hasDirectConnection() and edge.directConnectionType() == ConnectGraphEdge::DirectConnectionType::Bidirectional) {
        return {segment0, segment1};
    }

    // If one of the assembly paths is non-trivial, return it.
    if(length0 > 2) {
        return edge.assemblyPaths[0];
    }
    if(length1 > 2) {
        return edge.assemblyPaths[1];
    }

    // In all other cases, return an assembly path consisting of just the source and target segments.
    return {segment0, segment1};
}



// Use the ConnectGraph to update the AssemblyGraph.
void ReadFollower::updateAssemblyGraph(AssemblyGraph& assemblyGraph) const
{
    const bool debug = false;
    writeMemoryStatistics("ReadFollower::updateAssemblyGraph begin");

    // Create a disconnected version of each long Segment.
    std::map<Segment, Segment> longSegmentMap; // (oldSegment, newSegment) (They have the same id).
    BGL_FORALL_VERTICES(v, graph, ConnectGraph) {
        const Segment oldSegment = graph[v].segment;
        const Segment newSegment = createDisconnectedSegmentCopy(assemblyGraph, oldSegment);
        longSegmentMap.insert({oldSegment, newSegment});
    }



    // A short Segment usually appears in only one assembly path.
    // In that case, when making a disconnected copy of that Segment
    // we keep the same id.
    // However, occasionally a short Segment will appear in more than one
    // assembly path. In that case, the first copy keeps the same id,
    // but subsequent copies are given new ids.
    // This set keeps track of the short segments we already used.
    std::set<Segment> usedShortSegments;



    // Each edge of the ConnectGraph generates a linear chain between
    // the source and target segments.
    BGL_FORALL_EDGES(e, graph, ConnectGraph) {

        // Get the old Segments.
        const ConnectGraph::vertex_descriptor oldV0 = source(e, graph);
        const ConnectGraph::vertex_descriptor oldV1 = target(e, graph);
        const Segment oldSegment0 = graph[oldV0].segment;
        const Segment oldSegment1 = graph[oldV1].segment;

        // Get the corresponding new Segments.
        const Segment newSegment0 = longSegmentMap.at(oldSegment0);
        const Segment newSegment1 = longSegmentMap.at(oldSegment1);

        // Get the assembly path for this ConnectGraphEdge.
        const vector<Segment> oldAssemblyPath = graph.getAssemblyPath(e);

        if(debug) {
            cout << "Assembly path to connect " <<
                assemblyGraph[oldSegment0].id << " with " <<
                assemblyGraph[oldSegment1].id << ":" << endl;
            for(const Segment segment: oldAssemblyPath) {
                cout << assemblyGraph[segment].id << " ";
            }
            cout << endl;
        }

        // Sanity checks: the assembly path begins/ends at segment0/segment1.
        SHASTA2_ASSERT(oldAssemblyPath.size() >= 2);
        SHASTA2_ASSERT(oldAssemblyPath.front() == oldSegment0);
        SHASTA2_ASSERT(oldAssemblyPath.back() == oldSegment1);

        // Sanity check: the segments internal to the assembly path
        // are not long segment.
        for(uint64_t i=1; i<oldAssemblyPath.size()-1; i++) {
            const Segment oldSegment = oldAssemblyPath[i];
            SHASTA2_ASSERT(not longSegmentMap.contains(oldSegment));
        }

        // Generate the new Segments of this assembly path.
        vector<Segment> newAssemblyPath;
        newAssemblyPath.push_back(newSegment0);
        for(uint64_t i=1; i<oldAssemblyPath.size()-1; i++) {
            const Segment oldSegment = oldAssemblyPath[i];
            const Segment newSegment = createDisconnectedSegmentCopy(assemblyGraph, oldSegment);
            newAssemblyPath.push_back(newSegment);
            if(usedShortSegments.contains(oldSegment)) {
                assemblyGraph[newSegment].id = assemblyGraph.nextEdgeId++;
                cout << "Additional copy of " << assemblyGraph[oldSegment].id <<
                    " takes new id " << assemblyGraph[newSegment].id << endl;
            } else {
                usedShortSegments.insert(oldSegment);
            }
        }
        newAssemblyPath.push_back(newSegment1);



        // Now we have to connect adjacent segments.
        for(uint64_t i1=1; i1<newAssemblyPath.size(); i1++) {
            const uint64_t i0 = i1 - 1;
            const Segment newSegment0 = newAssemblyPath[i0];
            const Segment newSegment1 = newAssemblyPath[i1];

            const AssemblyGraph::vertex_descriptor newV0 = target(newSegment0, assemblyGraph);
            const AssemblyGraph::vertex_descriptor newV1 = source(newSegment1, assemblyGraph);

            const AnchorId anchorId0 = assemblyGraph[newV0].anchorId;
            const AnchorId anchorId1 = assemblyGraph[newV1].anchorId;

            // Create the new edge.
            // If the two anchors are the same, leave it empty without any steps.
            // Otherwise use the same process in Tangle1::addConnectPair.
            Segment newSegment;
            tie(newSegment, ignore) = add_edge(newV0, newV1, AssemblyGraphEdge(assemblyGraph.nextEdgeId++), assemblyGraph);
            AssemblyGraphEdge& newEdge = assemblyGraph[newSegment];
            if(anchorId0 != anchorId1) {

                // Create the RestrictedAnchorGraph, then:
                // - Remove vertices not accessible from anchorId0 and anchorId1.
                // - Remove cycles.
                // - Find the longest path.
                // - Add one step for each edge of the longest path of the RestrictedAnchorGraph.

                ostream html(0);
                const TangleMatrix1 tangleMatrix(
                    assemblyGraph,
                    vector<Segment>(1, newSegment0),
                    vector<Segment>(1, newSegment1),
                    html);

                try {
                    RestrictedAnchorGraph restrictedAnchorGraph(assemblyGraph.anchors, assemblyGraph.journeys, tangleMatrix, 0, 0, html);
                    vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                    // restrictedAnchorGraph.findLongestPath(longestPath);
                    restrictedAnchorGraph.findOptimalPath(anchorId0, anchorId1, longestPath);

                    for(const RestrictedAnchorGraph::edge_descriptor re: longestPath) {
                        const auto& rEdge = restrictedAnchorGraph[re];
                        if(rEdge.anchorPair.size() == 0) {
                            newEdge.clear();
                            SHASTA2_ASSERT(0);
                        }
                        newEdge.push_back(AssemblyGraphEdgeStep(rEdge.anchorPair, rEdge.offset));
                    }
                } catch(RestrictedAnchorGraph::NoTransitions&) {
                    cout << "Could not connect " << assemblyGraph[newSegment0].id <<
                        " with " << assemblyGraph[newSegment1].id << endl;
                    SHASTA2_ASSERT(0);
                }
            }
        }

    }



    // Remove the old copy of each long Segment.
    for(const auto& [oldSegment, newSegment]: longSegmentMap) {
        boost::remove_edge(oldSegment, assemblyGraph);
    }

    // Remove the old copy of short Segments we used.
    for(const Segment oldSegment: usedShortSegments) {
        boost::remove_edge(oldSegment, assemblyGraph);
    }

    writeMemoryStatistics("ReadFollower::updateAssemblyGraph end");
}



// Make a disconnected copy of a Segment. The copy keeps the same id,
// so the original Segment will have to be removed for the AssemblyGraph
// to remain valid.
Segment ReadFollower::createDisconnectedSegmentCopy(
    AssemblyGraph& assemblyGraph,
    Segment oldSegment) const
{

    // Get the assembly graph edge and vertices of the oldSegment.
    const AssemblyGraphEdge& oldEdge = assemblyGraph[oldSegment];
    const AssemblyGraph::vertex_descriptor oldV0 = source(oldSegment, assemblyGraph);
    const AssemblyGraph::vertex_descriptor oldV1 = target(oldSegment, assemblyGraph);
    const AssemblyGraphVertex& oldVertex0 = assemblyGraph[oldV0];
    const AssemblyGraphVertex& oldVertex1 = assemblyGraph[oldV1];

    // Create the new vertices.
    AssemblyGraphVertex newVertex0 = oldVertex0;
    AssemblyGraphVertex newVertex1 = oldVertex1;
    newVertex0.id = assemblyGraph.nextVertexId++;
    newVertex1.id = assemblyGraph.nextVertexId++;
    const AssemblyGraph::vertex_descriptor newV0 = add_vertex(newVertex0, assemblyGraph);
    const AssemblyGraph::vertex_descriptor newV1 = add_vertex(newVertex1, assemblyGraph);

    // Create the edge for the new Segment, keeping the same id.
    Segment newSegment;
    tie(newSegment, ignore)= add_edge(newV0, newV1, oldEdge, assemblyGraph);

    return newSegment;
}



// Checks that the SearchGraph is strand-symmetric.
void SearchGraph::check(const AssemblyGraph& assemblyGraph) const
{
    cout << "SearchGraph::check begins." << endl;
    const SearchGraph& searchGraph = *this;

    // For strand symatry, every edge must have an identical reverse complemented edge.
    BGL_FORALL_EDGES(e, searchGraph, SearchGraph) {

        // Get the segments of this edge.
        const vertex_descriptor v0 = source(e, searchGraph);
        const vertex_descriptor v1 = target(e, searchGraph);
        const Segment segment0 = searchGraph[v0].segment;
        const Segment segment1 = searchGraph[v1].segment;

        // Get the reverse complemented segments.
        const Segment segment0Rc = assemblyGraph[segment0].eRc;
        const Segment segment1Rc = assemblyGraph[segment1].eRc;
        SHASTA2_ASSERT(segment0Rc != assemblyGraphNullEdge);
        SHASTA2_ASSERT(segment1Rc != assemblyGraphNullEdge);

        // Get the corresponding SearchGraph vertices.
        const vertex_descriptor v0Rc = vertexMap.at(segment0Rc);
        const vertex_descriptor v1Rc = vertexMap.at(segment1Rc);

        // Get the reverse complented edge.
        auto[eRc, edgeExists] = edge(v1Rc, v0Rc, searchGraph);
        SHASTA2_ASSERT(edgeExists);

        // Check that they are the reverse complement of each other.
#if 0
        cout << "Checking " <<
            assemblyGraph[segment0].id << " " << assemblyGraph[segment1].id << " " <<
            assemblyGraph[segment0Rc].id << " " << assemblyGraph[segment1Rc].id << endl;
        if(not(searchGraph[e].isReverseComplement(searchGraph[eRc]))) {
            cout << "Strand symmetry check." << endl;
            cout << assemblyGraph[segment0].id << " " << assemblyGraph[segment1].id << " " <<
                searchGraph[e].commonCount << " " <<
                searchGraph[e].missingCount0 << " " <<
                searchGraph[e].missingCount1 << " " <<
                searchGraph[e].logP << " " <<
                searchGraph[e].weight << endl;
            cout << assemblyGraph[segment0Rc].id << " " << assemblyGraph[segment1Rc].id << " " <<
                searchGraph[eRc].commonCount << " " <<
                searchGraph[eRc].missingCount0 << " " <<
                searchGraph[eRc].missingCount1 << " " <<
                searchGraph[eRc].logP << " " <<
                searchGraph[eRc].weight << endl;
        }
#endif
        SHASTA2_ASSERT(searchGraph[e].isReverseComplement(searchGraph[eRc]));
    }

    cout << "SearchGraph::check ends." << endl;
}


bool SearchGraphEdge::isReverseComplement(const SearchGraphEdge& that) const
{
    if(commonCount != that.commonCount) {
        return false;
    }
    if(missingCount0 != that.missingCount1) {
        return false;
    }
    if(missingCount1 != that.missingCount0) {
        return false;
    }
    if(logP != that.logP) {
        return false;
    }
    if(weight != that.weight) {
        return false;
    }

    return true;
}



// Check that it is strand-symmetric.
void ConnectGraph::check(const AssemblyGraph& assemblyGraph) const
{
    cout << "ConnectGraph::check begins." << endl;
    const ConnectGraph& connectGraph = *this;

    // Every edge must have a corresponding reverse complemented edge.
    BGL_FORALL_EDGES(e, connectGraph, ConnectGraph) {

        // Find the Segments of this edge.
        const vertex_descriptor v0 = source(e, connectGraph);
        const vertex_descriptor v1 = target(e, connectGraph);
        const Segment segment0 = connectGraph[v0].segment;
        const Segment segment1 = connectGraph[v1].segment;

        // Find the reverse complemented segments and
        // the corresponding vertices.
        const Segment segment0Rc = assemblyGraph[segment0].eRc;
        const Segment segment1Rc = assemblyGraph[segment1].eRc;
        const vertex_descriptor v0Rc = vertexMap.at(segment0Rc);
        const vertex_descriptor v1Rc = vertexMap.at(segment1Rc);

        // Find the reverse complemented edge.
        auto[eRc, edgeExists] = edge(v1Rc, v0Rc, connectGraph);
        SHASTA2_ASSERT(edgeExists);
    }
    cout << "ConnectGraph::check ends." << endl;
}

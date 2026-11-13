#include "AssemblyGraph.hpp"
#include "deduplicate.hpp"
#include "findDescendants.hpp"
#include "Options.hpp"
#include "performanceLog.hpp"
#include "ReadFollowing5.hpp"
#include "ReadFollowing6.hpp"
#include "StrandSeparation1.hpp"
#include "SegmentGraph.hpp"
#include "Tangle.hpp"
#include "timestamp.hpp"
using namespace shasta2;
using namespace ReadFollowing5;
using namespace ReadFollowing6;



// This uses ReadFollowing5 code.
bool AssemblyGraph::readFollowingStrandSymmetric5(
    uint64_t tangleId,
    const Tangle& tangle,
    ostream& html,
    const string& debugOutputBaseName)
{
    SHASTA2_ASSERT(not tangle.entrances.empty());
    SHASTA2_ASSERT(not tangle.exits.empty());

    ReadFollowing5::Graph graph(*this, tangleId, tangle, html, debugOutputBaseName);

    if(graph.updateAssemblyGraph1()) {
        return true;
    } else {
        return graph.updateAssemblyGraph2();
    }
}



void AssemblyGraph::separateStrands(const string& debugOutputBaseName)
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t maxLengthStart = 10000;
    const uint64_t maxLengthEnd = 100000000;

    for(uint64_t maxLength=maxLengthStart; maxLength<=maxLengthEnd; maxLength*=2) {
        separateStrands1(maxLength, debugOutputBaseName);    
    }

    separateStrands1Global(debugOutputBaseName);
}



void AssemblyGraph::separateStrands1Global(const string& debugOutputBaseName)
{
    const uint64_t maxLength = std::numeric_limits<uint64_t>::max();
    separateStrands1(maxLength, debugOutputBaseName);    
}



void AssemblyGraph::separateStrands1(
    uint64_t maxLength, 
    const string& debugOutputBaseName)
{
    AssemblyGraph& assemblyGraph = *this;

    SHASTA2_ASSERT(countZeroLengthSegments() == 0);

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 begins: " <<
        debugOutputBaseName << endl;

    // Create tangles as connected components of the AssemblyGraph,
    // computed considering only short Segments.
    // The vertices of each tangle are returned sorted by id
    // so we can do binary searches in them.
    vector< vector<vertex_descriptor> > tangles;
    vector<uint64_t> tangleRc;
    createTanglesBySegmentLength(maxLength, tangles, tangleRc);

    // The self-complementary tangles are strand contacts.
    vector< vector<vertex_descriptor> > strandContacts;
    for(uint64_t tangleId=0; tangleId<tangles.size(); tangleId++) {
        if(tangleRc[tangleId] == tangleId) {
            strandContacts.push_back(tangles[tangleId]);
        }
    }

    if(maxLength == std::numeric_limits<uint64_t>::max()) {
        cout << "Strand separation at infinite maxLength " 
            " found " << strandContacts.size() << " strand contacts." << endl;
    } else {
        cout << "Strand separation at maxLength " << maxLength << 
            " found " << strandContacts.size() << " strand contacts." << endl;
    }

    // Write a csv file that can be imported into Bandage to see
    // the tangles.
    vector<uint64_t> strandContactsRc(strandContacts.size());
    std::iota(strandContactsRc.begin(), strandContactsRc.end(), 0);
    writeTangles(strandContacts, strandContactsRc,
        debugOutputBaseName + "-StrandContacts-Bandage.csv");

    // Process each strand contact separately.
    for(uint64_t strandContactId=0; strandContactId<strandContacts.size(); strandContactId++) {
        const StrandSeparation1::StrandContact strandContact(
            assemblyGraph,
            strandContacts[strandContactId],
            debugOutputBaseName,
            strandContactId);
    }

    removeZeroLengthSegmentsStrandSymmetric();
    removeIsolatedVertices();

    performanceLog << timestamp << "AssemblyGraph::separateStrands1 ends: " <<
        debugOutputBaseName << endl;
}



// This uses ReadFollowing6 code.
bool AssemblyGraph::readFollowingStrandSymmetric6(
    [[maybe_unused]] uint64_t tangleId,
    [[maybe_unused]] const shasta2::Tangle& shasta2Tangle,
    [[maybe_unused]] ostream& html,
    const string& debugOutputBaseName)
{
    SHASTA2_ASSERT(not shasta2Tangle.isSelfComplementary());

    ReadFollowing6::Tangle tangle(
        *this,
        shasta2Tangle.tangleVertices,
        html,
        tangleId,
        debugOutputBaseName);

    return false;
}


#if 0
void AssemblyGraph::findStrandContacts() const
{
    const AssemblyGraph& assemblyGraph = *this;
    const bool debug = true;

    // Create a graph with half the vertices and edges 
    // as the AssemblyGraph.
    // Each vertex corresponds to a pair of reverse complemented
    // vertices in the AssemblyGraph.
    // Each edge corresponds to a pair of reverse complemented
    // edges in the AssemblyGraph.
    // Names of vertices and edges of the AssemblyGraph
    // are prefixed with "a".

    class Vertex {
    public:
        // Ordered so id(v0) < id(v1).
        AssemblyGraph::vertex_descriptor av0;
        AssemblyGraph::vertex_descriptor av1;
    };

    class Edge {
    public:
        // Ordered so id(v0) < id(v1).
        AssemblyGraph::edge_descriptor ae0;
        AssemblyGraph::edge_descriptor ae1;
    };

    using Graph = boost::adjacency_list<
        boost::listS,
        boost::listS,
        boost::bidirectionalS,
        Vertex,
        Edge>;
    Graph graph;



    // Create the vertices.
    std::map<AssemblyGraph::vertex_descriptor, Graph::vertex_descriptor> vertexMap;
    BGL_FORALL_VERTICES(av, assemblyGraph, AssemblyGraph) {

        // Access this vertex and its reverse complement.
        const AssemblyGraphVertex aVertex = assemblyGraph[av];
        const AssemblyGraph::vertex_descriptor avRc = aVertex.vRc;
        SHASTA2_ASSERT(avRc != AssemblyGraph::null_vertex());
        const AssemblyGraphVertex aVertexRc = assemblyGraph[avRc];
        SHASTA2_ASSERT(aVertexRc.vRc == av);

        if(aVertex.id < aVertexRc.id) {
            const Graph::vertex_descriptor v = boost::add_vertex(Vertex({av, avRc}), graph);
            vertexMap.insert({av, v});
            vertexMap.insert({avRc, v});
        }
    }



    // Create the edges.
    std::map<AssemblyGraph::edge_descriptor, Graph::edge_descriptor> edgeMap;
    BGL_FORALL_EDGES(ae, assemblyGraph, AssemblyGraph) {

        // Access this edge and its reverse complement.
        const AssemblyGraphEdge aEdge = assemblyGraph[ae];
        const AssemblyGraph::edge_descriptor aeRc = aEdge.eRc;
        SHASTA2_ASSERT(aeRc != assemblyGraphNullEdge);
        const AssemblyGraphEdge aEdgeRc = assemblyGraph[aeRc];
        SHASTA2_ASSERT(aEdgeRc.eRc == ae);

        if(aEdge.id < aEdgeRc.id) {
            const AssemblyGraph::vertex_descriptor av0 = source(ae, assemblyGraph);
            const AssemblyGraph::vertex_descriptor av1 = target(ae, assemblyGraph);
            Graph::vertex_descriptor v0 = vertexMap.at(av0);
            Graph::vertex_descriptor v1 = vertexMap.at(av1);

            // Figure out if we need to swap v0 and v1.
            const AssemblyGraphVertex& aVertex0 = assemblyGraph[av0];
            const AssemblyGraphVertex& aVertex1 = assemblyGraph[av1];
            const AssemblyGraph::vertex_descriptor av0Rc = aVertex0.vRc;
            const AssemblyGraph::vertex_descriptor av1Rc = aVertex1.vRc;
            const AssemblyGraphVertex& aVertex0Rc = assemblyGraph[av0Rc];
            const AssemblyGraphVertex& aVertex1Rc = assemblyGraph[av1Rc];
            const bool swap0 = (aVertex0Rc.id < aVertex0.id);
            const bool swap1 = (aVertex1Rc.id < aVertex1.id);
            const bool swap01 = swap0 xor swap1;
            cout << "AAA " << aEdge.id << " " << int(swap0) << int(swap1) << int(swap01) << endl;
            if(swap01) {
                std::swap(v0, v1);
            }

            const auto[e, edgeWasAdded] = boost::add_edge(v0, v1, Edge({ae, aeRc}), graph);
            edgeMap.insert({ae, e});
            edgeMap.insert({aeRc, e});
        }
    }

    SHASTA2_ASSERT(num_vertices(assemblyGraph) == 2 * num_vertices(graph));
    SHASTA2_ASSERT(num_edges(assemblyGraph) == 2 * num_edges(graph));



    // Write the Graph in graphviz format.
    if(debug) {
        ofstream dot("FindStrandContacts.dot");
        dot << "digraph FindStrandContacts {\n";

        BGL_FORALL_VERTICES(v, graph, Graph) {
            const Vertex& vertex = graph[v];
            dot << id(vertex.av0) << "[label=\"" <<
                id(vertex.av0) << "\\n" << id(vertex.av1) <<
                "\"]\n";
        }

        BGL_FORALL_EDGES(e, graph, Graph) {
            const Edge& edge = graph[e];
            const vertex_descriptor v0 = source(e, graph);
            const vertex_descriptor v1 = target(e, graph);
            const Vertex& vertex0 = graph[v0];
            const Vertex& vertex1 = graph[v1];

            dot << id(vertex0.av0) << "->" << id(vertex1.av0) << "[label=\"" <<
                id(edge.ae0) << "\\n" << id(edge.ae1) <<
                "\"]\n";
        }
        dot << "}\n";
    }
}
#endif


#if 0
void AssemblyGraph::findStrandContacts() const
{
    const AssemblyGraph& assemblyGraph = *this;

    while(true) {


        cout << "Vertices with in-degree 0 are:";
        BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
            if(in_degree(v, assemblyGraph) == 0) {
                cout << " " << id(v);
            }
        }
        cout << endl;

        cout << "Enter a vertex id and maximum distance:" << endl;
        uint64_t vertexId;
        uint64_t maxDistance;
        cin >> vertexId >> maxDistance;
        if(not cin) {
            break;
        }

        vertex_descriptor v = null_vertex();
        BGL_FORALL_VERTICES(u, assemblyGraph, AssemblyGraph) {
            if(id(u) == vertexId) {
                v = u;
                break;
            }
        }
        if(v == null_vertex()) {
            cout << "A vertex with this id does not exist." << endl;
        }

        cout << "v is " << id(v) << endl;

    
        // x = vertices reachable from v, moving forward.
        vector<vertex_descriptor> x;
        findDescendants(assemblyGraph, v, maxDistance, x);
        std::ranges::sort(x, orderById);

        // y = the reverse complements of x.
        // Also, vertices reachable from the reverse complement of v, moving backward.
        vector<vertex_descriptor> y;
        for(const vertex_descriptor v: x) {
            y.push_back(assemblyGraph[v].vRc);
        }
        std::ranges::sort(y, orderById);

        // z = intersection of x and y.
        vector<vertex_descriptor> z;
        std::ranges::set_intersection(x, y, back_inserter(z), orderById);

        ofstream csv("findStrandContacts.csv");
        csv << "Segment,Color,\n";
        uint64_t segmentCount = 0;
        BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
            const vertex_descriptor v0 = source(segment, assemblyGraph);
            const vertex_descriptor v1 = target(segment, assemblyGraph);
            string color;
            if(
                std::ranges::binary_search(z, v0, orderById) and 
                std::ranges::binary_search(z, v1, orderById)) {
                color = "red";
                ++segmentCount;
            } else {
                color = "blue";
            }
            csv << id(segment) << "," << color << ",\n";
        }

        cout << "Found " << segmentCount << " segments." << endl;

    }
}
#endif


#if 0
void AssemblyGraph::findStrandContactsInteractive()
{
    AssemblyGraph& assemblyGraph = *this;
    vector<vertex_descriptor> strandContactVertices;


    while(true) {

        cout << 
            "Enter command:\n"
            "   f: Find a strand contact.\n"
            "   p: Process the current strand contact.\n"
            "   c: Compress the AssemblyGraph.\n"
            "   l: Reload AssemblyGraph.\n"
            "   s: Save and write AssemblyGraph.\n"
            "   q: Quit.\n";

        string s;
        cin >> s;
        if(not cin) {
            break;
        }
        if(s.size() != 1) {
            cout << "Enter a single character." << endl;
        }
        const char c = s[0];



        switch(c) 
        {
        default:
            break;


            
        case 'l':
        {
            cout << "Assembly stage to load?" << endl;
            string stage;
            cin >> stage;
            try {
                load(stage);
            } catch (std::exception& e) {
                cout << e.what() << endl;
            }
        }
        break;


            
        case 's':
        {
            cout << "Assembly stage for save?" << endl;
            string stage;
            cin >> stage;
            try {
                write(stage);
                save(stage);
            } catch (std::exception& e) {
                cout << e.what() << endl;
            }

        }
        break;


            
        case 'c':
        {
            cout << "Before compress, the AssemblyGraph has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " segments." << endl;
            strandSymmetricCompress();
            cout << "After compress, the AssemblyGraph has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " segments." << endl;
        }
        break;


            
        case 'f':
        {
            cout << "Enter start segment id and maximum distance:" << endl;
            uint64_t startSegmentId;
            uint64_t maxDistance;
            cin >> startSegmentId >> maxDistance;

            // Find the Segment with this id.
            Segment startSegment = assemblyGraphNullEdge;
            BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
                if(id(segment) == startSegmentId) {
                    startSegment = segment;
                    break;
                }
            }
            if(startSegment == assemblyGraphNullEdge) {
                cout << "A segment with this id does not exist." << endl;
            }

            // Our start vertex is the target vertex of this segment.
            const vertex_descriptor v = target(startSegment, assemblyGraph);
            // cout << "Start vertex is " << id(v) << endl;

            // x = vertices reachable from v, moving forward.
            vector<vertex_descriptor> x;
            findDescendants(assemblyGraph, v, maxDistance, x);
            x.push_back(v);
            std::ranges::sort(x, orderById);

            /*
            cout << "Descendants:" << endl;
            for(const vertex_descriptor v: x) {
                cout << id(v) << " ";
            }
            cout << endl;
            */

            // y = the reverse complements of x.
            // Also, vertices reachable from the reverse complement of v, moving backward.
            vector<vertex_descriptor> y;
            for(const vertex_descriptor v: x) {
                y.push_back(assemblyGraph[v].vRc);
            }
            std::ranges::sort(y, orderById);

            // The strand contact is the intersection of x and y.
            strandContactVertices.clear();
            std::ranges::set_intersection(x, y, back_inserter(strandContactVertices), orderById);

            // Write a file than can be loaded in Bandage to see the strand contact.
            ofstream csv("findStrandContacts.csv");
            csv << "Segment,Steps,Color,\n";
            uint64_t segmentCount = 0;
            BGL_FORALL_EDGES(segment, assemblyGraph, AssemblyGraph) {
                const vertex_descriptor v0 = source(segment, assemblyGraph);
                const vertex_descriptor v1 = target(segment, assemblyGraph);
                string color;
                if(
                    std::ranges::binary_search(strandContactVertices, v0, orderById) and 
                    std::ranges::binary_search(strandContactVertices, v1, orderById)) {
                    color = "LightCoral";
                    ++segmentCount;
                } else {
                    color = "Grey";
                }
                csv << id(segment) << "," << assemblyGraph[segment].size() << "," << color << ",\n";
            }

            cout << "The strand contact has " << 
                strandContactVertices.size() << " vertices and " <<
                segmentCount << " segments." << endl;
            cout << "To see it, load findStrandContacts.csv in Bandage." << endl;
        }
        break;


            
        case 'p':
        {
            cout << "Before processing this strand contact, the AssemblyGraph has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " segments." << endl;

            const StrandSeparation1::StrandContact strandContact(
                assemblyGraph,
                strandContactVertices,
                "Python",
                0);
            removeZeroLengthSegmentsStrandSymmetric();
            removeIsolatedVertices();
            
            cout << "After processing this strand contact, the AssemblyGraph has " <<
                num_vertices(assemblyGraph) << " vertices and " <<
                num_edges(assemblyGraph) << " segments." << endl;

        }
        break;


            
        case 'q':
            return;
        }
    }
}
#endif



// This passes all the vertices ot the strand contact.
void AssemblyGraph::findStrandContactsInteractive()
{
    AssemblyGraph& assemblyGraph = *this;

    vector<vertex_descriptor> strandContactVertices;
    BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
        strandContactVertices.push_back(v);
    }



    cout << "The initial AssemblyGraph has " <<
        num_vertices(assemblyGraph) << " vertices and " <<
        num_edges(assemblyGraph) << " segments." << endl;

    const StrandSeparation1::StrandContact strandContact(
        assemblyGraph,
        strandContactVertices,
        "Python",
        0);

    removeZeroLengthSegmentsStrandSymmetric();
    removeIsolatedVertices();
    strandSymmetricCompress();
    
    cout << "The final AssemblyGraph has " <<
        num_vertices(assemblyGraph) << " vertices and " <<
        num_edges(assemblyGraph) << " segments." << endl;

    write("Z");
    save("Z");
}

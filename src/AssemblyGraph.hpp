#pragma once

// In the AssemblyGraph, each edge is a linear chains of adjacent AnchorPairs.
// It is initially created from linear chains of edges in the AnchorGraph.
// Each AssemblyGraphEdge generates a gfa segment.

// Shasta.
#include "AssemblyGraphBaseClass.hpp"
#include "AnchorPair.hpp"
#include "Base.hpp"
#include "invalid.hpp"
#include "MappedMemoryOwner.hpp"
#include "MultithreadedObject.hpp"
#include "SuperbubbleChain.hpp"

// Boost libraries.
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

// Standard library.
#include <list>
#include "memory.hpp"
#include "tuple.hpp"


namespace shasta2 {


    class AssemblyGraphEdgeStep;
    const AssemblyGraphBaseClass::edge_descriptor assemblyGraphNullEdge = {0, 0, 0};

    class AnchorGraph;
    class Anchors;
    class Assembler;
    class Options;
    class Detangler;
    class Superbubble;
    class SuperbubbleChain;

}



// The AnchorId of a vertex is the last AnchorId of the last AnchorPair
// of all incoming edges and the first AnchorId of the first AnchorPair
// of each outgoing edge.
// When the AssemblyGraph is initially created from the AnchorGraph,
// there can be at most one vertex for each AnchorId.
// However that is no longer true after detangling.
class shasta2::AssemblyGraphVertex {
public:
    AnchorId anchorId = invalid<AnchorId>;
    uint64_t id = invalid<uint64_t>;

    // The reverse complement of this vertex.
    AssemblyGraphBaseClass::vertex_descriptor vRc = AssemblyGraphBaseClass::null_vertex();

    AssemblyGraphVertex(AnchorId anchorId, uint64_t id) :
        anchorId(anchorId), id(id) {}
    AssemblyGraphVertex() {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorId;
        ar & id;

        // We can't serialize vRc.
        // See AssemblyGraph::serialize, which serializes
        // reverse complement ids instead of vertex descriptors.
    }
};



class shasta2::AssemblyGraphEdgeStep {
public:
    AnchorPair anchorPair;
    uint64_t offset = invalid<uint64_t>;
    vector<Base> sequence;

    AssemblyGraphEdgeStep(
        const AnchorPair& anchorPair,
        uint64_t offset) :
        anchorPair(anchorPair),
        offset(offset)
    {}

    AssemblyGraphEdgeStep()
    {}

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & anchorPair;
        ar & offset;
        ar & sequence;

        // We can't serialize eRc.
        // See AssemblyGraph::serialize, which serializes
        // reverse complement ids instead of vertex descriptors.
    }
};



class shasta2::AssemblyGraphEdge : public vector<AssemblyGraphEdgeStep> {
public:
    uint64_t id = invalid<uint64_t>;
    bool wasAssembled = false;

    // The reverse complement of this edge.
    AssemblyGraphBaseClass::edge_descriptor eRc = assemblyGraphNullEdge;

    // The reverse complement of this edge.
    // AssemblyGraphBaseClass::edge_descriptor eRc;

    AssemblyGraphEdge(uint64_t id = invalid<uint64_t>) : id(id) {}

    void check(const Anchors&) const;

    uint64_t offset() const;
    uint64_t sequenceLength() const;
    void getSequence(vector<Base>&) const;

    // If wasAssembled is set, this returns sequenceLength().
    // Otherwise, it returns offset().
    uint64_t length() const;

    double averageCoverage() const;
    double lengthWeightedAverageCoverage() const;

    // The annotation is only used for testing/debugging.
    // It can be set via AssemblyGraph::setAnnotation, which is exposed to Python.
    string annotation;

    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        ar & boost::serialization::base_object< vector<AssemblyGraphEdgeStep> >(*this);
        ar & id;
        ar & wasAssembled;
        ar & annotation;

        // We can't serialize eRc.
        // See AssemblyGraph::serialize, which serializes
        // reverse complement ids instead of vertex descriptors.
    }

    void swapSteps(AssemblyGraphEdge& that)
    {
        vector<AssemblyGraphEdgeStep>::swap(that);
    }

    AnchorId firstAnchorId() const
    {
        return front().anchorPair.anchorIdA;
    }
    AnchorId lastAnchorId() const
    {
        return back().anchorPair.anchorIdB;
    }
};



class shasta2::AssemblyGraph :
    public AssemblyGraphBaseClass,
    public MappedMemoryOwner,
    public MultithreadedObject<AssemblyGraph> {
public:

    // Initial construction from the AnchorGraph.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const AnchorGraph&,
        const Options&);

    // Deserialize constructor.
    AssemblyGraph(
        const Anchors&,
        const Journeys&,
        const Options&,
        const string& stage);

    // Detangle, phase, assemble sequence, output.
    void simplifyAndAssemble();

    void setAnnotation(edge_descriptor, const string&);

    const Anchors& anchors;
    const Journeys& journeys;
    uint64_t nextVertexId = 0;
    uint64_t nextEdgeId = 0;
    const Options& options;

    void check(bool writeDetails = false) const;


    // Clear reverse complement information from all vertices and edges.
    // This needs to be done before operations that don't maintain
    // the vertices/edges vRc/eRc field.
    void clearReverseComplementInformation();

    // This creates a vertex vB, identical to the reverse complement
    // of vertex vA. It sets the vRc fields in vA and vB
    // and returns vB.
    vertex_descriptor createReverseComplementVertex(vertex_descriptor vA);

    // This creates an edge eB, identical to the reverse complement
    // of edge eA. It sets the eRc fields in eA and eB
    // and returns eB.
    // The vertices of eB must already exist.
    edge_descriptor createReverseComplementEdge(edge_descriptor eA);

private:


    // BUBBLE CLEANUP.
    // A bubble is a set of parallel edges in the AssemblyGraph.
    class Bubble {
    public:
        vertex_descriptor v0;
        vertex_descriptor v1;
        vector<edge_descriptor> edges;
    };

    // Find Bubbles.
    // The edges of each Bubble are sorted by id.
    void findBubbles(vector<Bubble>&) const;

public:
    uint64_t bubbleCleanup();
    uint64_t bubbleCleanupIterationMultithreaded(
        vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList,
        uint64_t threadCount);
    void bubbleCleanupIterationThreadFunction(uint64_t threadId);
    class BubbleCleanupIterationData {
    public:
        vector<Bubble> candidateBubbles;
        uint64_t modifiedCount = 0;
    };
    BubbleCleanupIterationData bubbleCleanupIterationData;
private:
    bool bubbleCleanup(const Bubble&);


    // Analyze a Bubble and finds pairs of "similar" branches.
    // See analyzeSimilarSequences.h for more information.
    bool analyzeBubble(
        const Bubble&,
        const vector<uint64_t> minRepeatCount,
        vector< pair<uint64_t, uint64_t> >& similarPairs
        ) const;

public:



    // BUBBLE PAIRS CLEANUP
    // The bubble cleanup functions above work on one bubble at a time
    // and don't guarantee strand symmetry.
    // The function below work on a pair of reverse complemented bubbles.
    // They require the AssemblyGraph to be strand-symmetric on input,
    // with the vertices vRc fields and edges eRc fields describing
    // the strand symmetry. On output, the AssemblyGraph is again in a
    // similar strand-symmetric state, as described by the vRc amd eRc field.

    // Find pairs of reverse complemented bubbles.
    // The edges of the first bubble in each pair are sorted by id.
    // The edges of the second bubble in each pair are sorted
    // consistently with the ones in the first pair,
    // that is, the reverse complement of p.first.edges[i] is p.second.edges[i].
    using BubblePair = pair<Bubble, Bubble>;
    void findBubblePairs(vector<BubblePair>&) const;

    uint64_t bubblePairCleanup();
    uint64_t bubblePairCleanupIterationMultithreaded(
        vector< pair<vertex_descriptor, vertex_descriptor> >& excludeList);
    void bubblePairCleanupIterationThreadFunction(uint64_t threadId);
    class BubblePairCleanupIterationData {
    public:
        vector<BubblePair> candidateBubblePairs;
        uint64_t modifiedCount = 0;
    };
    BubblePairCleanupIterationData bubblePairCleanupIterationData;
private:
    bool bubblePairCleanup(const BubblePair&);
public:


    // Compress linear chains of edges into a single edge.
    uint64_t compress();
    edge_descriptor compressLinearChain(const std::list<edge_descriptor>&);
    uint64_t compressDebugLevel = 0; // 1=minimal, 2=compact, 3=detailed.
    uint64_t strandSymmetricCompress();

    // Remove zero length segments by collapsing their vertices.
    void removeZeroLengthSegments();
    void collapseVertices(const vector<vertex_descriptor>&);

    // This cleans up linear chains by removing edges that have low
    // corrected Jaccard similarity with nearby edges and
    // replacing with new edges, constructed by connecting
    // the remaining edges.
    void cleanupLinearChains();
    void cleanupLinearChain(const vector<edge_descriptor>&);



    // Class to order vertices or edges by id.
    class OrderById {
    public:
        OrderById(const AssemblyGraph& assemblyGraph): assemblyGraph(assemblyGraph) {}
        const AssemblyGraph& assemblyGraph;
        bool operator()(vertex_descriptor x, vertex_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }
        bool operator()(edge_descriptor x, edge_descriptor y) const
        {
            return assemblyGraph[x].id < assemblyGraph[y].id;
        }

        // Also order pairs of edges.
        using EdgePair = pair<edge_descriptor, edge_descriptor>;
        bool operator()(const EdgePair& x,const EdgePair& y) const
        {
            if(assemblyGraph[x.first].id < assemblyGraph[y.first].id) {
                return true;
            }
            if(assemblyGraph[x.first].id > assemblyGraph[y.first].id) {
                return false;
            }
            return assemblyGraph[x.second].id < assemblyGraph[y.second].id;
        }
    };
    const OrderById orderById;



    bool hasSelfEdge(vertex_descriptor v) const
    {
        bool edgeExists = false;
        tie(ignore, edgeExists) = boost::edge(v, v, *this);
        return edgeExists;
    }



    // Return true if the two specified edges can be connected for assembly.
    // If necessary, this construct the RestrictedAnchorGraph between the two edges
    // and checks that all is good.
    bool canConnect(edge_descriptor, edge_descriptor) const;

    // The detangling process can generate empty edges (edges without steps).
    // This removes them by collapsing the vertices they join.
    void removeEmptyEdges();

    // Remove isolated vertices.
    void removeIsolatedVertices();

    // Remove connected components with a low N50.
    void removeLowN50Components();


    void connectDanglingSegments();


public:

    // Compute compressed journeys in the AssemblyGraph.
    // The compressed journey of an oriented read
    // is the sequence of assembly graph edges it visits.
    void computeJourneys();
    vector< vector<edge_descriptor> > compressedJourneys;



    // Superbubbles and SuperbubbleChains.
    // - A Superbubble is a Tangle in which all entrance edges are into
    //   a single vertex and all exit edges are from a single vertex.
    // - A SuperbubbleChain is a linear sequence of Superbubbles.
    //   In the sequence, the target vertex of a Superbubble is the same
    //   as the source vertex of the next Superbubble in the SuperbubbleChain.

    // This finds all Superbubbles seen using options.findSuperbubblesMaxDistance.
    // Some pairs of Superbubble can intersect (that is, they can have common edges).
    void findSuperbubbles(vector<Superbubble>&) const;

    // Remove Superbubbles that are entirely contained in a larger superbubble.
    void removeContainedSuperbubbles(vector<Superbubble>&) const;

    void writeSuperbubbles(const vector<Superbubble>&, const string& fileName) const;
    void writeSuperbubblesForBandage(const vector<Superbubble>&, const string& fileName) const;

    void findSuperbubbleChains(
        const vector<Superbubble>&,
        vector<SuperbubbleChain>&
        ) const;

    void writeSuperbubbleChains(const vector<SuperbubbleChain>&, const string& fileName) const;
    void writeSuperbubbleChainsForBandage(const vector<SuperbubbleChain>&, const string& fileName) const;



    // Phasing of SuperbubbleChains.
    void phaseSuperbubbleChains();
    void phaseSuperbubbleChainsThreadFunction(uint64_t threadId);

    // Strand-symmetric phasing of superbubble chains.
    void strandSymmetricPhaseSuperbubbleChains();
    void strandSymmetricPhaseSuperbubbleChainsThreadFunction(uint64_t threadId);

    // This phases the first SuperbubbleChain of a reverse complemented pair,
    // then copies it with reverse complement to replace the second SuperbubbleChain.
    void strandSymmetricPhase(
        const SuperbubbleChain&, uint64_t,
        const SuperbubbleChain&, uint64_t);



    class PhaseSuperbubbleChainsData {
    public:
        vector<SuperbubbleChain> superbubbleChains;

        // Pairs of reverse complemented SuperbubbleChains.
        // These are indexes in the superbubbleChains vector.
        vector< pair<uint64_t, uint64_t> > superbubbleChainPairs;
    };
    PhaseSuperbubbleChainsData phaseSuperbubbleChainsData;



    // Strongly connected components.

    // Find the non-trivial strongly connected components.
    // Each component is stored with vertices sorted to permit binary searches.
    void findStrongComponents(vector< vector<vertex_descriptor> >&) const;

    // This creates a csv file that can be loaded in bandage to see
    // the strongly connected components.
    void colorStrongComponents() const;

    // Pruning.
    void prune();
    uint64_t pruneIteration();



    // Read following.
    // Note assemlbyPaths are not necessarily paths in the AssemblyGraph.
    // There may be jumps, which are bridged using local assemblies.
    // Each assembly path generates a new linear chain of edges,
    // which is left in an uncompressed state.
    void readFollowing();


    // Simple connection of two segments (edges) without using
    // the RestrictedAnchorGraph.
    void simpleConnect(edge_descriptor, edge_descriptor);
    bool canSimpleConnect(edge_descriptor, edge_descriptor);
    void findOrientedReadIdsForSimpleConnect(edge_descriptor, edge_descriptor, vector<OrientedReadId>&) const;


    // Output.
    void write(const string& stage);
    void writeIntermediateStageIfRequested(const string& name);
    void writeFasta(const string& stage) const;
private:
    void writeGfa(const string& fileName) const;
    void writeGfa(ostream&) const;
    void writeGraphviz(const string& fileName) const;
    void writeGraphviz(ostream&) const;
    void writeCsv(const string& fileName) const;
    void writeCsv(ostream&) const;
    void writeDetailsCsv(const string& fileName) const;
    void writeDetailsCsv(ostream&) const;
    void writeSequenceLengthByCoverageCsv(const string& fileName) const;
    void writeSequenceLengthByCoverageCsv(ostream&) const;



    // Sequence assembly.

    // Assemble sequence for all edges.
public:
    void assembleAll();

    // Clear sequence from all steps of all edges.
    void clearSequence();
private:

    // Assemble sequence for the specified edge.
    void assemble(edge_descriptor);

    // Assemble sequence for step i of the specified edge.
    // This is the lowest level sequence assembly function and is not multithreaded.
    // It runs a LocalAssembly2 on the AnchorPair for that step.
    void assembleStep(edge_descriptor, uint64_t i);

    // Assemble sequence for all edges in the edgesToBeAssembled vector.
    // This fills in the stepsToBeAssembled with all steps of those edges,
    // then assembles each of the steps in parallel.
    void assemble();
    void assembleThreadFunction(uint64_t threadId);
    vector<edge_descriptor> edgesToBeAssembled;
    vector< pair<edge_descriptor, uint64_t> > stepsToBeAssembled;




    // Serialization.
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive& ar, unsigned int /* version */)
    {
        AssemblyGraph& assemblyGraph = *this;

        ar & boost::serialization::base_object<AssemblyGraphBaseClass>(*this);
        ar & nextVertexId;
        ar & nextEdgeId;

        // Serialize reverse complement vertices and edges.
        // Use ids instead of vertex_descriptors and edge descriptors
        // because vertex_descriptors and edge_descriptors cannot be serialized.
        if(Archive::is_saving::value) {

            // Map the id of a vertex to the id of its reverse complement,
            // then save this map.
            std::map<uint64_t, uint64_t> rcVertexMap;
            BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
                const AssemblyGraphVertex& vertex = assemblyGraph[v];
                const vertex_descriptor vRc = vertex.vRc;
                if(vRc != null_vertex()) {
                    rcVertexMap.insert(make_pair(vertex.id, assemblyGraph[vRc].id));
                }
            }
            ar & rcVertexMap;

            // Map the id of an edge to the id of its reverse complement,
            // then save this map.
            std::map<uint64_t, uint64_t> rcEdgeMap;
            BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
                const AssemblyGraphEdge& edge = assemblyGraph[e];
                const edge_descriptor eRc = edge.eRc;
                if(eRc != assemblyGraphNullEdge) {
                    rcEdgeMap.insert(make_pair(edge.id, assemblyGraph[eRc].id));
                }
            }
            ar & rcEdgeMap;

        } else {
            std::map<uint64_t, uint64_t> rcVertexMap;
            ar & rcVertexMap;
            std::map<uint64_t, vertex_descriptor> vertexMap;
            BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
                vertexMap.insert(make_pair(assemblyGraph[v].id, v));
            }
            BGL_FORALL_VERTICES(v, assemblyGraph, AssemblyGraph) {
                AssemblyGraphVertex& vertex = assemblyGraph[v];
                const auto it = rcVertexMap.find(vertex.id);
                if(it != rcVertexMap.end()) {
                    const uint64_t idRc = it->second;
                    vertex.vRc = vertexMap.at(idRc);
                }
            }


            std::map<uint64_t, uint64_t> rcEdgeMap;
            ar & rcEdgeMap;
            std::map<uint64_t, edge_descriptor> edgeMap;
            BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
                edgeMap.insert(make_pair(assemblyGraph[e].id, e));
            }
            BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
                AssemblyGraphEdge& edge = assemblyGraph[e];
                const auto it = rcEdgeMap.find(edge.id);
                if(it != rcEdgeMap.end()) {
                    const uint64_t idRc = it->second;
                    edge.eRc = edgeMap.at(idRc);
                }
            }
        }
    }



    void save(ostream&) const;
    void load(istream&);

    // These do save/load to/from mapped memory.
    // The file name is AssemblyGraph-Stage.
    void save(const string& stage) const;
    void load(const string& stage);
};

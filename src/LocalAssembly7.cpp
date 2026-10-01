// Shasta.
#include "LocalAssembly7.hpp"
#include "abpoaWrapper.hpp"
#include "Anchor.hpp"
#include "deduplicate.hpp"
#include "findReachableVertices.hpp"
#include "graphvizToHtml.hpp"
#include "msaRepair.hpp"
#include "Reads.hpp"
#include "theseusWrapper.hpp"
#include "tmpDirectory.hpp"
using namespace shasta2;

// Boost libraries.
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "chrono.hpp"
#include <cmath>
#include "fstream.hpp"
#include <iomanip>
#include "iostream.hpp"
#include <map>
#include <stack>



LocalAssembly7::LocalAssembly7(
    const Options& options,
    const Anchors& anchors,
    shared_ptr<const HomopolymerModel> homopolymerModelPointer,
    AnchorId anchorIdA,
    AnchorId anchorIdB,
    ostream& html,
    const vector<OrientedReadId>& orientedReadIds) :
    options(options),
    anchors(anchors),
    homopolymerModelPointer(homopolymerModelPointer),
    anchorIdA(anchorIdA),
    anchorIdB(anchorIdB),
    html(html)
{
    if(html) {
        html <<
            "<table>" <<
            "<tr><th class=left>Left anchor (anchor A)<td class=centered>" << anchorIdToString(anchorIdA) <<
            "<tr><th class=left>Right anchor (anchor B)<td class=centered>" << anchorIdToString(anchorIdB) <<
            "</table>";
    }

    // Fill in the orientedReads and the sequenceInfos.
    gatherOrientedReads(orientedReadIds);
    removeOutliers();
    estimateOffset();
    gatherSequences();
    writeOrientedReads();
    writeSequences();

    // Run the local assembly.
    run();

    writeSequence();
}



void LocalAssembly7::run()
{
    // Find total coverage constrained on both anchors.
    const uint64_t commonCoverage = getTotalCommonCoverage();


    // Try the fast path, if allowed and possible.
    if(commonCoverage >= options.commonCoverageThreshold) {
        runFastPath();
        if(success) {
            return;
        }
    }


    // If the fast path was not allowed or not possible,
    // use the requested Method.
    switch(options.method) {
    case Method::Adaptive:
        runAdaptive();
        break;
    case Method::Abpoa:
        runAbpoa();
        break;
    case Method::TheseusOnly:
        runTheseus(false);
        break;
    case Method::TheseusAll:
        runTheseus(true);
        break;
    case Method::DeBruijn:
        runDeBruijn();
        break;
    default:
        throw runtime_error("Invalid LocalAssembly7 Method.");
    }
}



// This returns the total coverage in sequences that are on both anchors.
uint64_t LocalAssembly7::getTotalCommonCoverage() const
{
    uint64_t commonCoverage = 0;
    for(const SequenceInfo& sequenceInfo: sequences) {
        if(sequenceInfo.isOnAnchorA and sequenceInfo.isOnAnchorB) {
            commonCoverage += sequenceInfo.coverage();
        }
    }
    return commonCoverage;
}



// This returns the maximul length of sequences that are on both anchors.
uint64_t LocalAssembly7::getMaxLengthCommon() const
{
    uint64_t maxLength = 0;
    for(const SequenceInfo& sequenceInfo: sequences) {
        if(sequenceInfo.isOnAnchorA and sequenceInfo.isOnAnchorB) {
            maxLength = max(maxLength, sequenceInfo.sequence.size());
        }
    }
    return maxLength;

}



void LocalAssembly7::runFastPath()
{
    if(not options.allowFastPath) {
        return;
    }

    // Can only use the fast path if the number of oriented reads on both anchors
    // is at least equal to commonThreshold.
    uint64_t commonCoverage = getTotalCommonCoverage();
    if(commonCoverage < options.commonCoverageThreshold) {
        return;
    }

    // Can only use fast path if, of those commonCoverage oriented reads,
    // at least a fraction equal to fastPathFractionThreshold
    // have the same sequence.
    // So we find the sequence with the most coverage
    // among the ones constrained on both sides.
    uint64_t sequenceIdBest = invalid<uint64_t>;
    uint64_t coverageBest = 0;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        if(sequenceInfo.isOnAnchorA and sequenceInfo.isOnAnchorB) {
            const uint64_t coverage = sequenceInfo.coverage();
            if(coverage > coverageBest) {
                coverageBest = coverage;
                sequenceIdBest = sequenceId;
            }
        }
    }
    const double fastPathFraction = double(coverageBest) / double(commonCoverage);
    if(fastPathFraction < options.fastPathFractionThreshold) {
        return;
    }

    // All is good. We can just store that best sequence as our consensus.
    sequence = sequences[sequenceIdBest].sequence;
    success = true;

    if(html) {
        html <<
            "<h3>Fast path </h3>"
            "<table>"
            "<tr><td class=left>Number of oriented reads on both anchors<td class=centered>" <<
            commonCoverage <<
            "<tr><td class=left>Sequence id of most frequent sequence of oriented "
            "reads on both anchors<td class=centered>" <<
            sequenceIdBest <<
            "<tr><td class=left>Number of oriented reads with that most frequent sequence<td class=centered>" <<
            coverageBest <<
            "<tr><td class=left>Fraction of oriented reads on both anchors "
            "that have that most frequent sequence<td class=centered>" <<
            std::setprecision(2) << std::fixed << fastPathFraction <<
            "</table>";
    }

}



// Construct the Dr Bruijn graph
// and, if successful, use it to assemble sequence.
void LocalAssembly7::runDeBruijn()
{
    for(SequenceInfo& sequenceInfo: sequences) {
        sequenceInfo.constructDeBruijnSequence(options.k);
        sequenceInfo.constructKmers(options.k);
    }
    gatherKmers();

    Graph graph;
    createGraph(options.k, graph);
    if(html) {
        html <<
            "<br>After removing unreachable vertices, the De Bruijn graph has " <<
            num_vertices(graph) <<
            " vertices and " << num_edges(graph) << " edges.";
    }


    // If the graph has cycles, give up, leaving success set to false.
    {
        // Create a vertexIndexMap, needed below.
        std::map<vertex_descriptor, uint64_t> vertexIndexMap;
        uint64_t vertexIndex = 0;
        BGL_FORALL_VERTICES(v, graph, Graph) {
            vertexIndexMap.insert(make_pair(v, vertexIndex++));
        }

        std::map<vertex_descriptor, uint64_t> componentMap;
        if(boost::strong_components(
            graph,
            boost::make_assoc_property_map(componentMap),
            vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)))
            != num_vertices(graph)) {
            if(html) {
                html << "<br>The De Bruijn graph contains cycles.";
                graph.writeVertices(kmers, "DeBruijnGraph.csv");
                writeKmerOccurrences(graph, "DeBruijnGraph-KmerOccurrences.csv");
                writeGraph(options.k, graph);
            }
            return;
        }
    }

    // Compute the optimal assembly path.
    graph.computeAssemblyPath();

    if(html) {
        graph.writeVertices(kmers, "DeBruijnGraph.csv");
        writeKmerOccurrences(graph, "DeBruijnGraph-KmerOccurrences.csv");
        writeAssemblyPath(graph);
    }
    writeGraph(options.k, graph);

    // Assemble sequence.
    assemble(options.k, graph);
    success = true;
}



void LocalAssembly7::gatherOrientedReads(
    const vector<OrientedReadId>& orientedReadIds)
{
    const Anchor anchorA = anchors[anchorIdA];
    const Anchor anchorB = anchors[anchorIdB];

    if(html) {
        html << "<h3>" << orientedReadIds.size() << " input oriented reads</h3><table>";
        for(const OrientedReadId orientedReadId: orientedReadIds) {
            html << "<tr><td class=centered>" << orientedReadId;
        }
        html << "</table>";
    }



    // Joint loop over the input OrientedReadIds and
    // the OrientedReadIds of the two Anchors.
    auto itA = anchorA.begin();
    const auto endA = anchorA.end();
    auto itB = anchorB.begin();
    const auto endB = anchorB.end();
    for(const OrientedReadId orientedReadId: orientedReadIds) {

        // Check if this OrientedReadId is on the two anchors.
        while((itA != endA) and (itA->orientedReadId < orientedReadId)) {
            ++itA;
        }
        while((itB != endB) and (itB->orientedReadId < orientedReadId)) {
            ++itB;
        }

        const bool isOnA = (itA != endA) and (itA->orientedReadId == orientedReadId);
        const bool isOnB = (itB != endB) and (itB->orientedReadId == orientedReadId);

        // If on neither anchor, this OrientedReadId cannot be used.
        if(not (isOnA or isOnB)) {
            continue;
        }

        // If on both anchors and negative offset, this OrientedReadId cannot be used.
        if(isOnA and isOnB and (itA->position > itB->position)) {
            continue;
        }

        // We can use this OrientedReadId for assembly.

        // Create an OrientedRead.
        OrientedRead& orientedRead = orientedReads.emplace_back();
        orientedRead.orientedReadId = orientedReadId;
        if(isOnA) {
            orientedRead.positionA = itA->position;
        }
        if(isOnB) {
            orientedRead.positionB = itB->position;
        }

    }
}



void LocalAssembly7::estimateOffset()
{
    uint64_t sum = 0;
    uint64_t n = 0;
    for(const OrientedRead& orientedRead: orientedReads) {
        if(orientedRead.isOnBothAnchors()) {
            sum += orientedRead.positionOffsetAB();
            ++n;
        }
    }

    SHASTA2_ASSERT(n > 0);
    offset = uint32_t(std::round(double(sum) / double(n)));
}



void LocalAssembly7::gatherSequences()
{
    // For reads fixed on one side only, we use a sequence length
    // equal to offset + aDrift * offset + bDrift.
    const uint32_t length = offset + uint32_t(std::round(
        options.aExtend * double(offset) + options.bExtend
        ));


    // Fill in the beginPosition and endPosition of each OrientedRead.
    // Those are the position ranges of the sequences that will be used for assembly.
    for(OrientedRead& orientedRead: orientedReads) {

        if(orientedRead.isOnBothAnchors()) {
            orientedRead.positionBegin = orientedRead.positionA;
            orientedRead.positionEnd  = orientedRead.positionB;
        } else if(orientedRead.isOnAnchorA()) {
            SHASTA2_ASSERT(not orientedRead.isOnAnchorB());
            const ReadId readId = orientedRead.orientedReadId.getReadId();
            const uint32_t readLength = uint32_t(anchors.reads.getReadSequenceLength(readId));
            orientedRead.positionBegin = orientedRead.positionA;
            orientedRead.positionEnd  = min(readLength, orientedRead.positionBegin + length);
        } else if(orientedRead.isOnAnchorB()) {
            SHASTA2_ASSERT(not orientedRead.isOnAnchorA());
            orientedRead.positionEnd = orientedRead.positionB;
            if(orientedRead.positionEnd > length) {
                orientedRead.positionBegin = orientedRead.positionEnd - length;
            } else {
                orientedRead.positionBegin  = 0;
            }
        } else {
            SHASTA2_ASSERT(0);
        }
    }



    // Now we can create the sequences.
    const Reads& reads = anchors.reads;
    vector<Base> sequence;
    for(OrientedRead& orientedRead: orientedReads) {
        const OrientedReadId orientedReadId = orientedRead.orientedReadId;

        // Create the sequence of this read (portion to be used for this local assembly).
        sequence.clear();
        for(uint32_t position=orientedRead.positionBegin; position!=orientedRead.positionEnd; position++) {
            sequence.push_back(reads.getOrientedReadBase(orientedReadId, position));
        }

        // See if we already have this sequence in this table.
        bool found = false;
        for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
            SequenceInfo& sequenceInfo = sequences[sequenceId];
            if(
                (orientedRead.isOnAnchorA() == sequenceInfo.isOnAnchorA) and
                (orientedRead.isOnAnchorB() == sequenceInfo.isOnAnchorB) and
                (sequence == sequenceInfo.sequence)) {
                orientedRead.sequenceId = sequenceId;
                sequenceInfo.orientedReadIds.push_back(orientedReadId);
                found = true;
                break;
            }
        }
        if(not found) {
            orientedRead.sequenceId = sequences.size();
            sequences.emplace_back(orientedRead.isOnAnchorA(), orientedRead.isOnAnchorB(), orientedReadId, sequence);
        }
    }
}



void LocalAssembly7::writeOrientedReads() const
{
    if(not html) {
        return;
    }

    html << "<h3>" << orientedReads.size() << " usable oriented reads</h3>";

    html << "<table><tr>"
        "<th>Oriented<br>read id"
        "<th>On A"
        "<th>On B"
        "<th>PositionA"
        "<th>PositionB"
        "<th>Length<br>(bases)"
        "<th>PositionAB<br>offset"
        "<th>Assembly<br>position<br>begin"
        "<th>Assembly<br>position<br>end"
        "<th>Assembly<br>position<br>length"
        "<th>Sequence<br>id"
        ;


    uint64_t commonCount = 0;
    for(const OrientedRead& orientedRead: orientedReads) {
        if(orientedRead.isOnBothAnchors()) {
            ++commonCount;
        }
        html <<
            "<tr>"
            "<th class=centered>" << orientedRead.orientedReadId <<
            "<td class=centered>" << (orientedRead.isOnAnchorA() ? "&check;" : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorB() ? "&check;" : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorA() ? to_string(orientedRead.positionA) : "") <<
            "<td class=centered>" << (orientedRead.isOnAnchorB() ? to_string(orientedRead.positionB) : "") <<
            "<td class=centered>" << anchors.reads.getReadSequenceLength(orientedRead.orientedReadId.getReadId()) <<
            "<td class=centered>" << (orientedRead.isOnBothAnchors() ? to_string(orientedRead.positionOffsetAB()) : "") <<
            "<td class=centered>" << orientedRead.positionBegin <<
            "<td class=centered>" << orientedRead.positionEnd <<
            "<td class=centered>" << orientedRead.positionOffsetForAssembly() <<
            "<td class=centered>" << orientedRead.sequenceId
            ;
    }

    html << "</table>";

    html << "<br>Estimated offset using " << commonCount <<
        " oriented reads common to the left and right anchors is " << offset << " bases.";


}



void LocalAssembly7::writeSequences() const
{
    if(not html) {
        return;
    }

    html << "<h3>Oriented read sequences used for assembly</h3>"
        "<table><tr>"
        "<th>Sequence<br>id"
        "<th>On A"
        "<th>On B"
        "<th>Length"
        "<th>Coverage"
        "<th class=left>Sequence";

    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        html <<
            "<tr>"
            "<td class=centered>" << sequenceId <<
            "<td class=centered>" << (sequenceInfo.isOnAnchorA ? "&check;" : "") <<
            "<td class=centered>" << (sequenceInfo.isOnAnchorB ? "&check;" : "") <<
            "<td class=centered>" << sequenceInfo.sequence.size() <<
            "<td class=centered>" << sequenceInfo.orientedReadIds.size() <<
            "<td class=left style='font-family:monospace;white-space: nowrap'>";
        std::ranges::copy(sequenceInfo.sequence, ostream_iterator<Base>(html));
    }

    html << "</table>";
}



void LocalAssembly7::SequenceInfo::constructDeBruijnSequence(uint64_t k)
{
    deBruijnSequence.clear();

    if(isOnAnchorA) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnSequence.push_back(Base::fromInteger(uint8_t(10)));
        }
    }
    std::ranges::copy(sequence, back_inserter(deBruijnSequence));
    if(isOnAnchorB) {
        for(uint64_t i=0; i<k; i++) {
            deBruijnSequence.push_back(Base::fromInteger(uint8_t(20)));
        }
    }
}



void LocalAssembly7::SequenceInfo::constructKmers(uint64_t k)
{
    kmers.clear();

    // This can be made faster by shifting the k-mers instead of constructing
    // all of them from scratch.
    Kmer kmer;
    kmer.reserve(k);
    for(uint64_t position=0; position+k<=deBruijnSequence.size(); position++) {
        kmer.clear();
        for(uint64_t i=0; i<k; i++) {
            kmer.push_back(deBruijnSequence[position + i]);
        }
        kmers.push_back(kmer);
    }
}



void LocalAssembly7::gatherKmers()
{
    std::map<Kmer, vector<KmerOccurrence> > kmersMap;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        vector<Kmer>& sequenceKmers = sequences[sequenceId].kmers;
        for(uint64_t position=0; position<sequenceKmers.size(); position++) {
            const Kmer& kmer = sequenceKmers[position];
            kmersMap[kmer].push_back(KmerOccurrence(sequenceId, position));
        }
    }

    kmers.clear();
    std::ranges::copy(kmersMap, back_inserter(kmers));

    if(html) {
        html << "<br>There are " << kmers.size() << " De Bruijn k-mers.";
    }
}



void LocalAssembly7::createGraph(uint64_t k, Graph& graph)
{

    // Vectors to store the vertex corresponding to each position of each sequence.
    vector< vector<vertex_descriptor> > vertexTable;
    for(const SequenceInfo& sequenceInfo: sequences) {
        vector<vertex_descriptor>& v = vertexTable.emplace_back();
        v.resize(sequenceInfo.deBruijnSequence.size() - k + 1, Graph::null_vertex());
    }



    // Each distinct k-mer generates a vertex.
    for(uint64_t kmerId=0; kmerId<kmers.size(); kmerId++) {
        const auto& [kmer, occurrences] = kmers[kmerId];

        uint64_t coverage = 0;
        for(const KmerOccurrence& occurrence: occurrences) {
            coverage += sequences[occurrence.sequenceId].coverage();
        }
        const Graph::vertex_descriptor v = boost::add_vertex(
            Vertex(graph.nextVertexId++, kmerId, occurrences, coverage), graph);
        for(const KmerOccurrence& occurrence: occurrences) {
            vertexTable[occurrence.sequenceId][occurrence.position] = v;
        }
    }



    // Flag the A and B vertices.
    // Use the first and last vertex of a sequence fixed on both sides.
    graph.vA = Graph::null_vertex();
    graph.vB = Graph::null_vertex();
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequence = sequences[sequenceId];
        if(sequence.isOnAnchorA and sequence.isOnAnchorB) {
            graph.vA = vertexTable[sequenceId].front();
            graph.vB = vertexTable[sequenceId].back();
        }
    }
    SHASTA2_ASSERT(graph.vA != Graph::null_vertex());
    SHASTA2_ASSERT(graph.vB != Graph::null_vertex());
    graph[graph.vA].isAVertex = true;
    graph[graph.vB].isBVertex = true;



    // Now create the edges.
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const uint64_t coverage = sequences[sequenceId].coverage();
        const vector<vertex_descriptor>& sequenceVertices = vertexTable[sequenceId];
        for(uint64_t position1=1; position1<sequenceVertices.size(); position1++) {
            const uint64_t position0 = position1 - 1;
            const vertex_descriptor v0 = sequenceVertices[position0];
            const vertex_descriptor v1 = sequenceVertices[position1];
            auto [e, edgeExists] = boost::edge(v0, v1, graph);
            if(edgeExists) {
                graph[e].coverage += coverage;
            } else {
                auto[e, ignore] = boost::add_edge(v0, v1, graph);
                graph[e].coverage = coverage;
            }
        }
    }
    if(html) {
        html << "<br>The initial De Bruijn graph has " <<
            num_vertices(graph) << " vertices and " << num_edges(graph) << " edges.";
    }

    graph.removeUnreachableVertices();
    if(html) {
        html << "<br>After removing inaccessible vertices, the De Bruijn graph has " <<
            num_vertices(graph) << " vertices and " << num_edges(graph) << " edges.";
    }

    // Compute edge weights.
    BGL_FORALL_EDGES(e, graph, Graph) {
        Edge& edge = graph[e];
        const double logP = options.logPCoefficient * double(edge.coverage);
        edge.weight = std::pow(10., -0.1 * logP);
    }

}



void LocalAssembly7::writeGraph(uint64_t k, const Graph& graph)
{
    if(not html) {
        return;
    }

    // Write it in graphviz format.
    const string dotFileName = "DeBruijnGraph-" + to_string(k) + ".dot";
    graph.writeGraphviz(dotFileName, kmers);

    // Display it in html in svg format.
    const double timeout = 10.;
    const string options = "-Nshape=rectangle -Gnslimit=1 -Gnslimit1=1 -Gmclimit=1 -Gsearchsize=10 -Gsplines=false";
    html << "<br>";
    try {
        graphvizToHtml(dotFileName, "dot", timeout, options, html, true);
    } catch(std::exception&) {
        html << "Unable to display the graph." << endl;
    }
}



void LocalAssembly7::Graph::writeGraphviz(const string& fileName, const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers) const
{
    ofstream dot(fileName);
    writeGraphviz(dot, kmers);
}



void LocalAssembly7::Graph::writeGraphviz(ostream& dot, const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers) const
{
    const Graph& graph = *this;

    dot << "digraph DeBruijn {\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        dot << vertex.vertexId << " [";
        dot << "label=\"v" << vertex.vertexId << "\\nk" << vertex.kmerId << "\\n" << vertex.coverage;
        const auto&[kmer, ignore] = kmers[vertex.kmerId];
        for(uint64_t i=0; i<kmer.size(); i++) {
            if((i%10) == 0) {
                dot << "\\n";
            }
            dot << kmer[i];
        }
        dot << "\"";
        if(graph[v].isAVertex) {
            dot << " style=filled fillcolor=Pink";
        } else if(graph[v].isBVertex) {
            dot << " style=filled fillcolor=LightBlue";
        } else if(graph[v].isOnAssemblyPath) {
            dot << " style=filled fillcolor=LightGreen";
        }
        dot << "];\n";
    }



    dot << std::fixed << std::setprecision(1);
    BGL_FORALL_EDGES(e, graph, Graph) {
        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);
        const Vertex& vertex0 = graph[v0];
        const Vertex& vertex1 = graph[v1];
        dot << vertex0.vertexId << "->" << vertex1.vertexId <<
            "[penwidth=" << 0.2*double(graph[e].coverage) <<
            " tooltip=\"" << graph[e].coverage << "\"";

        if(vertex0.isOnAssemblyPath and vertex1.isOnAssemblyPath) {
            dot << " color=green";
        }

        dot <<
            "]"
            ";\n";
    }

    dot << "}\n";
}



void LocalAssembly7::Graph::writeVertices(
    const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
    const string& fileName) const
{
    ofstream csv(fileName);
    writeVertices(kmers, csv);
}



void LocalAssembly7::Graph::writeVertices(
    const vector< pair<Kmer, vector<KmerOccurrence> > >& kmers,
    ostream& csv) const
{
    const Graph& graph = *this;
    csv << "VertexId,KmerId,coverage,Kmer,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        if((in_degree(v, graph) == 0) and (out_degree(v, graph) == 0)) {
            continue;
        }
        const Vertex& vertex = graph[v];
        const Kmer& kmer = kmers[vertex.kmerId].first;
        csv << vertex.vertexId << ",";
        csv << vertex.kmerId << ",";
        csv << vertex.coverage << ",";
        std::ranges::copy(kmer, ostream_iterator<Base>(csv));
        csv << ",";
        csv << "\n";
    }
}



void LocalAssembly7::writeAssemblyPath(const Graph& graph) const
{

    ofstream csv("AssemblyPath.csv");
    csv << "Position,VertexId,Coverage,Kmer,\n";

    for(uint64_t position=0; position<graph.assemblyPath.size(); position++) {
        const vertex_descriptor v = graph.assemblyPath[position];
        const Vertex& vertex = graph[v];
        const Kmer& kmer = kmers[vertex.kmerId].first;

        csv << position << ",";
        csv << vertex.vertexId << ",";
        csv << vertex.coverage << ",";
        std::ranges::copy(kmer, ostream_iterator<Base>(csv));
        csv << ",";
        csv << "\n";

    }
}



void LocalAssembly7::Graph::removeUnreachableVertices()
{
    Graph& graph = *this;

    // Remove the vertices that are not reachable from vA moving forward.
    std::set<vertex_descriptor> reachableVertices;
    findReachableVertices(graph, vA, 0, reachableVertices);
    vector<vertex_descriptor> verticesToBeRemoved;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

    // Disconnect the vertices that are not reachable from vB moving backward.
    findReachableVertices(graph, vB, 1, reachableVertices);
    verticesToBeRemoved.clear();
    BGL_FORALL_VERTICES(v, graph, Graph) {
        if(not reachableVertices.contains(v)) {
            verticesToBeRemoved.push_back(v);
        }
    }
    for(const vertex_descriptor v: verticesToBeRemoved) {
        boost::clear_vertex(v, graph);
        boost::remove_vertex(v, graph);
    }

}



void LocalAssembly7::Graph::computeAssemblyPath()
{
    Graph& graph = *this;

    // Create a vertexIndexMap, needed below.
    std::map<vertex_descriptor, uint64_t> vertexIndexMap;
    uint64_t vertexIndex = 0;
    BGL_FORALL_VERTICES(v, graph, Graph) {
        vertexIndexMap.insert(make_pair(v, vertexIndex++));
    }

    // Find a shortest path tree starting at vA.
    std::map<vertex_descriptor, vertex_descriptor> predecessorMap;
    dijkstra_shortest_paths(graph, vA,
       vertex_index_map(boost::make_assoc_property_map(vertexIndexMap)).
       weight_map(boost::get(&Edge::weight, graph)).
       predecessor_map(boost::make_assoc_property_map(predecessorMap))
       );

    // Walk back the predecessorMap staring at vB.
    vertex_descriptor v = vB;
    while(true) {
        assemblyPath.push_back(v);
        if(v == vA) {
            break;
        }
        v = predecessorMap.at(v);
    }
    std::ranges::reverse(assemblyPath);

    for(const vertex_descriptor v: assemblyPath) {
        graph[v].isOnAssemblyPath = true;
    }
}




// This uses the assembly path to assemble sequence.
// The assembly path begins at vA, which contains invalid sequence (Base::fromInteger(10)),
// and ends at vB, which also contains invalid sequence (Base::fromInteger(20)).
// If the assembly path has N vertices, the base offset between vA and vB is N-1.
// The true assembled sequence begins k bases after vA and ends at the base before vB,
// so its length is n = N-1-k.

// For each assembled base position, we can choose among k vertices the one
// from which we get the base, which will be the same for all of the k vertices.
void LocalAssembly7::assemble(uint64_t k, Graph& graph)
{

    const uint64_t N = graph.assemblyPath.size();
    const uint64_t n = N - 1 - k;
    sequence.clear();
    sequence.resize(n);
    if(html) {
        html << "<br>The assembly path contains " << N << " vertices.";
        html << "<br>Assembled sequence is " << n << " bases long.";

        uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
        for(const vertex_descriptor v: graph.assemblyPath){
            minCoverage = min(minCoverage, graph[v].coverage);
        }
        html << "<br>Minimum coverage on the assembly path is " << minCoverage << ".";
    }




    // Look over all positions of assembled sequence.
    for(uint64_t position=0; position<n; position++) {

        // Loop over the k vertices we can use to get the
        // base at this position.
        Base b;
        for(uint64_t i=0; i<k; i++){
            const uint64_t positionInPath = position + k - i;
            const uint64_t positionInVertexKmer = i;
            const vertex_descriptor v = graph.assemblyPath[positionInPath];
            const Vertex& vertex = graph[v];
            const Kmer& kmer = kmers[vertex.kmerId].first;
            if(i == 0) {
                b = kmer[positionInVertexKmer];
            } else {
                SHASTA2_ASSERT(b == kmer[positionInVertexKmer]);
            }
        }
        sequence[position] = b;
    }
}



void LocalAssembly7::writeConsensus(const vector< pair<Base, uint64_t> >& consensus) const
{
    if(not html) {
        return;
    }

    uint64_t minCoverage = std::numeric_limits<uint64_t>::max();
    uint64_t maxCoverage = 0;
    for(const auto& [ignore, coverage] : consensus) {
        minCoverage = min(minCoverage, coverage);
        maxCoverage = max(maxCoverage, coverage);
    }

    html <<
        "<h3>Assembled sequence with coverage</h3>"
        "<table>"
        "<tr><th class=left>Length<td class=left>" << consensus.size();

    // Write a line with the sequence.
    html <<
        "<tr><th class=left>Sequence<td class=left style='font-family:monospace;white-space:nowrap'>";
    for(uint64_t position=0; position<consensus.size(); position++) {
        const auto& [base, coverage] = consensus[position];
        html << "<span title='Position " << position <<
            " coverage " << coverage << "'";
        html << ">";
        html << base;
        html << "</span>";

    }

    // Write a line with coverage at each position.
    std::map<char, uint64_t> coverageLegend;
    html <<
        "<tr><th class=left>Coverage<td class=left style='font-family:monospace;white-space:nowrap'>";
    for(uint64_t position=0; position<consensus.size(); position++) {
        const auto& [ignore, coverage] = consensus[position];
        const char coverageCharacter = getCoverageCharacter(coverage);
        if((coverage > 9) and (coverage < 36)) {
            coverageLegend.insert(make_pair(coverageCharacter, coverage));
        }
        html << "<span title='Position " << position <<
            " coverage " << coverage << "'";
        html << ">";
        html << coverageCharacter;
        html << "</span>";
    }
    html << "</table>";

    // Write coverage legend.
    html << "<h3>Coverage legend</h3>"
        "<table><tr><th>Character<th>Coverage";
    if(minCoverage < 10) {
        html << "<tr><td class=centered>1-9<td class=centered>As displayed";
    }
    for(const auto& [character, coverage]: coverageLegend) {
        html << "<tr><td class=centered>" << character <<
            "<td class=centered>" << coverage;
    }
    if(maxCoverage >= 36) {
        html << "<tr><td class=centered>*<td class=centered>&ge;36";
    }
    html << "</table>";
}



char LocalAssembly7::getCoverageCharacter(uint64_t coverage)
{
    if(coverage < 10) {
        return char(coverage + '0');
    } else if(coverage < 36) {
        return char(coverage - 10 + 'A');
    } else {
        return '*';
    }
}


void LocalAssembly7::writeAlignment(
    const vector< vector<AlignedBase> >& alignment,
    const vector<AlignedBase>& alignedConsensus,
    const vector< pair<Base, uint64_t> >& consensus,
    const vector< pair<uint64_t, uint64_t> > sequenceIdsWithWeight
    ) const
{
     if(not html) {
        return;
    }

    const uint64_t n = alignment.size();
    SHASTA2_ASSERT(sequenceIdsWithWeight.size() == n);

    html <<
        "<h3>Alignment</h3>"
        "<table>"
        "<tr><th>Sequence<br>id<th>Weight<th>Length<th class=left>Sequence";

    uint64_t totalWeight = 0;
    for(uint64_t i=0; i<n; i++) {
        const vector<AlignedBase>& alignmentRow = alignment[i];
        SHASTA2_ASSERT(alignmentRow.size() == alignedConsensus.size());
        const auto& [sequenceId, weight] = sequenceIdsWithWeight[i];
        totalWeight += weight;
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        html <<
            "<tr>"
            "<td class=centered>" << sequenceId <<
            "<td class=centered>" << weight <<
            "<td class=centered>" << sequenceInfo.sequence.size() <<
            "<td class=left style='font-family:monospace;white-space:nowrap'>";
        for(uint64_t position=0; position<alignmentRow.size(); position++) {
            const AlignedBase b = alignmentRow[position];
            const AlignedBase bConsensus = alignedConsensus[position];
            const bool isDiscordant = (b != bConsensus);
            if(isDiscordant) {
                html << "<span style='background-color:Pink'>";
            }
            html << b;
            if(isDiscordant) {
                html << "</span>";
            }
        }
        std::ranges::copy(alignmentRow, ostream_iterator<AlignedBase>(html));
    }

    // Add a row with aligned consensus.
    html <<
        "<tr style='background-color:LightCyan'>"
        "<th class=left>Consensus" <<
        "<td class=centered>" << totalWeight <<
        "<td class=centered>" << consensus.size() <<
        "<td class=left style='font-family:monospace;white-space:nowrap'>";
    uint64_t positionInConsensus = 0;
    for(uint64_t position=0; position<alignedConsensus.size(); position++) {
        const AlignedBase b = alignedConsensus[position];
        if(b.isGap()) {
            html << "<span style='background-color:Pink'>-</span>";
        } else {
            const uint64_t coverage = consensus[positionInConsensus].second;
            html << "<span title='Position " << positionInConsensus <<
                ", coverage " << coverage << "'";
            html << ">" << b << "</span>";
            ++positionInConsensus;
        }
    }

    // Add a row with  aligned consensus coverage.
    html <<
        "<tr>"
        "<th class=left>Coverage" <<
        "<td class=centered>" << totalWeight <<
        "<td class=centered>" << consensus.size() <<
        "<td class=left style='font-family:monospace;white-space:nowrap'>";
    positionInConsensus = 0;
    for(uint64_t position=0; position<alignedConsensus.size(); position++) {
        const AlignedBase b = alignedConsensus[position];
        if(b.isGap()) {
            html << "<span style='background-color:Pink'>-</span>";
        } else {
            const uint64_t coverage = consensus[positionInConsensus].second;
            const char coverageCharacter = getCoverageCharacter(coverage);
            html << "<span title='Position " << positionInConsensus <<
                ", coverage " << coverage << "'";
            html << ">" << coverageCharacter << "</span>";
            ++positionInConsensus;
        }
    }

    html << "</table>";
}




void LocalAssembly7::writeSequence() const
{
    if(not html) {
        return;
    }

    if(success) {

        // Write the consensus.
        html <<
            "<h3>Assembled sequence</h3>"
            "<table>"
            "<tr><th class=left>Length<td class=left>" << sequence.size() <<
            "<tr><th class=left>Sequence<td class=left style='font-family:monospace;white-space:nowrap'>";
        std::ranges::copy(sequence, ostream_iterator<Base>(html));
        html << "</table>";

    } else {

        html << "<br>The local assembly failed." << endl;
    }

}



void LocalAssembly7::removeOutliers()
{
    // Gather the offsets.
    vector<pair<uint64_t, uint64_t > > offsetTable; // (index in orientedReadInfos, offset).
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        const OrientedRead& orientedRead = orientedReads[i];
        if(orientedRead.isOnBothAnchors()) {
            const uint64_t offset = orientedRead.positionOffsetAB();
            offsetTable.push_back({i, offset});
        }
    }
    sort(offsetTable.begin(), offsetTable.end(), OrderPairsBySecondOnly<uint64_t, uint64_t>());

    // Find places where there is an unreasonably jump in the offset.
    vector<uint64_t> violations(1, 0);
    for(uint64_t i1=1; i1<offsetTable.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t offset0 = offsetTable[i0].second;
        const uint64_t offset1 = offsetTable[i1].second;
        if(not checkOffsets(offset0, offset1)) {
            violations.push_back(i1);
            // cout << "Violation " << offset0 << " " << offset1 << " " << i1 << endl;
        }
    }
    violations.push_back(offsetTable.size());

    // If no violations were found, keep all the OrientedReadInfos.
    // This is the most common case.
    if(violations.size() == 2) {
        return;
    }

#if 0
    cout << "violations vector ";
    std::ranges::copy(violations, ostream_iterator<uint64_t>(cout, " "));
    cout << endl;
#endif

    // Find the largest interval between violations.
    uint64_t keepBegin = 0;
    uint64_t keepEnd = 0;
    for(uint64_t i1=1; i1<violations.size(); i1++) {
        const uint64_t i0 = i1 - 1;
        const uint64_t violation0 = violations[i0];
        const uint64_t violation1 = violations[i1];
        if(violation1 - violation0 > keepEnd - keepBegin) {
            keepBegin = violation0;
            keepEnd = violation1;
        }
    }
    // cout << "keepBegin " << keepBegin << ", keepEnd " << keepEnd << endl;

    // Only keep OrientedReads that are at positions [keepBegin, keepEnd)
    // in the offset table.
    std::set<uint64_t> discard;
    for(uint64_t i=0; i<keepBegin; i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReads[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }
    for(uint64_t i=keepEnd; i<offsetTable.size(); i++) {
        const uint64_t j = offsetTable[i].first;
        discard.insert(j);
        if(html) {
            html << "<br>Discarding " << orientedReads[j].orientedReadId <<
                " due to inconsistent offsets.";
        }
    }

    vector<OrientedRead> newOrientedReads;
    for(uint64_t i=0; i<orientedReads.size(); i++) {
        if(not discard.contains(i)) {
            newOrientedReads.push_back(orientedReads[i]);
        }
    }
    orientedReads.swap(newOrientedReads);
}



bool LocalAssembly7::checkOffsets(uint64_t offset0, uint64_t offset1)
{
    // EXPOSE WHEN CODE STABILIZES.
    const double aDrift = 0.02;
    const double bDrift = 100.;

    if(offset1 == offset0) {
        return true;
    }

    SHASTA2_ASSERT(offset1 > offset0);

    const double average = 0.5 * double(offset0 + offset1);
    const uint64_t difference = offset1 - offset0;

    const double acceptableDifference = aDrift * average + bDrift;

    return difference < uint64_t(std::round(acceptableDifference));

}



void LocalAssembly7::writeKmerOccurrences(const Graph& graph, const string& fileName) const
{
    ofstream csv(fileName);
    writeKmerOccurrences(graph, csv);
}



void LocalAssembly7::writeKmerOccurrences(const Graph& graph, ostream& csv) const
{
    csv << "VertexId,KmerId,Coverage,SequenceId,Position in sequence,\n";

    BGL_FORALL_VERTICES(v, graph, Graph) {
        const Vertex& vertex = graph[v];
        for(const KmerOccurrence& kmerOccurrence: vertex.occurrences) {
            const uint64_t positionInSequence = kmerOccurrence.position;

            csv << vertex.vertexId << ",";
            csv << vertex.kmerId << ",";
            csv << vertex.coverage << ",";
            csv << kmerOccurrence.sequenceId << ",";
            csv << positionInSequence << ",";
            csv << "\n";
        }
    }
}



void LocalAssembly7::runAbpoa()
{
    // Get the sequenceIds to be used, sorted in order of decreasing coverage.
    vector<uint64_t> sequenceIds;
    getSequencesOnBothAnchors(sequenceIds);

    if(html) {
        html <<
            "<h3>Local assembly with abpoa</h3>"
            "The local assembly will use the following "
            "sequences of oriented reads on both anchors, "
            "presented to abpoa in this order of decreasing coverage."
            "<br>Each sequence is presented to abpoa a number of times "
            "equal to its coverage, with (implicit) weight 1."
            "<br><br><table>"
            "<tr><th>Sequence<br>id<th>Coverage<th>Length";
        for(const uint64_t sequenceId: sequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            html <<
                "<tr>"
                "<td class=centered>" << sequenceId <<
                "<td class=centered>" << sequenceInfo.coverage() <<
                "<td class=centered>" << sequenceInfo.sequence.size();
        }

        html << "</table>";
    }


    // Abpoa doesn't support weights, so we have to enter each sequence
    // a number of times equal to its coverage.
    vector< vector<Base> > msaSequences;
    vector< pair<uint64_t, uint64_t> > msaSequenceIdsWithWeight;
    for(const uint64_t sequenceId: sequenceIds) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        for(uint64_t i=0; i<sequenceInfo.coverage(); i++) {
            msaSequences.push_back(sequenceInfo.sequence);
            msaSequenceIdsWithWeight.push_back(make_pair(sequenceId, 1));
        }
    }
    if(html) {
        html << "<br>Total coverage for abpoa is " << msaSequences.size() << ".";
    }

    // Run abpoa.
    vector< pair<Base, uint64_t> > consensus;
    vector< vector<AlignedBase> > alignment;
    vector<AlignedBase> alignedConsensus;
    const auto t0 = steady_clock::now();
    // The alignment is normally computed only for the html display. It is
    // also needed when the repair below is requested.
    const bool computeAlignment = bool(html) or options.useMsaRepair;
    abpoa(msaSequences, consensus, alignment, alignedConsensus, computeAlignment);
    const auto t1 = steady_clock::now();
    SHASTA2_ASSERT(alignment.size() == msaSequenceIdsWithWeight.size());

    // Repair the bad homopolymer regions of the alignment, if requested (see
    // Options::useMsaRepair). Every row here spans the whole alignment - abpoa
    // takes no anchoring information - and each entered sequence
    // already stands for one unit of coverage, so every row is anchored on
    // both sides and votes with weight 1.
    uint64_t repairedRegionCount = 0;
    const auto t2 = steady_clock::now();
    if(options.useMsaRepair) {
        const vector<uint64_t> weights(alignment.size(), 1);

        // Each row is one read, in the order entered above.
        vector< array<uint64_t, 2> > strandWeights;
        if(homopolymerModelPointer) {
            for(const uint64_t sequenceId: sequenceIds) {
                for(const OrientedReadId orientedReadId: sequences[sequenceId].orientedReadIds) {
                    array<uint64_t, 2> strandWeight = {0, 0};
                    strandWeight[orientedReadId.getStrand()] = 1;
                    strandWeights.push_back(strandWeight);
                }
            }
        }

        MsaRepairOptions msa1Options;
        msa1Options.homopolymerModelPointer = homopolymerModelPointer;
        repairedRegionCount = msaRepair(alignment, alignedConsensus, consensus, weights, {},
            msa1Options, strandWeights);
    }
    const auto t3 = steady_clock::now();

    if(html) {
        html << "<br>Abpoa completed in " << seconds(t1-t0) << " seconds.";
        if(options.useMsaRepair) {
            html << "<br>Msa repair completed in " << seconds(t3-t2) <<
                " seconds and rebuilt " << repairedRegionCount << " region(s) of the alignment.";
        }
        writeAlignment(alignment, alignedConsensus, consensus, msaSequenceIdsWithWeight);
        writeConsensus(consensus);
    }

    // Store the sequence.
    for(const auto& [b, ignore]: consensus) {
        sequence.push_back(b);
    }
    success = true;
}



// If useAll is true, this uses oriented reads that
// appear only in the left or right anchor.
void LocalAssembly7::runTheseus(bool useAll)
{
    // Get the sequenceIds to be used, sorted in order of decreasing coverage.
    vector<uint64_t> bothSidesFixedSequenceIds;
    getSequencesOnBothAnchors(bothSidesFixedSequenceIds);
    vector<uint64_t> leftFixedSequenceIds;
    vector<uint64_t> rightFixedSequenceIds;
    if(useAll) {
        getSequencesOnAnchorA(leftFixedSequenceIds);
        getSequencesOnAnchorB(rightFixedSequenceIds);
    }

    if(html) {
        html <<
            "<h3>Local assembly with Theseus</h3>"
            "The local assembly will use the following "
            "sequences of oriented reads on both anchors, "
            "presented to Theseus in this order."
            "<br><br><table>"
            "<tr><th>Sequence<br>id<th>On<br>A<th>On<br>B<th>Coverage<th>Length";
        for(const uint64_t sequenceId: bothSidesFixedSequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            html <<
                "<tr>"
                "<td class=centered>" << sequenceId <<
                "<td class=centered>&check;" <<
                "<td class=centered>&check;" <<
                "<td class=centered>" << sequenceInfo.coverage() <<
                "<td class=centered>" << sequenceInfo.sequence.size();
        }
        for(const uint64_t sequenceId: leftFixedSequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            html <<
                "<tr>"
                "<td class=centered>" << sequenceId <<
                "<td class=centered>&check;" <<
                "<td class=centered>" <<
                "<td class=centered>" << sequenceInfo.coverage() <<
                "<td class=centered>" << sequenceInfo.sequence.size();
        }
        for(const uint64_t sequenceId: rightFixedSequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            html <<
                "<tr>"
                "<td class=centered>" << sequenceId <<
                "<td class=centered>" <<
                "<td class=centered>&check;" <<
                "<td class=centered>" << sequenceInfo.coverage() <<
                "<td class=centered>" << sequenceInfo.sequence.size();
        }

        html << "</table>";
    }



    // Gather the sequences to be passed to theseus.
    uint64_t totalWeight = 0;
    vector< pair<uint64_t, uint64_t> > msaSequenceIdsWithWeight;
    vector< pair<vector<Base>, uint64_t> > bothSidesFixedSequences;
    for(const uint64_t sequenceId: bothSidesFixedSequenceIds) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        const uint64_t coverage = sequenceInfo.coverage();
        bothSidesFixedSequences.push_back(make_pair(sequenceInfo.sequence, coverage));
        msaSequenceIdsWithWeight.push_back(make_pair(sequenceId, coverage));
        totalWeight += coverage;
    }
    vector< pair<vector<Base>, uint64_t> > leftFixedSequences;
    vector< pair<vector<Base>, uint64_t> > rightFixedSequences;
    if(useAll) {
        for(const uint64_t sequenceId: leftFixedSequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            const uint64_t coverage = sequenceInfo.coverage();
            leftFixedSequences.push_back(make_pair(sequenceInfo.sequence, coverage));
            msaSequenceIdsWithWeight.push_back(make_pair(sequenceId, coverage));
            totalWeight += coverage;
        }
        for(const uint64_t sequenceId: rightFixedSequenceIds) {
            const SequenceInfo& sequenceInfo = sequences[sequenceId];
            const uint64_t coverage = sequenceInfo.coverage();
            rightFixedSequences.push_back(make_pair(sequenceInfo.sequence, coverage));
            msaSequenceIdsWithWeight.push_back(make_pair(sequenceId, coverage));
            totalWeight += coverage;
        }
    }
    if(html) {
        html << "<br>Total coverage for Theseus is " << totalWeight << ".";
        theseusWriteFile(bothSidesFixedSequences, leftFixedSequences, rightFixedSequences,
            "Pericles.fasta");
    }

    // Run Theseus. The alignment is normally computed only for the html
    // display. It is also needed when the repair below is requested.
    vector< pair<Base, uint64_t> > consensus;
    vector<AlignedBase> alignedConsensus;
    vector< vector<AlignedBase> > alignment;
    const bool computeAlignment = bool(html) or options.useMsaRepair;
    const auto t0 = steady_clock::now();
    theseus(
        bothSidesFixedSequences, leftFixedSequences, rightFixedSequences,
        consensus, alignment, alignedConsensus, computeAlignment);
    const auto t1 = steady_clock::now();
    if(computeAlignment) {
        SHASTA2_ASSERT(alignment.size() == msaSequenceIdsWithWeight.size());
    }

    // Repair the bad homopolymer regions of the alignment, if requested (see
    // Options::useMsa1).
    uint64_t repairedRegionCount = 0;
    const auto t2 = steady_clock::now();
    if(options.useMsaRepair) {

        // Theseus returns the rows in the order the groups were given, so how
        // each row is anchored is known here rather than guessed from its
        // gaps. That matters: a read fixed on one side only does not reach
        // across the whole alignment, and Theseus pads the part it does not
        // reach with the same '-' used for a deletion. Counting that padding
        // as a deletion would delete bases every read covering them agrees on.
        vector<Anchoring> anchoring;
        anchoring.reserve(alignment.size());
        anchoring.insert(anchoring.end(), bothSidesFixedSequences.size(),
            Anchoring::BothSides);
        anchoring.insert(anchoring.end(), leftFixedSequences.size(),
            Anchoring::LeftOnly);
        anchoring.insert(anchoring.end(), rightFixedSequences.size(),
            Anchoring::RightOnly);
        SHASTA2_ASSERT(anchoring.size() == alignment.size());

        // Unlike abpoa, Theseus takes a weight per sequence, so each read
        // appears once and votes with its coverage.
        vector<uint64_t> weights;
        weights.reserve(alignment.size());
        for(const auto& [sequenceId, weight]: msaSequenceIdsWithWeight) {
            weights.push_back(weight);
        }

        // A row stands for all the reads with its sequence, which can be on
        // either strand.
        vector< array<uint64_t, 2> > strandWeights;
        if(homopolymerModelPointer) {
            for(const auto& [sequenceId, weight]: msaSequenceIdsWithWeight) {
                array<uint64_t, 2> strandWeight = {0, 0};
                for(const OrientedReadId orientedReadId: sequences[sequenceId].orientedReadIds) {
                    ++strandWeight[orientedReadId.getStrand()];
                }
                strandWeights.push_back(strandWeight);
            }
        }

        MsaRepairOptions msa1Options;
        msa1Options.homopolymerModelPointer = homopolymerModelPointer;
        repairedRegionCount = msaRepair(alignment, alignedConsensus, consensus, weights, anchoring,
            msa1Options, strandWeights);
    }
    const auto t3 = steady_clock::now();

    if(html) {
        html << "<br>Theseus completed in " << seconds(t1-t0) << " seconds.";
        if(options.useMsaRepair) {
            html << "<br>Msa repair completed in " << seconds(t3-t2) <<
                " seconds and rebuilt " << repairedRegionCount << " region(s) of the alignment.";
        }
        writeAlignment(alignment, alignedConsensus, consensus, msaSequenceIdsWithWeight);
        writeConsensus(consensus);
    }

    // Store the sequence.
    for(const auto& [b, ignore]: consensus) {
        sequence.push_back(b);
    }
    success = true;
}



void LocalAssembly7::runAdaptive()
{
    // Try fast path first, if allowed.
    if(options.allowFastPath) {
        runFastPath();
        if(success) {
            return;
        }
    }


    // See if we have enough oriented reads on both anchors.
    if(getTotalCommonCoverage() >= options.commonCoverageThreshold) {

        // We have enough oriented reads on both anchors.
        // Only use those oriented reads on both anchors.

        if(getMaxLengthCommon() <= options.maxAbpoaLength) {
            // The MSA is not too long, use abpoa.
            runAbpoa();
        } else {
            // The MSA is long, use theseus,
            // using only oriented reads on both anchors.
            runTheseus(false);
        }

    } else {

        // We don'thave enough oriented reads on both anchors,
        // so we also use oriented reads that are on just one anchor.
        runTheseus(true);

    }
}



void LocalAssembly7::Options::setMethod(const string& s)
{
    if(s == "Adaptive") {
        method = Method::Adaptive;
    } else if(s == "Adaptive") {
        method = Method::Adaptive;
    } else if(s == "Abpoa") {
        method = Method::Abpoa;
    } else if(s == "TheseusOnly") {
        method = Method::TheseusOnly;
    } else if(s == "TheseusAll") {
        method = Method::TheseusAll;
    } else if(s == "DeBruijn") {
        method = Method::DeBruijn;
    } else {
        method = Method::Invalid;
    }
}



// Get  pairs(sequenceId, coverage) for the SequenceInfos on both anchors,
// sorted by decreasing coverage.
// These are used for assembly with abpoa or theseus.
// SequenceId is the index in the sequences vector.
void LocalAssembly7::getSequencesOnBothAnchors(
    vector<uint64_t>& v) const
{
    // Get pairs (sequenceId, coverage) and sort them by decreasing coverage.
    vector< pair<uint64_t, uint64_t> > x;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        if(sequenceInfo.isOnAnchorA and sequenceInfo.isOnAnchorB) {
            x.push_back(make_pair(sequenceId, sequenceInfo.coverage()));
        }
    }
    sort(x.begin(), x.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Costruct the vector with the sequenceIds only.
    v.clear();
    for(const auto& [sequenceId, ignore]: x) {
        v.push_back(sequenceId);
    }
}



// Same, but for SequenceInfos on one anchor only.
// These are used for local assembly with theseus.
void LocalAssembly7::getSequencesOnAnchorA(
    vector<uint64_t>& v) const
{
    // Get pairs (sequenceId, coverage) and sort them by decreasing coverage.
    vector< pair<uint64_t, uint64_t> > x;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        if(sequenceInfo.isOnAnchorA and (not sequenceInfo.isOnAnchorB)) {
            x.push_back(make_pair(sequenceId, sequenceInfo.coverage()));
        }
    }
    sort(x.begin(), x.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Costruct the vector with the sequenceIds only.
    v.clear();
    for(const auto& [sequenceId, ignore]: x) {
        v.push_back(sequenceId);
    }
}



void LocalAssembly7::getSequencesOnAnchorB(
    vector<uint64_t>& v) const
{
    // Get pairs (sequenceId, coverage) and sort them by decreasing coverage.
    vector< pair<uint64_t, uint64_t> > x;
    for(uint64_t sequenceId=0; sequenceId<sequences.size(); sequenceId++) {
        const SequenceInfo& sequenceInfo = sequences[sequenceId];
        if(not(sequenceInfo.isOnAnchorA) and sequenceInfo.isOnAnchorB) {
            x.push_back(make_pair(sequenceId, sequenceInfo.coverage()));
        }
    }
    sort(x.begin(), x.end(), OrderPairsBySecondOnlyGreater<uint64_t, uint64_t>());

    // Costruct the vector with the sequenceIds only.
    v.clear();
    for(const auto& [sequenceId, ignore]: x) {
        v.push_back(sequenceId);
    }

}

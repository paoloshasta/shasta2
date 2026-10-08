// Shasta2.
#include "StrandSeparation1.hpp"
#include "AssemblyGraph.hpp"
#include "color.hpp"
#include "computeLayout.hpp"
#include "deduplicate.hpp"
#include "DisjointSets.hpp"
#include "graphvizToHtml.hpp"
#include "html.hpp"
#include "weightedShuffle.hpp"
using namespace shasta2;
using namespace StrandSeparation1;

// Standard library.
#include <iomanip>
#include <random>


// StrandContact constructor.
// The strandContactVerticesmust be sorted by id.
// The debugOutputBaseName and strandContactId are only used for debug output.
StrandContact::StrandContact(
    AssemblyGraph& assemblyGraph,
    const vector<AssemblyGraph::vertex_descriptor>& strandContactVertices,
    const string& debugOutputBaseName,
    uint64_t strandContactId
    ) :
    assemblyGraph(assemblyGraph),
    strandContactVertices(strandContactVertices),
    debugOutputBaseName(debugOutputBaseName),
    strandContactId(strandContactId)
{
    const bool debug = false;
    if(debug) {
        html.open(debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + ".html");
        cout << "Working on strand contact " << strandContactId << endl;
        writeHtmlBegin(html, "Strand contact " + to_string(strandContactId));
        writeMakeAllTablesCopyable(html);
        html <<
            "</head>"
            "<body onload='makeAllTablesCopyable()'>";
        html << "<h1>Strand contact " << strandContactId << "</h1>";
    }

    // Create the BipartiteGraph.
    gatherSegmentPairs();
    countReadOccurrences();
    createBipartiteGraph();

    // Strand separation in the BipartiteGraph.
    Split split;
    computeSplit(split);

    if(html) {
        html << "<h2>Best strand separation found</h2>";
        writeSplitSummary(split);
        writeSplitDetails(split);
        writeBipartiteGraphGraphviz(split);
        writeBipartiteGraphCustom(split);
    }

    storeSegmentInformation(split);
    updateAssemblyGraph();

    if(html) {
        writeHtmlEnd(html);
        cout << "Done working on strand contact " << strandContactId << endl;
    }
}


uint64_t StrandContact::id(AssemblyGraphBaseClass::vertex_descriptor v) const
{
    return assemblyGraph.id(v);
}



// Store component information in the SegmentPairs.
void StrandContact::storeSegmentInformation(const Split& split)
{
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        SegmentPair& segmentPair = segmentPairs[segmentPairId];
        for(uint64_t segmentIndexInPair=0; segmentIndexInPair<2; segmentIndexInPair++) {
            SegmentInfo& segmentInfo = segmentPair.segmentInfos[segmentIndexInPair];
            const BipartiteGraph::vertex_descriptor v = segmentInfo.v;
            segmentInfo.componentId = split.vertexComponent[v];
        }
        SHASTA2_ASSERT(segmentPair.segmentInfos[0].componentId == (segmentPair.segmentInfos[1].componentId ^ 1));

        const BipartiteGraphVertexStatistics statistics =
            bipartiteGraph.getVertexStatistics(split, segmentPair.segmentInfos[0].v);
        segmentPair.crossStrandEdgeFrequencyRatio = statistics.crossStrandEdgeFrequencyRatio();
    }


    if(not html) {
        return;
    }

    html <<
        std::fixed << std::setprecision(3) <<
        "<h2>Segment information for the best strand separation found</h2>"
        "<table>"
        "<tr>"
        "<th>Segment<br>pair id"
        "<th>Segment 0"
        "<th>Component 0"
        "<th>Segment 1"
        "<th>Component 1"
        "<th>Cross-strand<br>edge<br>frequency<br>ratio"
        "<th>Ambiguous";
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        const SegmentPair& segmentPair = segmentPairs[segmentPairId];
        const bool isAmbiguous = (segmentPair.crossStrandEdgeFrequencyRatio > maxCrossStrandFrequencyRatio);
        const string color0 = (isAmbiguous ? ambiguousColor() : componentColor(segmentPair.segmentInfos[0].componentId));
        const string color1 = (isAmbiguous ? ambiguousColor() : componentColor(segmentPair.segmentInfos[1].componentId));
        html <<
            "<tr>"
            "<td class=centered>" << segmentPairId <<
            "<td class=centered style='background-color:" << color0 << "'>" << segmentPair.segmentInfos[0].id <<
            "<td class=centered>" << segmentPair.segmentInfos[0].componentId <<
            "<td class=centered style='background-color:" << color1 << "'>" << segmentPair.segmentInfos[1].id <<
            "<td class=centered>" << segmentPair.segmentInfos[1].componentId <<
            "<td class=centered>";
        if(segmentPair.crossStrandEdgeFrequencyRatio > 0.) {
            html << segmentPair.crossStrandEdgeFrequencyRatio;
        }
        html << "<td class=centered>";
        if(isAmbiguous) {
            html << "&check;";
        }
    }

    html << "</table>";



    // Write a csv file that can be loaded in Bandage to show the Split.
    const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-Split-Bandage.csv";
    ofstream csv(fileName);
    csv << "Segment,Component,Color\n";
    for(const SegmentPair& segmentPair: segmentPairs) {
        const bool isAmbiguous = (segmentPair.crossStrandEdgeFrequencyRatio > maxCrossStrandFrequencyRatio);
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            const uint64_t componentId = segmentInfo.componentId;
            csv << segmentInfo.id << ",";

            if(isAmbiguous) {
                csv << "?";
            } else {
                csv << componentId;
            }
            csv << ",";

            if(isAmbiguous) {
                csv << ambiguousColor();
            } else {
                csv << componentColor(componentId);
            }
            csv << ",";

            csv << "\n";
        }
    }



    // Write a gfa file containing only the segments that will
    // be used in the rest of the strand separation process
    // for this StrandContact.
    // These are segments that belong to an even number component
    // or are flagged as ambiguous.
    vector<Segment> segmentsForOutput;
    for(const SegmentPair& segmentPair: segmentPairs) {
        const bool isAmbiguous = (segmentPair.crossStrandEdgeFrequencyRatio > maxCrossStrandFrequencyRatio);
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            const uint64_t componentId = segmentInfo.componentId;
            if(isAmbiguous or ((componentId % 2) == 0)) {
                segmentsForOutput.push_back(segmentInfo.segment);
            }
        }
    }
    std::ranges::sort(segmentsForOutput, assemblyGraph.orderById);

    const string gfaFileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-Split-Bandage.gfa";
    assemblyGraph.writeGfa(gfaFileName, segmentsForOutput);

}



uint64_t StrandContact::id(Segment segment) const
{
    return assemblyGraph.id(segment);
}



void StrandContact::gatherSegmentPairs()
{
    SHASTA2_ASSERT(std::ranges::is_sorted(strandContactVertices, assemblyGraph.orderById));

    // Out-edges of the strandContactVertices give us Segments
    // internal to the StrandContact plus the exits.
    for(const AssemblyGraph::vertex_descriptor v0: strandContactVertices) {
        BGL_FORALL_OUTEDGES(v0, segment, assemblyGraph, AssemblyGraph) {
            const Segment segmentRc = assemblyGraph[segment].eRc;
            SHASTA2_ASSERT(segmentRc != segment);

            if(id(segment) < id(segmentRc)) {
                const AssemblyGraph::vertex_descriptor v1 = target(segment, assemblyGraph);
                const bool isExit = not std::ranges::binary_search(strandContactVertices, v1, assemblyGraph.orderById);
                SegmentPair& segmentPair = segmentPairs.emplace_back();

                SegmentInfo& segmentInfo0 = segmentPair.segmentInfos[0];
                SegmentInfo& segmentInfo1 = segmentPair.segmentInfos[1];

                segmentInfo0.segment = segment;
                segmentInfo0.id = id(segment);
                segmentInfo0.isExit = isExit;
                segmentInfo1.segment = segmentRc;
                segmentInfo1.id = id(segmentRc);
                segmentInfo1.isEntrance = isExit;

                segmentPair.length = assemblyGraph[segment].length();
                SHASTA2_ASSERT(segmentPair.length == assemblyGraph[segmentRc].length());

                segmentPair.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
                SHASTA2_ASSERT(segmentPair.coverage == assemblyGraph[segmentRc].lengthWeightedAverageCoverage());
            }
        }
    }



    // In-edges of the strandContactVertices give us the entrances.
    for(const AssemblyGraph::vertex_descriptor v0: strandContactVertices) {
        BGL_FORALL_INEDGES(v0, segment, assemblyGraph, AssemblyGraph) {
            const Segment segmentRc = assemblyGraph[segment].eRc;
            SHASTA2_ASSERT(segmentRc != segment);

            if(id(segment) < id(segmentRc)) {
                const AssemblyGraph::vertex_descriptor v1 = source(segment, assemblyGraph);
                const bool isEntrance = not std::ranges::binary_search(strandContactVertices, v1, assemblyGraph.orderById);

                if(isEntrance) {
                    SegmentPair& segmentPair = segmentPairs.emplace_back();

                    SegmentInfo& segmentInfo0 = segmentPair.segmentInfos[0];
                    SegmentInfo& segmentInfo1 = segmentPair.segmentInfos[1];

                    segmentInfo0.segment = segment;
                    segmentInfo0.id = id(segment);
                    segmentInfo0.isEntrance = true;
                    segmentInfo1.segment = segmentRc;
                    segmentInfo1.id = id(segmentRc);
                    segmentInfo1.isExit = true;

                    segmentPair.length = assemblyGraph[segment].length();
                    SHASTA2_ASSERT(segmentPair.length == assemblyGraph[segmentRc].length());

                    segmentPair.coverage = assemblyGraph[segment].lengthWeightedAverageCoverage();
                    SHASTA2_ASSERT(segmentPair.coverage == assemblyGraph[segmentRc].lengthWeightedAverageCoverage());
                }
            }
        }
    }

    std::ranges::sort(segmentPairs, {}, &SegmentPair::id0);
    writeSegmentPairs();
}



void StrandContact::writeSegmentPairs()
{
    if(not html) {
        return;
    }

    html << std::fixed << std::setprecision(1);
    html << "<h2>Segment pairs</h2>"
        "<table>"
        "<tr>"
        "<th>Segment<br>pair id"
        "<th>Segment 0"
        "<th>Segment 1"
        "<th>Length"
        "<th>Coverage"
        "<th>Segment 0<br>is entrance"
        "<th>Segment 0<br>is exit"
        "<th>Segment 1<br>is entrance"
        "<th>Segment 1<br>is exit";

    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        const SegmentPair& segmentPair = segmentPairs[segmentPairId];
        html <<
            "<tr>"
            "<td class=centered>" << segmentPairId <<
            "<td class=centered>" << segmentPair.segmentInfos[0].id <<
            "<td class=centered>" << segmentPair.segmentInfos[1].id <<
            "<td class=centered>" << segmentPair.length <<
            "<td class=centered>" << segmentPair.coverage;

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[0].isEntrance) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[0].isExit) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[1].isEntrance) {
            html << "&check;";
        }

        html << "<td class=centered>";
        if(segmentPair.segmentInfos[1].isExit) {
            html << "&check;";
        }
    }

    html << "</table>";



    // Write a csv file that can be loaded in Bandage to show this StrandContact
    // with its entrances ane exits.
    const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-Bandage.csv";
    ofstream csv(fileName);
    csv << "Segment,Classification,Color\n";
    for(const SegmentPair& segmentPair: segmentPairs) {
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            csv << segmentInfo.id << ",";
            if(segmentInfo.isEntrance) {
                SHASTA2_ASSERT(not segmentInfo.isExit);
                csv << "Entrance,";
                csv << hslToRgbString(0.333, 0.5, 0.6) << ",";  // Green
            } else if(segmentInfo.isExit) {
                SHASTA2_ASSERT(not segmentInfo.isEntrance);
                csv << "Exit,";
                csv << hslToRgbString(0., 0.5, .6) << ",";      // Red
            } else {
                csv << "Internal,";
                csv << hslToRgbString(0.6, 0.5, .6) << ",";     // Blue
            }
            csv << endl;
        }
    }

}



void StrandContact::countReadOccurrences()
{
    // Count occurrences of reads in the first Segment of each SegmentPair.
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        const SegmentPair& segmentPair = segmentPairs[segmentPairId];

        // Use the first Segment of the SegmentPair.
        const Segment segment = segmentPair.segmentInfos[0].segment;

        for(const AssemblyGraphEdgeStep& step: assemblyGraph[segment]) {
            for(const OrientedReadId orientedReadId: step.anchorPair.orientedReadIds) {
                const ReadId readId = orientedReadId.getReadId();
                const Strand strand = orientedReadId.getStrand();
                readOccurrenceMap[readId].emplace_back(ReadOccurrence(segmentPairId, strand));
            }
        }
    }

    // Deduplicate and count the occurrences for each ReadId.
    vector<uint64_t> count;
    for(auto&[readId, occurrences]: readOccurrenceMap) {
        deduplicateAndCount(occurrences, count);
        for(uint64_t i=0; i<occurrences.size(); i++) {
            occurrences[i].frequency = count[i];
        }
    }

    // Remove from the map reads that occur in just one Segment.
    for(auto it=readOccurrenceMap.begin(); it!=readOccurrenceMap.end(); /* Increment later */) {
        auto itNext = it;
        ++itNext;
        if(it->second.size() == 1) {
            readOccurrenceMap.erase(it);
        }
        it = itNext;
    }

    for(const auto&[readId, occurrences]: readOccurrenceMap) {
        SHASTA2_ASSERT(occurrences.size() > 1);
    }

    writeReadOccurrences();

}



void StrandContact::writeReadOccurrences()
{
    if(not html) {
        return;
    }

    const string fileName = debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-ReadOccurrences.csv";
    ofstream csv(fileName);
    csv << "ReadId,Strand,Segment,Frequency\n";
    for(const auto&[readId, occurrences]: readOccurrenceMap) {
        for(const auto& occurrence: occurrences) {
            const Segment segment = segmentPairs[occurrence.segmentPairId].segmentInfos[0].segment;
            csv << readId << ",";
            csv << occurrence.strand << ",";
            csv << id(segment) << ",";
            csv << occurrence.frequency << "\n";
        }
    }

    html << "<h2>Read occurrences</h2>"
        "For details of occurrences of reads in the segments, see "
        "<a href='"<< fileName << "'>" << fileName << "</a>.";
}



void StrandContact::createBipartiteGraph()
{
    // Create the vertices corresponding to Segments.
    for(uint64_t segmentPairId=0; segmentPairId<segmentPairs.size(); segmentPairId++) {
        SegmentPair& segmentPair = segmentPairs[segmentPairId];
        for(uint64_t segmentIndexInPair=0; segmentIndexInPair<2; segmentIndexInPair++) {
            SegmentInfo& segmentInfo = segmentPair.segmentInfos[segmentIndexInPair];
            segmentInfo.v = bipartiteGraph.addVertex(segmentPairId, segmentIndexInPair);
        }
    }

    // Create the vertices corresponding to OrientedReadIds.
    for(const auto&[readId, ignore]: readOccurrenceMap) {
        for(Strand strand=0; strand<2; strand++) {
            const OrientedReadId orientedReadId(readId, strand);
            const BipartiteGraph::vertex_descriptor v = bipartiteGraph.addVertex(orientedReadId);
            orientedReadIdVertexMap.insert({orientedReadId, v});
        }
    }



    // Now generate the edges.
    // Each ReadOccurrence generates a pair of reverse complemented edges.
    for(const auto&[readId, occurrences]: readOccurrenceMap) {


        // Now loop over the occurrences of this ReadId.
        for(const auto& occurrence: occurrences) {
            SHASTA2_ASSERT(occurrences.size() > 1);

            // Get the two OrientedReadIds of this ReadId.
            const OrientedReadId orientedReadId0(readId, occurrence.strand);;
            const OrientedReadId orientedReadId1(readId, 1 - occurrence.strand);;

            // Get corresponding BipartiteGraph vertices
            const BipartiteGraph::vertex_descriptor vOrientedRead0 =
                orientedReadIdVertexMap.at(orientedReadId0);
            const BipartiteGraph::vertex_descriptor vOrientedRead1 =
                orientedReadIdVertexMap.at(orientedReadId1);

            // The SegmentPair contains the two BipartiteGraph vertices
            // for this SegmentPair.
            const SegmentPair& segmentPair = segmentPairs[occurrence.segmentPairId];

            const BipartiteGraph::vertex_descriptor vSegment0 = segmentPair.segmentInfos[0].v;
            auto[e, ignore] = boost::add_edge(vOrientedRead0, vSegment0,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);

            const BipartiteGraph::vertex_descriptor vSegment1 = segmentPair.segmentInfos[1].v;
            auto[eRc, ignoreRc] = boost::add_edge(vOrientedRead1, vSegment1,
                BipartiteGraphEdge(occurrence.frequency), bipartiteGraph);

            bipartiteGraph.edgePairs.push_back({e, eRc, occurrence.frequency});

        }
    }

    // Sort the EdgePairs by decreasing frequency.
    std::ranges::sort(
        bipartiteGraph.edgePairs,
        std::greater<uint64_t>(),
        &BipartiteGraph::EdgePair::frequency);
}



// Add a vertex representing a Segment.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(uint64_t segmentPairId, uint64_t segmentIndexInPair)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(segmentPairId, segmentIndexInPair), bipartiteGraph);
}




void StrandContact::writeBipartiteGraphGraphviz(const Split& split)
{
    if(not html) {
        return;
    }

    const string dotFileName = debugOutputBaseName + "-StrandContact-" +
        to_string(strandContactId) + ".dot";
    bipartiteGraph.writeGraphviz(dotFileName, segmentPairs, assemblyGraph, split);

    const double timeout = 30.;
    const string options = "-Nshape=point -Epenwidth=0.2 -Gratio=expand -Gsize=15";
    html << "<h2>Bipartite graph</h2>"
        "<br>In the bipartite graph, each vertex represents a segment or "
        "an oriented read. Oriented reads are displayed as small dots. ";

    html << "<br><br><table>"
        "<tr><th class=left>Total number of vertices<td class=centered>" << num_vertices(bipartiteGraph) <<
        "<tr><th class=left>Number of vertices representing segments<td class=centered>" << 2 * segmentPairs.size() <<
        "<tr><th class=left>Number of vertices representing oriented reads<td class=centered>" << 2 * readOccurrenceMap.size() <<
        "<tr><th class=left>Total number of edges<td class=centered>" << num_edges(bipartiteGraph) <<
        "</table>";

    html << "<br>" << dotFileName << "<br>";

    try {
        graphvizToHtml(dotFileName, "sfdp", timeout, options, html, true);
    } catch (std::exception&) {
        html << "The bipartite graph is too complex to display.";
    }
}



void StrandContact::writeBipartiteGraphCustom(const Split& split)
{
    if(not html) {
        return;
    }

    // If command "customLayout" is not available, don't do anything.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA2_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    if(returnCode != 0) {
        return;
    }

    // Create a map containing the desired length for each edge.
    std::map<BipartiteGraph::edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraph::vertex_descriptor v0 = source(e, bipartiteGraph);
        const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
        const bool isSameComponent = (split.vertexComponent[v0] == split.vertexComponent[v1]);
        const double length = (isSameComponent ? 1. : 10.);
        edgeLengthMap.insert({e, length});
    }

    // Compute the layout.
    std::map<BipartiteGraph::vertex_descriptor, array<double, 2> > positionMap;
    const int quality = 2;
    const double timeout = 30.;
    const auto layoutReturnCode = computeLayoutCustom(bipartiteGraph, edgeLengthMap, positionMap, quality, timeout);
    if(layoutReturnCode != ComputeLayoutReturnCode::Success) {
        html << "<br>The custom layout of the bipartite graph cannot be displayed.";
        return;
    }

    // Compute the bounding box of the layout.
    double xMin = std::numeric_limits<double>::max();
    double xMax = std::numeric_limits<double>::min();
    double yMin = xMin;
    double yMax = xMax;
    for(const auto& p: positionMap) {
        const array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        xMin = min(xMin, x);
        xMax = max(xMax, x);
        yMin = min(yMin, y);
        yMax = max(yMax, y);
    }

    // Enlarge the bounding box a bit.
    const double extend = 0.05 * max(xMax-xMin, yMax-yMin);
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;

    // Make it square,
    if((xMax - xMin) > (yMax - yMin)) {
        const double delta = ((xMax - xMin) - (yMax - yMin)) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = ((yMax - yMin) - (xMax - xMin)) / 2.;;
        xMin -= delta;
        xMax += delta;
    }


    // Begin the svg.
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    const uint64_t sizePixels = 900;
    html << std::scientific;
    const string svgId = "BipartiteGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  sizePixels <<
        "' height='" << sizePixels <<
        "' viewbox='" << xMin << " " << yMin << " " <<
        xMax-xMin << " " <<
        yMax-yMin << "'"
        " style='background-color:#f0f0f0'"
        ">\n";



    // Write the edges first so they don't obscure the vertices.
    const double edgeThicknessFactor = (xMax - xMin) * 0.0001;
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraph::vertex_descriptor v0 = source(e, bipartiteGraph);
        const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
        const auto&[x0, y0] = positionMap.at(v0);
        const auto&[x1, y1] = positionMap.at(v1);
        const bool isSameComponent = (split.vertexComponent[v0] == split.vertexComponent[v1]);
        const string color = (isSameComponent ? "Black" : "Red");

        const uint64_t frequency = bipartiteGraph[e].frequency;
        const double thickness = edgeThicknessFactor * (1. + 6. * std::log10(frequency));

        html <<
            "<line x1='" << x0 << "' y1='" << y0 <<
            "' x2='" << x1 << "' y2='" << y1 <<
            "' stroke='" << color <<
            "' stroke-width='" << thickness <<
            "' />";
    }



    // Write the vertices.
    const double segmentRadius = (xMax - xMin) * 0.006;
    const double orientedReadRadius = (xMax - xMin) * 0.002;
    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const string color = componentColor(split.vertexComponent[v]);
        const double radius = (bipartiteGraph[v].isSegment ? segmentRadius : orientedReadRadius);
        const auto&[x, y] = positionMap.at(v);
        html << "<circle cx='" << x << "' cy='" << y <<
            "' fill='" << color <<
            "' r='" << radius <<
            "' />";
    }


    // Finish the svg.
    html << "</svg></div>";
}



void BipartiteGraph::writeGraphviz(
    const string& fileName,
    const vector<SegmentPair>& segmentPairs,
    const AssemblyGraph& assemblyGraph,
    const Split& split) const
{
    const BipartiteGraph& bipartiteGraph = *this;

    ofstream dot(fileName);

    dot << "graph BipartiteGraph {\n";



    // Vertices.
    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex = bipartiteGraph[v];
        const string color = StrandContact::componentColor(split.vertexComponent[v]);

        if(vertex.isSegment) {
            const uint64_t segmentPairId = vertex.segmentPairId;
            const uint64_t segmentIndexInPair = vertex.segmentIndexInPair;
            const SegmentPair& segmentPair = segmentPairs[segmentPairId];
            const SegmentInfo& segmentInfo = segmentPair.segmentInfos[segmentIndexInPair];
            const Segment segment = segmentInfo.segment;
            dot << assemblyGraph.id(segment);
            dot << " [width=0.1";
        } else {
            dot << "\"" << vertex.orientedReadId << "\"";
            dot << " [width=0.02";
        }
        dot << " color=\"" << color << "\"";
        dot << "]";
        dot << ";\n";
    }



    // Edges.
    for(uint64_t edgePairIndex=0; edgePairIndex<edgePairs.size(); edgePairIndex++) {
        const EdgePair& edgePair = edgePairs[edgePairIndex];
        const array<edge_descriptor, 2> edgePairEdges = {edgePair.e, edgePair.eRc};

        for(const edge_descriptor e: edgePairEdges) {

            const vertex_descriptor v0 = source(e, bipartiteGraph);
            const vertex_descriptor v1 = target(e, bipartiteGraph);
            const BipartiteGraphVertex& vertex0 = bipartiteGraph[v0];
            const BipartiteGraphVertex& vertex1 = bipartiteGraph[v1];

            const uint64_t frequency = bipartiteGraph[e].frequency;
            const double thickness = 0.1 * (1. + std::log10(frequency));

            if(vertex0.isSegment) {
                const uint64_t segmentPairId0 = vertex0.segmentPairId;
                const uint64_t segmentIndexInPair0 = vertex0.segmentIndexInPair;
                const SegmentPair& segmentPair0 = segmentPairs[segmentPairId0];
                const SegmentInfo& segmentInfo0 = segmentPair0.segmentInfos[segmentIndexInPair0];
                const Segment segment0 = segmentInfo0.segment;
                dot << assemblyGraph.id(segment0);
            } else {
                dot << "\"" << vertex0.orientedReadId << "\"";
            }

            dot << "--";

            if(vertex1.isSegment) {
                const uint64_t segmentPairId1 = vertex1.segmentPairId;
                const uint64_t segmentIndexInPair1 = vertex1.segmentIndexInPair;
                const SegmentPair& segmentPair1 = segmentPairs[segmentPairId1];
                const SegmentInfo& segmentInfo1 = segmentPair1.segmentInfos[segmentIndexInPair1];
                const Segment segment1 = segmentInfo1.segment;
                dot << assemblyGraph.id(segment1);
            } else {
                dot << "\"" << vertex1.orientedReadId << "\"";
            }

            dot << "[";
            dot << "penwidth=\"" << thickness << "\"";

            if(split.isCrossStrandEdgePair(edgePairIndex)) {
                dot << " color=red";
            }

            dot << "]";

            dot << ";\n";
        }

    }

    dot << "}\n";
}



// Add a vertex representing an OrientedReadId.
BipartiteGraph::vertex_descriptor BipartiteGraph::addVertex(OrientedReadId orientedReadId)
{
    BipartiteGraph& bipartiteGraph = *this;
    return boost::add_vertex(BipartiteGraphVertex(orientedReadId), bipartiteGraph);

}



// At the first iteration, use the EdgePairs sorted by decreasing
// frequency. At subsequent iterations, use random shuffles of
// the EdgePairs, weighted by their frequency.
void StrandContact::computeSplit(Split& bestSplit) const
{
    // EXPOSE WHEN CODE STABILIZES.
    const uint64_t iterationCount = 1000;

    // Gather the weight of each EdgePair.
    vector<double> weights;
    for(const auto& edgePair: bipartiteGraph.edgePairs) {
        weights.push_back(double(edgePair.frequency));
    }

    std::mt19937 generator;
    vector<uint64_t> shuffle;

    vector<uint64_t> edgePairsIndexes(bipartiteGraph.edgePairs.size());

    if(html) {
        html << "<h2>Strand separation iterations</h2>"
            "<table><tr><th>Iteration<th>Cross-strand<br>frequency";
    }

    Split split;
    for(uint64_t iteration=0; iteration<iterationCount; iteration++) {
        split.clear();
        if(iteration == 0) {
            std::ranges::iota(edgePairsIndexes, 0);
        } else {
            weightedShuffle(weights, generator, edgePairsIndexes);
        }
        bipartiteGraph.computeSplit(edgePairsIndexes, split);

        const bool isBestSoFar = (iteration == 0) or (split.crossStrandFrequency < bestSplit.crossStrandFrequency);
        if(isBestSoFar) {
            bestSplit = split;
        }

        if(html) {
            html << "<tr";
            if(isBestSoFar) {
                html << " style='background-color:Pink'";
            }
            html << ">";
            html <<
                "<td class=centered>" << iteration <<
                "<td class=centered>" << 2 * split.crossStrandFrequency;
        }
    }

    if(html) {
        html << "</table>";
    }
}



// Return the reverse complement of a vertex.
// Because vertices are added in reverse complemented pairs,
// pairs of reverse complemented vertices have consecutive vertex_descriptors.
BipartiteGraph::vertex_descriptor BipartiteGraph::reverseComplement(vertex_descriptor v) const
{
    return v ^ 1;
}



// This processes the EdgePairs in the order described by
// the edgePairsIndexes.
// It stores the components in the components vector
// and also fills in the componentIndex in all the vertices.
// Reverse complemented components are numbered consecutively.
void BipartiteGraph::computeSplit(
    const vector<uint64_t>& edgePairsIndexes,
    Split& split) const
{
    const BipartiteGraph& bipartiteGraph = *this;
    split.clear();

    DisjointSets disjointSets(num_vertices(bipartiteGraph));

    for(const uint64_t edgePairIndex: edgePairsIndexes) {
        const auto& edgePair = edgePairs[edgePairIndex];
        const auto eA = edgePair.e;
        const auto eB = edgePair.eRc;

        const auto v0A = source(eA, bipartiteGraph);
        const auto v1A = target(eA, bipartiteGraph);
        const auto v0B = source(eB, bipartiteGraph);
        const auto v1B = target(eB, bipartiteGraph);

        const auto v0ARc = reverseComplement(v0A);
        const auto v1ARc = reverseComplement(v1A);
        const auto v0BRc = reverseComplement(v0B);
        const auto v1BRc = reverseComplement(v1B);

        const bool strandViolationA = (disjointSets.findSet(v1A) == disjointSets.findSet(v0ARc));
        const bool strandViolationB = (disjointSets.findSet(v1B) == disjointSets.findSet(v0BRc));
        const bool strandViolationARc = (disjointSets.findSet(v0A) == disjointSets.findSet(v1ARc));
        const bool strandViolationBRc = (disjointSets.findSet(v0B) == disjointSets.findSet(v1BRc));

        const bool strandViolation = strandViolationA;
        SHASTA2_ASSERT(strandViolationB == strandViolation);
        SHASTA2_ASSERT(strandViolationARc == strandViolation);
        SHASTA2_ASSERT(strandViolationBRc == strandViolation);

        if(strandViolation) {
            split.crossStrandEdgePairIndexes.push_back(edgePairIndex);
            split.crossStrandFrequency += edgePair.frequency;
        } else {
            disjointSets.unionSet(v0A, v1A);
            disjointSets.unionSet(v0B, v1B);
        }
    }

    std::ranges::sort(split.crossStrandEdgePairIndexes);

    // Gather the components.
    disjointSets.gatherComponents(1, split.components);

    // Reorder the components so pairs of
    // reverse complemented components are numbered consecutively.
    // Because of the way vertices are created, this can be done simply by
    // sorting by the first index of each component.
    class SortByFirstElement {
    public:
    public:
         bool operator()(const vector<vertex_descriptor>& v0, const vector<vertex_descriptor>& v1) const
        {
             SHASTA2_ASSERT(not v0.empty());
             SHASTA2_ASSERT(not v1.empty());
             return v0.front() < v1.front();
        }
    };
    std::ranges::sort(split.components, SortByFirstElement());

    // Store the componentId of the vertices.
    split.vertexComponent.resize(num_vertices(bipartiteGraph));
    for(uint64_t componentId=0; componentId<split.components.size(); componentId++) {
        const vector<vertex_descriptor>& component = split.components[componentId];
        for(const vertex_descriptor v: component) {
            split.vertexComponent[v] = componentId;
        }
    }
}



void StrandContact::writeSplitSummary(const Split& split) const
{
    if(not html) {
        return;
    }

    uint64_t totalEdgeCount = 0;
    uint64_t totalEdgeFrequency = 0;
    uint64_t crossStrandEdgeCount = 0;
    uint64_t crossStrandEdgeFrequency = 0;
    BGL_FORALL_EDGES(e, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraph::vertex_descriptor v0 = source(e, bipartiteGraph);
        const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);

        const uint64_t frequency = bipartiteGraph[e].frequency;
        ++totalEdgeCount;
        totalEdgeFrequency += frequency;

        if(split.vertexComponent[v0] != split.vertexComponent[v1]) {
            ++crossStrandEdgeCount;
            crossStrandEdgeFrequency += frequency;
        }
    }

    html <<
        std::setprecision(6) <<
        "<table>"
        "<tr><th><th>Total<th>Cross-strand<th>Ratio"
        "<tr><th class=left>Number of edges"
        "<td class=centered>" << totalEdgeCount <<
        "<td class=centered>" << crossStrandEdgeCount <<
        "<td class=centered>" << double(crossStrandEdgeCount) / double(totalEdgeCount) <<
        "<tr><th class=left>Edge frequency"
        "<td class=centered>" << totalEdgeFrequency <<
        "<td class=centered>" << crossStrandEdgeFrequency <<
        "<td class=centered>" << double(crossStrandEdgeFrequency) / double(totalEdgeFrequency) <<
        "</table>";
}



BipartiteGraphVertexStatistics BipartiteGraph::getVertexStatistics(
    const Split& split,
    vertex_descriptor v0) const
{
    const BipartiteGraph& bipartiteGraph = *this;
    BipartiteGraphVertexStatistics statistics;

    BGL_FORALL_OUTEDGES(v0, e, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraph::vertex_descriptor v1 = target(e, bipartiteGraph);
        const uint64_t frequency = bipartiteGraph[e].frequency;
        ++statistics.totalEdgeCount;
        statistics.totalEdgeFrequency += frequency;
        if(split.vertexComponent[v0] != split.vertexComponent[v1]) {
            ++statistics.crossStrandEdgeCount;
            statistics.crossStrandEdgeFrequency += frequency;
        }
    }

    return statistics;
}



double BipartiteGraphVertexStatistics::crossStrandEdgeRatio() const
{
    return double(crossStrandEdgeCount) / double(totalEdgeCount);
}



double BipartiteGraphVertexStatistics::crossStrandEdgeFrequencyRatio() const
{
    return double(crossStrandEdgeFrequency) / double(totalEdgeFrequency);
}



uint64_t BipartiteGraphVertexStatistics::sameStrandEdgeCount() const
{
    return totalEdgeCount - crossStrandEdgeCount;
}



uint64_t BipartiteGraphVertexStatistics::sameStrandEdgeFrequency() const
{
    return totalEdgeFrequency - crossStrandEdgeFrequency;
}



int64_t BipartiteGraphVertexStatistics::crossStrandEdgeFrequencyExcess() const
{
    return int64_t(crossStrandEdgeFrequency) - int64_t(sameStrandEdgeFrequency());
}



void StrandContact::writeSplitDetails(const Split& split) const
{
    if(not html) {
        return;
    }

    const string csvFileName =
        debugOutputBaseName + "-StrandContact-" + to_string(strandContactId) + "-EvaluateStrandSeparation.csv";
    ofstream csv(csvFileName);

    csv <<
        "Vertex,Segment,OrientedReadId,"
        "Total edge count,Cross-strand edge count,Cross-strand edge ratio,"
        "Total edge frequency,Same-strand edge frequency,Cross-strand edge frequency,"
        "Cross-strand edge frequency excess,"
        "Cross-strand edge frequency ratio,\n";

    BGL_FORALL_VERTICES(v, bipartiteGraph, BipartiteGraph) {
        const BipartiteGraphVertex& vertex = bipartiteGraph[v];
        const BipartiteGraphVertexStatistics statistics = bipartiteGraph.getVertexStatistics(split, v);

        csv << v << ",";

        if(vertex.isSegment) {
            csv << segmentPairs[vertex.segmentPairId].segmentInfos[vertex.segmentIndexInPair].id;
        }
        csv << ",";

        if(not vertex.isSegment) {
            csv << vertex.orientedReadId;
        }
        csv << ",";

        csv << statistics.totalEdgeCount << ",";

        if(statistics.crossStrandEdgeCount) {
            csv << statistics.crossStrandEdgeCount;
        }
        csv << ",";

        if(statistics.crossStrandEdgeCount) {
            csv << statistics.crossStrandEdgeRatio();
        }
        csv << ",";

        csv << statistics.totalEdgeFrequency << ",";
        csv << statistics.sameStrandEdgeFrequency() << ",";

        if(statistics.crossStrandEdgeFrequency) {
            csv << statistics.crossStrandEdgeFrequency;
        }
        csv << ",";

        csv << statistics.crossStrandEdgeFrequencyExcess() << ",";


        if(statistics.crossStrandEdgeFrequency) {
            csv << statistics.crossStrandEdgeFrequencyRatio();
        }
        csv << ",";

        csv << "\n";

    }
    html << "<br>See <a href='" << csvFileName << "'>" << csvFileName <<
        "</a> for detailed evaluation of strand separation.";


}



void Split::clear()
{
    components.clear();
    vertexComponent.clear();
    crossStrandEdgePairIndexes.clear();
    crossStrandFrequency = 0;
}



bool Split::isCrossStrandEdgePair(uint64_t i) const
{
    return std::ranges::binary_search(crossStrandEdgePairIndexes, i);
}



string StrandContact::componentColor(uint64_t componentId)
{
    return randomHslColor(componentId + 1000, 0.5, 0.6);
}



string StrandContact::ambiguousColor()
{
    return hslToRgbString(0.0, 0., 0.4);
}



// We update the AssemblyGraph as follows:
// - We make identical copies of Segments that are not entrances or exits
//   and that either:
//   * Are classified as ambiguous
//   or
//   * Are assigned an evenly number component.
//   For the ambiguous case, it may be necessary to remove from the copy
//   oriented reads that don't belong to the same evenly number component.
// - Except for interface vertices (targets of entrances and source of exits),
//   the copies of the Segments use newly generated vertices.
//   So they are only connected to each other and to the entrances and exits.
// - As we create new vertices and Segments, we also
//   create the corresponding reverse complemented copies.
void StrandContact::updateAssemblyGraph()
{
    using vertex_descriptor = AssemblyGraph::vertex_descriptor;

    // Find the interface vertices.
    vector<vertex_descriptor> interfaceVertices;
    for(const SegmentPair& segmentPair: segmentPairs) {
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            if(segmentInfo.isEntrance) {
                const Segment segment = segmentInfo.segment;
                const vertex_descriptor v = target(segment, assemblyGraph);
                interfaceVertices.push_back(v);
            }
            if(segmentInfo.isExit) {
                const Segment segment = segmentInfo.segment;
                const vertex_descriptor v = source(segment, assemblyGraph);
                interfaceVertices.push_back(v);
            }
        }
    }
    deduplicate(interfaceVertices);

    // A map that gives the new vertex corresponding to an old vertex
    // in the replicated copies of segments.
    // For interface vertices, the new vertex is the same as the old vertex.
    std::map<vertex_descriptor, vertex_descriptor> newVertexMap;
    for(const vertex_descriptor v: interfaceVertices) {
        newVertexMap.insert({v, v});
    }



    // Replicate Segments as described at the beginning of this function.
    std::map<Segment, Segment> newSegmentMap;
    for(const SegmentPair& segmentPair: segmentPairs) {
        const bool isAmbiguous = (segmentPair.crossStrandEdgeFrequencyRatio > maxCrossStrandFrequencyRatio);
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            if(segmentInfo.isEntrance) {
                continue;
            }
            if(segmentInfo.isExit) {
                continue;
            }
            const bool isEvenComponent = (segmentInfo.componentId & 1) == 0;
            if(isAmbiguous or isEvenComponent) {

                // Ok, we are going to make a copy of this Segment.
                const Segment segment = segmentInfo.segment;

                // But first we need to make sure we have the vertices.
                // Get the new vertices for the copy of this segment, creating them if necessary.
                // Also create their reverse complements.
                const vertex_descriptor v0Old = source(segment, assemblyGraph);
                const vertex_descriptor v1Old = target(segment, assemblyGraph);
                const auto it0 = newVertexMap.find(v0Old);
                const auto it1 = newVertexMap.find(v1Old);
                vertex_descriptor v0New;
                vertex_descriptor v1New;
                if(it0 == newVertexMap.end()) {
                    const AnchorId anchorId0 = assemblyGraph[v0Old].anchorId;
                    v0New = add_vertex(AssemblyGraphVertex(anchorId0, assemblyGraph.nextVertexId++), assemblyGraph);
                    newVertexMap.insert({v0Old, v0New});
                    assemblyGraph.createReverseComplementVertex(v0New);
                } else {
                    v0New = it0->second;
                }
                if(it1 == newVertexMap.end()) {
                    const AnchorId anchorId1 = assemblyGraph[v1Old].anchorId;
                    v1New = add_vertex(AssemblyGraphVertex(anchorId1, assemblyGraph.nextVertexId++), assemblyGraph);
                    newVertexMap.insert({v1Old, v1New});
                    assemblyGraph.createReverseComplementVertex(v1New);
                } else {
                    v1New = it1->second;
                }

                // Make an exact copy of this Segment, between these new vertices.
                // For the ambiguous case, it may be necessary to remove from the copy
                // oriented reads that don't belong to the same evenly number component.
                auto[newSegment, wasAdded] = add_edge(v0New, v1New, assemblyGraph[segment], assemblyGraph);
                SHASTA2_ASSERT(wasAdded);
                newSegmentMap.insert({segment, newSegment});

                // Also create the reverse complement edge.
                assemblyGraph.createReverseComplementEdge(newSegment);
            }
        }
    }



    // Recursively prune new segments that are dangling and that are a copy of an ambiguous segment.
    while(true) {
        vector<Segment> originalCopiesOfSegmentsToBeRemoved;
        for(const SegmentPair& segmentPair: segmentPairs) {
            const bool isAmbiguous = (segmentPair.crossStrandEdgeFrequencyRatio > maxCrossStrandFrequencyRatio);
            for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
                if(segmentInfo.isEntrance) {
                    continue;
                }
                if(segmentInfo.isExit) {
                    continue;
                }
                if(not isAmbiguous) {
                    continue;
                }
                const Segment segment = segmentInfo.segment;
                const auto it = newSegmentMap.find(segment);
                if(it == newSegmentMap.end()) {
                    // We already removed it.
                    continue;
                }
                const Segment newSegment = it->second;
                const vertex_descriptor v0 = source(newSegment, assemblyGraph);
                const vertex_descriptor v1 = target(newSegment, assemblyGraph);
                const bool isDangling = ((in_degree(v0, assemblyGraph) == 0) or(out_degree(v1, assemblyGraph) == 0));
                if(isDangling) {
                    originalCopiesOfSegmentsToBeRemoved.push_back(segment);
                }
            }
        }

        if(originalCopiesOfSegmentsToBeRemoved.empty()) {
            break;
        }

        for(const Segment segment: originalCopiesOfSegmentsToBeRemoved) {
            const Segment newSegment = newSegmentMap.at(segment);
            const Segment newSegmentRc = assemblyGraph[newSegment].eRc;
            if(html) {
                html << "<br>Removing dangling ambiguous segment " << id(newSegment) << flush;
            }
            boost::remove_edge(newSegment, assemblyGraph);
            boost::remove_edge(newSegmentRc, assemblyGraph);
            newSegmentMap.erase(segment);
        }
    }



    // Now we can remove all the Segments of this StrandContact,
    // except for entrances or exits.
    for(const SegmentPair& segmentPair: segmentPairs) {
        for(const SegmentInfo& segmentInfo: segmentPair.segmentInfos) {
            if(segmentInfo.isEntrance) {
                continue;
            }
            if(segmentInfo.isExit) {
                continue;
            }
            boost::remove_edge(segmentInfo.segment, assemblyGraph);
        }
    }

    // This leaves some isolated vertices that will be removed later.
}

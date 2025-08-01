// Shasta.
#include "LocalAnchorGraph.hpp"
#include "AnchorGraph.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "computeLayout.hpp"
#include "html.hpp"
#include "HttpServer.hpp"
#include "MurmurHash2.hpp"
#include "tmpDirectory.hpp"
#include "runCommandWithTimeout.hpp"
using namespace shasta;

// Boost libraries.
#include <boost/graph/iteration_macros.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

// Standard library.
#include "fstream.hpp"
#include <queue>



// New constructor from the global AnchorGraph.
// vertex and edge descriptors use a L (local) suffix for the LocalAnchorGraph
// and G (global) for the AnchorGraph.
LocalAnchorGraph::LocalAnchorGraph(
    const Anchors& anchors,
    const AnchorGraph& anchorGraph,
    const vector<AnchorId>& anchorIds,
    uint64_t maxDistance,
	uint64_t minCoverage,
    bool edgesMarkedForAssembly) :
    anchors(anchors),
    anchorGraphPointer(&anchorGraph),
    maxDistance(maxDistance)
{
    LocalAnchorGraph& localAnchorGraph = *this;

    // Initialize a BFS from these AnchorIds.
    std::queue<vertex_descriptor> q;
    for(const AnchorId anchorId: anchorIds) {
        SHASTA_ASSERT(not vertexMap.contains(anchorId));
        const vertex_descriptor vL = boost::add_vertex(LocalAnchorGraphVertex(anchorId, 0), localAnchorGraph);
        vertexMap.insert({anchorId, vL});
        q.push(vL);
    }



    // BFS to find the vertices. We will add the edges later.
    while(not q.empty()) {

        // Dequeue a vertex.
        vertex_descriptor v0L = q.front();
        q.pop();
        const LocalAnchorGraphVertex& vertex0L = localAnchorGraph[v0L];
        const AnchorId anchorId0 = vertex0L.anchorId;

        const uint64_t distance0 = vertex0L.distance;
        const uint64_t distance1 = distance0 + 1;

        // Get the corresponding global AnchorGraph vertex.
        // In the global AnchorGraph, vertex descriptors are AnchorIds.
        const AnchorGraph::vertex_descriptor v0G = anchorId0;

        // Loop over out-edges.
        BGL_FORALL_OUTEDGES(v0G, eG, anchorGraph, AnchorGraph) {
            if(edgesMarkedForAssembly and (not anchorGraph[eG].useForAssembly)) {
                continue;
            }
            if(anchorGraph[eG].anchorPair.size() < minCoverage) {
            	continue;
            }
            const AnchorGraph::vertex_descriptor v1G = target(eG, anchorGraph);
            const AnchorId anchorId1 = v1G;

            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            const vertex_descriptor v1L = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), localAnchorGraph);
            vertexMap.insert({anchorId1, v1L});
            if(distance1 < maxDistance) {
                q.push(v1L);
            }
        }

        // Loop over in-edges.
        BGL_FORALL_INEDGES(v0G, eG, anchorGraph, AnchorGraph) {
            if(edgesMarkedForAssembly and (not anchorGraph[eG].useForAssembly)) {
                continue;
            }
            if(anchorGraph[eG].anchorPair.size() < minCoverage) {
            	continue;
            }
            const AnchorGraph::vertex_descriptor v1G = source(eG, anchorGraph);
            const AnchorId anchorId1 = v1G;

            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                continue;
            }

            const vertex_descriptor v1L = boost::add_vertex(LocalAnchorGraphVertex(anchorId1, distance1), localAnchorGraph);
            vertexMap.insert({anchorId1, v1L});
            if(distance1 < maxDistance) {
                q.push(v1L);
            }
        }
    }



    // Now add the edges.
    BGL_FORALL_VERTICES(v0L, localAnchorGraph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex0L = localAnchorGraph[v0L];
        const AnchorId anchorId0 = vertex0L.anchorId;

        // Get the corresponding global AnchorGraph vertex.
        // In the global AnchorGraph, vertex descriptors are AnchorIds.
        const AnchorGraph::vertex_descriptor v0G = anchorId0;

        // Loop over out-edges.
        BGL_FORALL_OUTEDGES(v0G, eG, anchorGraph, AnchorGraph) {
            if(edgesMarkedForAssembly and (not anchorGraph[eG].useForAssembly)) {
                continue;
            }
            if(anchorGraph[eG].anchorPair.size() < minCoverage) {
            	continue;
            }
            const AnchorGraph::vertex_descriptor v1G = target(eG, anchorGraph);
            const AnchorId anchorId1 = v1G;

            auto it1 = vertexMap.find(anchorId1);
            if(it1 != vertexMap.end()) {
                const vertex_descriptor v1L = it1->second;

                LocalAnchorGraphEdge edge;
                edge.eG = eG;
                add_edge(v0L, v1L, edge, localAnchorGraph);
            }

        }
    }
}



void LocalAnchorGraph::writeGraphviz(
    const string& fileName,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer) const
{
    ofstream file(fileName);
    writeGraphviz(file, options, assemblyGraph3Pointer);
}



void LocalAnchorGraph::writeGraphviz(
    ostream& s,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer) const
{
    const LocalAnchorGraph& graph = *this;
    SHASTA_ASSERT(anchorGraphPointer);
    const AnchorGraph& anchorGraph = *anchorGraphPointer;

    if(options.vertexColoring == "byAssemblyAnnotations") {
        SHASTA_ASSERT(assemblyGraph3Pointer);
    }

    AnchorId referenceAnchorId = invalid<AnchorId>;
    if(options.vertexColoring == "byReadComposition") {
        referenceAnchorId = anchorIdFromString(options.referenceAnchorIdString);
        if((referenceAnchorId == invalid<AnchorId>) or (referenceAnchorId >= anchors.size())) {
            throw runtime_error("Invalid reference anchor id " + options.referenceAnchorIdString +
                ". Must be a number between 0 and " +
                to_string(anchors.size() / 2 - 1) + " followed by + or -.");
        }
    }
    const uint64_t referenceAnchorIdCoverage = anchors[referenceAnchorId].coverage();

    s << "digraph LocalAnchorGraph {\n";

    // Write the vertices.
    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex = graph[v];
        const AnchorId anchorId = vertex.anchorId;
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = anchors[anchorId].coverage();

        // Get annotation information, if needed.
        bool hasVertexAnnotation = false;
        vector<AssemblyGraph::edge_descriptor> annotationEdges;
        if(options.vertexColoring == "byAssemblyAnnotations") {
            hasVertexAnnotation = assemblyGraph3Pointer->hasVertexAnnotation(anchorId);
            if(not hasVertexAnnotation) {
                assemblyGraph3Pointer->findAnnotationEdges(anchorId, annotationEdges);
            }
        }

        // Annotation text.
        string annotationText;
        if(options.vertexColoring == "byAssemblyAnnotations") {
            if(hasVertexAnnotation) {
                annotationText = "Vertex";
            } else if(annotationEdges.size() > 1) {
                annotationText = "Multiple";
            } else if(annotationEdges.size() == 1) {
                const uint64_t segmentId = (*assemblyGraph3Pointer)[annotationEdges.front()].id;
                annotationText = to_string(segmentId);
            }
        }

        // Vertex name.
        s << "\"" << anchorIdString << "\"";

        // Begin vertex attributes.
        s << "[";

        // URL
        s << "URL=\"exploreAnchor?anchorIdString=" << HttpServer::urlEncode(anchorIdString) << "\"";

        // Tooltip.
        s << " tooltip=\"" << anchorIdString << " " << coverage;
        if(not annotationText.empty()) {
            s << " " << annotationText;
        }
        s << "\"";

        // Label.
        if(options.vertexLabels) {
            s << " label=\"" << anchorIdString << "\\n" << coverage;
            if(not annotationText.empty()) {
                s << "\\n" << annotationText;
            }
            s << "\"";
        }



        // Color.
        if(vertex.distance == 0) {
            s << " fillcolor=LightBlue style=filled";
        } else if(vertex.distance == maxDistance) {
            s << " fillcolor=Cyan style=filled";
        } else {

            // Color by similarity of read composition with the reference Anchor.
            if(options.vertexColoring == "byReadComposition") {
                AnchorPairInfo info;
                anchors.analyzeAnchorPair(referenceAnchorId, anchorId, info);

                double hue = 1.;    // 0=red, 1=green.
                if(options.similarityMeasure == "commonCount") {
                    // By common count.
                    hue = double(info.common) / double(referenceAnchorIdCoverage);

                } else if(options.similarityMeasure == "jaccard") {
                    // By Jaccard similarity.
                    hue = info.jaccard();
                } else {
                    // By corrected Jaccard similarity.
                    hue = info.correctedJaccard();
                 }

                const string colorString = "\"" + to_string(hue / 3.) + " 1 1\"";
                if(options.vertexLabels) {
                    s << " style=filled fillcolor=" << colorString;
                } else {
                    s << " color=" << colorString;
                    s << " fillcolor=" << colorString;
                }
            } else if(options.vertexColoring == "byAssemblyAnnotations") {
                if(hasVertexAnnotation) {
                    s << " style=filled fillcolor=Red";
                } else {
                    if(annotationEdges.size() > 1) {
                        s << " style=filled fillcolor=Green";    // Multiple segments
                    } else if(annotationEdges.size() == 1) {
                        const uint64_t segmentId = (*assemblyGraph3Pointer)[annotationEdges.front()].id;
                        const uint32_t hashValue = MurmurHash2(&segmentId, sizeof(segmentId), 759);
                        const uint32_t hue = hashValue % 1000;
                        s << "style=filled fillcolor=\"" << double(hue / 1000.) << " .6 .9\"";
                    }
                }
            }        }



        // Size.
        if(not options.vertexLabels) {
            const double displaySize =
                (options.vertexSizeByCoverage ?
                options.vertexSize * sqrt(0.1 * double(coverage)) :
                options.vertexSize
                ) / 72.;
            s << " width=" << displaySize ;
            s << " penwidth=" << 0.5 * displaySize;
        }

        // End vertex attributes.
        s << "]";

        // End the line for this vertex.
        s << ";\n";
    }



    // Write the edges.
    BGL_FORALL_EDGES(e, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edge = graph[e];

        const vertex_descriptor v0 = source(e, graph);
        const vertex_descriptor v1 = target(e, graph);

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const LocalAnchorGraphVertex& vertex1 = graph[v1];

        const AnchorId anchorId0 = vertex0.anchorId;
        const AnchorId anchorId1 = vertex1.anchorId;

        const string anchorId0String = anchorIdToString(anchorId0);
        const string anchorId1String = anchorIdToString(anchorId1);

        const AnchorGraph::edge_descriptor eG = edge.eG;
        const AnchorGraphEdge& edgeG = anchorGraph[eG];
        const AnchorPair& anchorPair = edgeG.anchorPair;
        SHASTA_ASSERT(anchorPair.anchorIdA == anchorId0);
        SHASTA_ASSERT(anchorPair.anchorIdB == anchorId1);

        const uint64_t coverage = anchorPair.orientedReadIds.size();
        const uint64_t offset = edgeG.offset;

        string color = "Black";

        if(options.edgeColoring == "byFlags") {
			if(edgeG.isParallelEdge) {
				color = "Orange";
			} else if(edgeG.addedAtDeadEnd) {
				color = "Magenta";
			} else if(not edgeG.useForAssembly) {
				color = "GreenYellow";
			}
        }

        else if(options.edgeColoring == "random") {
            // To decide the color, hash the AnchorIds.
            // This way we always get the same color for the same edge.
            const auto p = make_pair(anchorId0, anchorId1);
            const uint32_t hashValue = MurmurHash2(&p, sizeof(p), 759);
            const uint32_t hue = hashValue % 360;
            color = to_string(hue / 360.) + " 1. 1.";
        }


        s << "\"" << anchorId0String << "\"->";
        s << "\"" << anchorId1String << "\"";

        // Begin edge attributes.
        s << " [";

        // URL
        s << "URL=\"" + anchorPair.url() + "\"";

        // Tooltip.
        s << " tooltip="
            "\"" << anchorId0String << " to "
            << anchorId1String <<
            ", coverage " << coverage <<
            ", offset " << offset << "\"";

        // Label.
        if(options.edgeLabels) {
            s << " label=\"" <<
                coverage <<
                "\\nOffset " << offset << "\"";
        }

        // Color.
        s << " color=\"" << color << "\"";

        // Thickness.
        s << " penwidth=" << 0.5 * options.edgeThickness * double(coverage);

        // Arrow size.
        s << " arrowsize=" << 0.5 * options.arrowSize;

        // Length. Only use by fdp and neato layouts.
        const double displayLength =
            (options.minimumEdgeLength +
                options.additionalEdgeLengthPerKb * 0.001 * double(offset)) / 72.;
        s << " len=" << displayLength;



        // End edge attributes.
        s << "]";

        // End the line for this edge.
        s << ";\n";
    }


    s << "}\n";
}




LocalAnchorGraphDisplayOptions::LocalAnchorGraphDisplayOptions(const vector<string>& request)
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);


    sizePixels = 600;
    HttpServer::getParameterValue(request, "sizePixels", sizePixels);

    layoutMethod = (customLayoutIsAvailable ? "custom" : "sfdp");
    HttpServer::getParameterValue(request, "layoutMethod", layoutMethod);

    vertexColoring = "black";
    HttpServer::getParameterValue(request, "vertexColoring", vertexColoring);

    similarityMeasure = "commonCount";
    HttpServer::getParameterValue(request, "similarityMeasure", similarityMeasure);

    referenceAnchorIdString = "";
    HttpServer::getParameterValue(request, "referenceAnchorId", referenceAnchorIdString);

    assemblyStage = "";
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    edgeColoring = "random";
    HttpServer::getParameterValue(request, "edgeColoring", edgeColoring);

    vertexSize =  1.;
    HttpServer::getParameterValue(request, "vertexSize", vertexSize);

    string vertexSizeByCoverageString;
    vertexSizeByCoverage = HttpServer::getParameterValue(request,
        "vertexSizeByCoverage", vertexSizeByCoverageString);

    string vertexLabelsString;
    vertexLabels = HttpServer::getParameterValue(request,
        "vertexLabels", vertexLabelsString);

    minimumEdgeLength = 1.;
    HttpServer::getParameterValue(request, "minimumEdgeLength", minimumEdgeLength);

    additionalEdgeLengthPerKb = 1.;
    HttpServer::getParameterValue(request, "additionalEdgeLengthPerKb", additionalEdgeLengthPerKb);

    edgeThickness = 1.;
    HttpServer::getParameterValue(request, "edgeThickness", edgeThickness);

    arrowSize = 1.;
    HttpServer::getParameterValue(request, "arrowSize", arrowSize);

    string edgeLabelsString;
    edgeLabels = HttpServer::getParameterValue(request,
        "edgeLabels", edgeLabelsString);
}



void LocalAnchorGraphDisplayOptions::writeForm(ostream& html) const
{
    // Figure out if command "customLayout" is available.
    const int commandStatus = std::system("which customLayout > /dev/null");
    SHASTA_ASSERT(WIFEXITED(commandStatus));
    const int returnCode = WEXITSTATUS(commandStatus);
    const bool customLayoutIsAvailable = (returnCode == 0);

    html <<
        "<tr>"
        "<th class=left title='Graphics size in pixels. "
        "Changing this works better than zooming. Make it larger if the graph is too crowded."
        " Ok to make it much larger than screen size.'>"
        "Graphics size in pixels"
        "<td class=centered><input type=text required name=sizePixels size=8 style='text-align:center'" <<
        " value='" << sizePixels <<
        "'>";

    html <<
        "<tr>"
        "<th>Layout method"
        "<td class=left>"
        "<input type=radio required name=layoutMethod value='sfdp'" <<
        (layoutMethod == "sfdp" ? " checked=on" : "") <<
        ">sfdp"
        "<br><input type=radio required name=layoutMethod value='fdp'" <<
        (layoutMethod == "fdp" ? " checked=on" : "") <<
        ">fdp"
        "<br><input type=radio required name=layoutMethod value='neato'" <<
        (layoutMethod == "neato" ? " checked=on" : "") <<
        ">neato"
        "<br><input type=radio required name=layoutMethod value='dot'" <<
        (layoutMethod == "dot" ? " checked=on" : "") <<
        ">dot";

    // If command "customLayout" is available, add an option for that.
    if(customLayoutIsAvailable) {
        html <<
            "<br><input type=radio required name=layoutMethod value='custom'" <<
            (layoutMethod == "custom" ? " checked=on" : "") <<
            ">custom";
    }

    html <<
        "<tr>"
        "<th>Vertices"
        "<td class=left>"
        "<input type=text name=vertexSize style='text-align:center' required size=6 value=" <<
        vertexSize << "> Vertex size (arbitrary units)"
        "<br><input type=checkbox name=vertexSizeByCoverage" <<
        (vertexSizeByCoverage ? " checked" : "") <<
        "> Size proportional to coverage"

        "<hr>"
        "<input type=checkbox name=vertexLabels" <<
        (vertexLabels ? " checked" : "") <<
        "> Labels (dot layout only)"

        "<hr>"
        "<b>Vertex coloring</b>"

        "<br><input type=radio required name=vertexColoring value='black'" <<
        (vertexColoring == "black" ? " checked=on" : "") << ">Black"
        "<br><input type=radio required name=vertexColoring value='byReadComposition'" <<
        (vertexColoring == "byReadComposition" ? " checked=on" : "") <<
        "> By similarity of read composition using similarity measure:"

        "<div style='padding-left:50px'>"
        "<input type=radio required name=similarityMeasure value='commonCount'" <<
        (similarityMeasure == "commonCount" ? " checked=on" : "") << ">Number of common oriented reads"
        "<br><input type=radio required name=similarityMeasure value='jaccard'" <<
        (similarityMeasure == "jaccard" ? " checked=on" : "") << ">Jaccard similarity"
        "<br><input type=radio required name=similarityMeasure value='correctedJaccard'" <<
        (similarityMeasure == "correctedJaccard" ? " checked=on" : "") << ">Corrected Jaccard similarity"

        "<br><input type=text name=referenceAnchorId size=6 style='text-align:center'";
    if(not referenceAnchorIdString.empty()) {
        html << " value='" << referenceAnchorIdString + "'";
    }
    html << "> Reference anchor id</div>";

    // Vertex coloring using assembly annotations.
    html <<
        "<input type=radio required name=vertexColoring value='byAssemblyAnnotations'" <<
        (vertexColoring == "byAssemblyAnnotations" ? " checked=on" : "") <<
        "> By annotations on assembly stage "
        "<input type=text name=assemblyStage style='text-align:center'";
        if(not assemblyStage.empty()) {
            html << " value='" << assemblyStage + "'";
        }
        html << " size=10>";


    html <<
        "<tr>"
        "<th>Edges"
        "<td class=left>"

        "<b>Edge coloring</b>"
		"<br><input type=radio required name=edgeColoring value='black'" <<
		(edgeColoring == "black" ? " checked=on" : "") << "> Black"
		"<br><input type=radio required name=edgeColoring value='byFlags'" <<
		(edgeColoring == "byFlags" ? " checked=on" : "") << "> By flags"
        "<br><input type=radio required name=edgeColoring value='random'" <<
        (edgeColoring == "random" ? " checked=on" : "") << "> Random"
        "<br><input type=radio required name=edgeColoring value='byCoverageLoss'" <<
        (edgeColoring == "byCoverageLoss" ? " checked=on" : "") << "> By coverage loss"
        "<hr>"

        "<b>Edge graphics</b>"

        "<br><input type=text name=edgeThickness style='text-align:center' required size=6 value=" <<
        edgeThickness << "> Thickness (arbitrary units)"

        "<br><input type=text name=minimumEdgeLength style='text-align:center' required size=6 value=" <<
        minimumEdgeLength << "> Minimum edge length (arbitrary units)"

        "<br><input type=text name=additionalEdgeLengthPerKb style='text-align:center' required size=6 value=" <<
        additionalEdgeLengthPerKb << "> Additional edge length per Kb (arbitrary units)"

        "<br><input type=text name=arrowSize style='text-align:center' required size=6 value=" <<
        arrowSize << "> Arrow size (arbitrary units)"

        "<hr>"
        "<input type=checkbox name=edgeLabels" <<
        (edgeLabels ? " checked" : "") <<
        "> Labels (dot layout only)";

}



void LocalAnchorGraph::writeHtml(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer)
{
    if((options.layoutMethod == "dot") and (options.vertexLabels or options.edgeLabels)) {

        // Use svg output from graphviz.
        writeHtml1(html, options, assemblyGraph3Pointer);

    } else {

        // Compute graph layout and use it to generate svg.
        writeHtml2(html, options, assemblyGraph3Pointer);

    }
}



// This is the code that uses svg output from graphviz.
void LocalAnchorGraph::writeHtml1(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer) const
{


        // Write it out in graphviz format.
        const string uuid = to_string(boost::uuids::random_generator()());
        const string dotFileName = tmpDirectory() + uuid + ".dot";
        writeGraphviz(dotFileName, options, assemblyGraph3Pointer);

        // Use graphviz to compute the layout.
        const string svgFileName = dotFileName + ".svg";
        const string shape = options.vertexLabels ? "rectangle" : "point";
        string command =
            options.layoutMethod +
            " -T svg " + dotFileName + " -o " + svgFileName +
            " -Nshape=" + shape +
            " -Gsize=" + to_string(options.sizePixels/72) + " -Gratio=expand ";
        if(options.vertexLabels) {
            command += " -Goverlap=false";
        }
        // cout << "Running command: " << command << endl;
        const int timeout = 30;
        bool timeoutTriggered = false;
        bool signalOccurred = false;
        int returnCode = 0;
        runCommandWithTimeout(command, timeout, timeoutTriggered, signalOccurred, returnCode);
        if(signalOccurred) {
            html << "Error during graph layout. Command was<br>" << endl;
            html << command;
            return;
        }
        if(timeoutTriggered) {
            html << "Timeout during graph layout." << endl;
            return;
        }
        if(returnCode!=0 ) {
            html << "Error during graph layout. Command was<br>" << endl;
            html << command;
            return;
        }
        std::filesystem::remove(dotFileName);



        // Write the svg to html.
        html << "<p><div style='border:solid;display:inline-block'>";
        ifstream svgFile(svgFileName);
        html << svgFile.rdbuf();
        svgFile.close();
        html << "</div>";

        // Remove the .svg file.
        std::filesystem::remove(svgFileName);

        // Add drag and zoom.
        addSvgDragAndZoom(html);
    }



// This is the code that computes the graph layout,
// then creates the svg.
void LocalAnchorGraph::writeHtml2(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer)
{
    // Use scientific notation because svg does not accept floating points
    // ending with a decimal point.
    html << std::scientific;

    computeLayout(options);
    computeLayoutBoundingBox();

    Box viewportBox = boundingBox;
    viewportBox.extend(0.05);
    viewportBox.makeSquare();

    // Begin the svg.
    const string svgId = "LocalAnchorGraph";
    html <<
        "\n<br><div style='display:inline-block;vertical-align:top;'>"
        "<svg id='" << svgId <<
        "' width='" <<  options.sizePixels <<
        "' height='" << options.sizePixels <<
        "' viewbox='" << viewportBox.xMin << " " << viewportBox.yMin << " " <<
        viewportBox.xSize() << " " <<
        viewportBox.ySize() << "'"
        " style='background-color:#f0f0f0'"
        ">\n";

    // Write the edges first so they don't obscure the vertices.
    writeEdges(html, options);

    // Write the vertices.
    writeVertices(html, options, assemblyGraph3Pointer);

    // Finish the svg.
    html << "</svg></div>";

    // Side panel.
    html << "<div style='display:inline-block;margin-left:20px'>";
    writeSvgControls(html, options);
    html << "</div>";
}



void LocalAnchorGraph::computeLayout(const LocalAnchorGraphDisplayOptions& options)
{
    const LocalAnchorGraph& graph = *this;
    SHASTA_ASSERT(anchorGraphPointer);
    const AnchorGraph& anchorGraph = *anchorGraphPointer;


    // Create a map containing the desired length for each edge.
    std::map<edge_descriptor, double> edgeLengthMap;
    BGL_FORALL_EDGES(eL, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edgeL = graph[eL];
        const AnchorGraphEdge& edgeG = anchorGraph[edgeL.eG];
        const uint64_t offset = edgeG.offset;

        const double displayLength =
            options.minimumEdgeLength +
            options.additionalEdgeLengthPerKb * 0.001 * double(offset);
        edgeLengthMap.insert({eL, displayLength});
    }

    // Compute the graph layout.
    layout.clear();
    const double timeout = 30.;
    if(options.layoutMethod == "custom") {
        const int quality = 2;
        computeLayoutCustom(
            graph,
            edgeLengthMap,
            layout,
            quality,
            timeout);
    } else {
        const string additionalOptions = "";
        computeLayoutGraphviz(
            graph,
            options.layoutMethod,
            timeout,
            layout,
            additionalOptions,
            &edgeLengthMap);
        // If the layout is dot, reverse the y coordinates so the arrows point down.
        if(options.layoutMethod == "dot") {
            for(auto& p: layout) {
                auto& y = p.second[1];
                y = -y;
            }
        }
    }
}



void LocalAnchorGraph::computeLayoutBoundingBox()
{

    boundingBox.xMin = std::numeric_limits<double>::max();
    boundingBox.xMax = std::numeric_limits<double>::min();
    boundingBox.yMin = boundingBox.xMin;
    boundingBox.yMax = boundingBox.xMax;
    for(const auto& p: layout) {
        const array<double, 2>& xy = p.second;
        const double x = xy[0];
        const double y = xy[1];
        boundingBox.xMin = min(boundingBox.xMin, x);
        boundingBox.xMax = max(boundingBox.xMax, x);
        boundingBox.yMin = min(boundingBox.yMin, y);
        boundingBox.yMax = max(boundingBox.yMax, y);
    }

}



void LocalAnchorGraph::Box::makeSquare()
{
    if(xSize() > ySize()) {
        const double delta = (xSize() - ySize()) / 2.;
        yMin -= delta;
        yMax += delta;
    } else {
        const double delta = (ySize() - xSize()) / 2.;
        xMin -= delta;
        xMax += delta;
    }
}



void LocalAnchorGraph::Box::extend(double factor)
{
    const double extend = factor * max(xSize(), ySize());
    xMin -= extend;
    xMax += extend;
    yMin -= extend;
    yMax += extend;
}



void LocalAnchorGraph::writeVertices(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options,
    const AssemblyGraphPostprocessor* assemblyGraph3Pointer) const
{
    const LocalAnchorGraph& graph = *this;

    if(options.vertexColoring == "byAssemblyAnnotations") {
        SHASTA_ASSERT(assemblyGraph3Pointer);
    }

    const double scalingFactor =
        (options.layoutMethod == "sfdp") ? 0.002 : 0.01;

    // Get the reference anchor, if needed.
    AnchorId referenceAnchorId = invalid<AnchorId>;
    if(options.vertexColoring == "byReadComposition") {
        referenceAnchorId = anchorIdFromString(options.referenceAnchorIdString);
        if((referenceAnchorId == invalid<AnchorId>) or (referenceAnchorId >= anchors.size())) {
            throw runtime_error("Invalid reference anchor id " + options.referenceAnchorIdString +
                ". Must be a number between 0 and " +
                to_string(anchors.size() / 2 - 1) + " followed by + or -.");
        }
    }
    const uint64_t referenceAnchorIdCoverage = anchors[referenceAnchorId].coverage();

    html << "\n<g id='vertices' style='stroke:none'>";

    BGL_FORALL_VERTICES(v, graph, LocalAnchorGraph) {
        const LocalAnchorGraphVertex& vertex = graph[v];
        const AnchorId anchorId = vertex.anchorId;
        const string anchorIdString = anchorIdToString(anchorId);
        const uint64_t coverage = anchors[anchorId].coverage();

        // Get the position of this vertex in the computed layout.
        const auto it = layout.find(v);
        SHASTA_ASSERT(it != layout.end());
        const auto& p = it->second;
        const double x = p[0];
        const double y = p[1];

        AnchorPairInfo info;
        if(options.vertexColoring == "byReadComposition") {
            anchors.analyzeAnchorPair(referenceAnchorId, anchorId, info);
        }

        // Get annotation information, if needed.
        bool hasVertexAnnotation = false;
        vector<AssemblyGraph::edge_descriptor> annotationEdges;
        if(options.vertexColoring == "byAssemblyAnnotations") {
            hasVertexAnnotation = assemblyGraph3Pointer->hasVertexAnnotation(anchorId);
            if(not hasVertexAnnotation) {
                assemblyGraph3Pointer->findAnnotationEdges(anchorId, annotationEdges);
            }
        }


        // Choose the color for this vertex.
        string color;
        if(vertex.distance == maxDistance) {
            color = "Cyan";
        } else if(vertex.distance == 0) {
            color = "Blue";
        } else if(options.vertexColoring == "byReadComposition") {

            // Color by similarity of read composition with the reference Anchor.
            if(options.vertexColoring == "byReadComposition") {

                double hue = 1.;    // 0=red, 1=green.
                if(options.similarityMeasure == "commonCount") {
                    // By common count.
                    hue = double(info.common) / double(referenceAnchorIdCoverage);

                } else if(options.similarityMeasure == "jaccard") {
                    // By Jaccard similarity.
                    hue = info.jaccard();
                } else {
                    // By corrected Jaccard similarity.
                    hue = info.correctedJaccard();
                }
                color = "hsl(" + to_string(uint32_t(std::round(hue * 120.))) +
                    ",100%,50%)";

            }
        } else if(options.vertexColoring == "byAssemblyAnnotations") {
            if(hasVertexAnnotation) {
                color = "Red";
            } else {
                if(annotationEdges.size() > 1) {
                    color = "Green";    // Multiple segments
                } else if(annotationEdges.size() == 1) {
                    const uint64_t segmentId = (*assemblyGraph3Pointer)[annotationEdges.front()].id;
                    const uint32_t hashValue = MurmurHash2(&segmentId, sizeof(segmentId), 759);
                    const uint32_t hue = hashValue % 360;
                    color = "hsl(" + to_string(hue) + ",50%,50%)";
                } else {
                    color = "Black";
                }
            }
        } else {
            color = "Black";
        }



        // Hyperlink.
        html << "\n<a href='exploreAnchor?anchorIdString=" <<
            HttpServer::urlEncode(anchorIdString) << "'>";

        // Write the vertex.
        html << "<circle cx='" << x << "' cy='" << y <<
            "' fill='" << color <<
            "' r='" << options.vertexSize * (scalingFactor * double(coverage)) <<
            "' id='" << anchorIdString << "'>"
            "<title>" << anchorIdString << ", coverage " << coverage;
        if(options.vertexColoring == "byReadComposition") {
            html << ", common " << info.common << ", J " <<
                std::fixed << std::setprecision(2) << info.jaccard() <<
                ", J' " << info.correctedJaccard();
            if(info.common > 0) {
                html << ", offset " << info.offsetInBases;
            }
        }
        if(options.vertexColoring == "byAssemblyAnnotations") {
            if(hasVertexAnnotation) {
                html << ", assembly graph vertex";
            } else if(annotationEdges.size() > 1) {
                html << ", multiple segments";
            } else if(annotationEdges.size() == 1) {
                const uint64_t segmentId = (*assemblyGraph3Pointer)[annotationEdges.front()].id;
                html << ", segment " << segmentId;
            }
        }
        html << "</title></circle>";

        // End the hyperlink.
        html << "</a>";
    }
    html << "\n</g>";
}




void LocalAnchorGraph::writeEdges(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& options) const
{
    const LocalAnchorGraph& graph = *this;
    SHASTA_ASSERT(anchorGraphPointer);
    const AnchorGraph& anchorGraph = *anchorGraphPointer;

    const double scalingFactor =
        (options.layoutMethod == "sfdp") ? 0.001 : 0.005;

    html << "\n<g id=edges>";

    BGL_FORALL_EDGES(eL, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edgeL = graph[eL];
        const AnchorGraphEdge& edgeG = anchorGraph[edgeL.eG];
        const AnchorPair& anchorPair = edgeG.anchorPair;
        const uint64_t coverage = anchorPair.orientedReadIds.size();
        const uint64_t offset = edgeG.offset;

        const vertex_descriptor v0 = source(eL, graph);
        const vertex_descriptor v1 = target(eL, graph);

        // Get the position of these vertices in the computed layout.
        const auto it0 = layout.find(v0);
        SHASTA_ASSERT(it0 != layout.end());
        const auto& p0 = it0->second;
        const double x0 = p0[0];
        const double y0 = p0[1];
        const auto it1 = layout.find(v1);
        SHASTA_ASSERT(it1 != layout.end());
        const auto& p1 = it1->second;
        const double x1 = p1[0];
        const double y1 = p1[1];

        const LocalAnchorGraphVertex& vertex0 = graph[v0];
        const AnchorId anchorId0 = vertex0.anchorId;
        const string anchorIdString0 = anchorIdToString(anchorId0);
        const LocalAnchorGraphVertex& vertex1 = graph[v1];
        const AnchorId anchorId1 = vertex1.anchorId;
        const string anchorIdString1 = anchorIdToString(anchorId1);

        string color = "Black";

        if(options.edgeColoring == "byFlags") {
			if(edgeG.isParallelEdge) {
				color = "Orange";
			} else if(edgeG.addedAtDeadEnd) {
				color = "Magenta";
			} else if(not edgeG.useForAssembly) {
				color = "GreenYellow";
			}
        }

        else if(options.edgeColoring == "random") {
            // To decide the color, hash the AnchorIds.
            // This way we always get the same color for the same edge.
            const auto p = make_pair(anchorId0, anchorId1);
            const uint32_t hashValue = MurmurHash2(&p, sizeof(p), 759);
            const uint32_t hue = hashValue % 360;
            color = "hsl(" + to_string(hue) + ",50%,50%)";
        }

        // Hyperlink.
        html << "\n<a href='" << anchorPair.url() << "'>";
        /*
        html << "\n<a href='exploreAnchorPair?"
            "anchorIdAString=" << HttpServer::urlEncode(anchorIdString0) << "&"
            "anchorIdBString=" << HttpServer::urlEncode(anchorIdString1) << "'>";
        */

        html <<
            "\n<line x1='" << x0 << "' y1='" << y0 <<
            "' x2='" << x1 << "' y2='" << y1 <<
            "' stroke='" << color <<
            "' stroke-width='" << scalingFactor * options.edgeThickness * double(coverage) <<
            "'>"
            "<title>" <<
            anchorIdString0 << " to " << anchorIdString1 <<
            ", coverage " << coverage << ", offset " << offset;
        html << "</title>""</line>";

        // End the hyperlink.
        html << "</a>";
    }
    html << "</g>";



    // Write the "arrows" to show edge directions.
    // They are just short lines near the target vertex of each edge.
    html << "\n<g id=arrows";
    if((options.edgeColoring == "black") or (options.edgeColoring == "byFlags")) {
        html << " stroke=white";
    } else {
        html << " stroke=black";
    }
    html << ">";
    BGL_FORALL_EDGES(eL, graph, LocalAnchorGraph) {
        const LocalAnchorGraphEdge& edgeL = graph[eL];
        const AnchorGraphEdge& edgeG = anchorGraph[edgeL.eG];
        const AnchorPair& anchorPair = edgeG.anchorPair;
        const uint64_t coverage = anchorPair.orientedReadIds.size();

        const vertex_descriptor v0 = source(eL, graph);
        const vertex_descriptor v1 = target(eL, graph);

        // Get the position of these vertices in the computed layout.
        const auto it0 = layout.find(v0);
        SHASTA_ASSERT(it0 != layout.end());
        const auto& p0 = it0->second;
        const double x0 = p0[0];
        const double y0 = p0[1];
        const auto it1 = layout.find(v1);
        SHASTA_ASSERT(it1 != layout.end());
        const auto& p1 = it1->second;
        const double x1 = p1[0];
        const double y1 = p1[1];

        const double relativeArrowLength = 0.3;
        const double x2 = (1. - relativeArrowLength) * x1 + relativeArrowLength * x0;
        const double y2 = (1. - relativeArrowLength) * y1 + relativeArrowLength * y0;

        html <<
            "\n<line x1='" << x1 << "' y1='" << y1 <<
            "' x2='" << x2 << "' y2='" << y2 <<
            "' stroke-width='" << 0.2 * scalingFactor * options.edgeThickness * double(coverage) <<
            "' />";

    }
    html << "</g>";
}



void LocalAnchorGraph::writeSvgControls(
    ostream& html,
    const LocalAnchorGraphDisplayOptions& /* options */) const
{
    html <<
        "<p><table>";

    // Add drag and zoom.
    addSvgDragAndZoom(html);



    // Buttons to change vertex size.
    html << R"stringDelimiter(
    <tr><th class=left>Vertex size<td>
    <button type='button' onClick='changeVertexSize(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeVertexSize(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeVertexSize(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeVertexSize(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeVertexSize(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeVertexSize(10.)' style='width:3em'>+++</button>
        <script>
        function changeVertexSize(factor)
        {
            var vertexGroup = document.getElementById('vertices');
            var vertices = vertexGroup.getElementsByTagName('circle');
            for(i=0; i<vertices.length; i++) {
                v = vertices[i];
                v.setAttribute('r', factor * v.getAttribute('r'));
            }
        }
        </script>
        )stringDelimiter";



    // Buttons to change edge thickness.
    html << R"stringDelimiter(
    <tr><th class=left>Edge thickness<td>
    <button type='button' onClick='changeThickness(0.1)' style='width:3em'>---</button>
    <button type='button' onClick='changeThickness(0.5)' style='width:3em'>--</button>
    <button type='button' onClick='changeThickness(0.8)' style='width:3em'>-</button>
    <button type='button' onClick='changeThickness(1.25)' style='width:3em'>+</button>
    <button type='button' onClick='changeThickness(2.)' style='width:3em'>++</button>
    <button type='button' onClick='changeThickness(10.)' style='width:3em'>+++</button>
        <script>
        function changeThickness(factor)
        {
            var edgeGroup = document.getElementById('edges');
            var edges = edgeGroup.getElementsByTagName('line');
            for(i=0; i<edges.length; i++) {
                e = edges[i];
                e.setAttribute('stroke-width', factor * e.getAttribute('stroke-width'));
            }

            var arrowsGroup = document.getElementById('arrows');
            var arrows = arrowsGroup.getElementsByTagName('line');
            for(i=0; i<arrows.length; i++) {
                a = arrows[i];
                a.setAttribute('stroke-width', factor * a.getAttribute('stroke-width'));
            }
        }
        </script>
        )stringDelimiter";



    // Zoom buttons.
    html << R"stringDelimiter(
        <tr title='Or use the mouse wheel.'><th class=left>Zoom<td>
        <button type='button' onClick='zoomSvg(0.1)' style='width:3em'>---</button>
        <button type='button' onClick='zoomSvg(0.5)' style='width:3em'>--</button>
        <button type='button' onClick='zoomSvg(0.8)' style='width:3em'>-</button>
        <button type='button' onClick='zoomSvg(1.25)' style='width:3em'>+</button>
        <button type='button' onClick='zoomSvg(2.)' style='width:3em'>++</button>
        <button type='button' onClick='zoomSvg(10.)' style='width:3em'>+++</button>
    )stringDelimiter";



    // Buttons to highlight an anchor and zoom to an anchor.
    html << R"stringDelimiter(
        <tr><td colspan=2>
        <button onClick='highlightAnchor()'>Highlight</button>
        <button onClick='zoomToAnchor()'>Zoom to</button>anchor
        <input id=selectedAnchorId type=text size=10 style='text-align:center'>
    <script>
    function zoomToAnchor()
    {
        // Get the anchor id from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        zoomToGivenAnchor(anchorId);
    }
    function zoomToGivenAnchor(anchorId)
    {
        var element = document.getElementById(anchorId);
        // Find the bounding box and its center.
        var box = element.getBBox();
        var xCenter = box.x + 0.5 * box.width;
        var yCenter = box.y + 0.5 * box.height;

        // Change the viewbox of the svg to be a bit larger than a square
        // containing the bounding box.
        var enlargeFactor = 5.;
        var size = enlargeFactor * Math.max(box.width, box.height);
        var factor = size / width;
        width = size;
        height = size;
        x = xCenter - 0.5 * size;
        y = yCenter - 0.5 * size;
        var svg = document.querySelector('svg');
        svg.setAttribute('viewBox', `${x} ${y} ${size} ${size}`);
        ratio = size / svg.getBoundingClientRect().width;
        svg.setAttribute('font-size', svg.getAttribute('font-size') * factor);
    }
    function highlightAnchor()
    {
        // Get the anchor id  from the input field.
        var anchorId = document.getElementById("selectedAnchorId").value;
        var element = document.getElementById(anchorId);

        element.style.fill = "Magenta";
    }
    </script>
    )stringDelimiter";


    html << "</table>";

    // Scroll down to the svg.
    const string svgId = "LocalAnchorGraph";
    html <<
        "<script>"
        "document.getElementById('" << svgId << "').scrollIntoView({block:'center'});"
        "</script>";

    html <<
        "<p>Use Ctrl+Click to pan."
        "<p>Use Ctrl-Wheel or the above buttons to zoom.";
}

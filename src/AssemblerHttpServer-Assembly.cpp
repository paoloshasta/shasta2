// Shasta.
#include "Assembler.hpp"
#include "Anchor.hpp"
#include "areSimilarSequences.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "deduplicate.hpp"
#include "findConvergingVertex.hpp"
#include "GTest.hpp"
#include "LocalAssembly6.hpp"
#include "LocalAssembly7.hpp"
#include "Markers.hpp"
#include "RestrictedAnchorGraph.hpp"
#include "SegmentStepSupport.hpp"
#include "Superbubble.hpp"
#include "TangleMatrix1.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/algorithm/string.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/tokenizer.hpp>

// Standard library.
#include "fstream.hpp"
#include "tuple.hpp"



void Assembler::exploreSegments(
    const vector<string>& request,
    ostream& html)
{
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    html << "<h2>Assembly graph segments</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "</table>"
        "<input type=submit value='Get information'>"
        "</form>";

    if(assemblyStage.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    html <<
        "<h2>Assembly graph at stage " << assemblyStage << " </h2>"
        "<p>The assembly graph at stage " << assemblyStage <<
        " has " << num_edges(assemblyGraph) << " edges (segments)." << endl;

    html << "<table><tr><th>Vertex<br>(segment)<br>id<th>Number<br>of<br>steps"
        "<th>Average<br>coverage<th>Estimated<br>length<th>Actual<br>length";

    BGL_FORALL_EDGES(e, assemblyGraph, AssemblyGraph) {
        const AssemblyGraphEdge& edge = assemblyGraph[e];
        const uint64_t coverage = uint64_t(std::round(edge.averageCoverage()));
        const string url = "exploreSegmentSteps?assemblyStage=" + assemblyStage + "&segmentName=" + to_string(edge.id);
        html <<
            "<tr>"
            "<td class=centered><a href='" << url << "'>" << edge.id << "</a>"
            "<td class=centered>" << edge.size() <<
            "<td class=centered>" << coverage <<
            "<td class=centered>" << edge.offset() <<
            "<td class=centered>";
        if(edge.wasAssembled) {
            html << edge.sequenceLength();
        }
    }

    html << "</table>";
}



void Assembler::exploreSegmentSequence(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    uint64_t sequenceBegin = 0;
    const bool sequenceBeginIsPresent = getParameterValue(request, "sequenceBegin", sequenceBegin);

    uint64_t sequenceEnd = 0;
    const bool sequenceEndIsPresent = getParameterValue(request, "sequenceEnd", sequenceEnd);


    // Start the form.
    html << "<h2>Assembled segment sequence</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage size=8 style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=8>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName size=8 style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";

    html << "<tr><th>Begin position<td class=centered>"
        "<input type=text name=sequenceBegin style='text-align:center' size=8";
    if(sequenceBeginIsPresent) {
        html << " value=" << sequenceBegin;
    }
    html << ">";

    html << "<tr><th>End position<td class=centered>"
        "<input type=text name=sequenceEnd style='text-align:center' size=8";
    if(sequenceEndIsPresent) {
        html << " value=" << sequenceEnd;
    }
    html << ">";
\
    // End the form.
    html <<
        "</table>"
        "<p><input type=submit value='Get segment sequence'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];


    if(not edge.wasAssembled) {
        html << "<p>Assembled sequence for segment " << segmentId <<
            " at assembly stage " << assemblyStage << " is not available.";
        return;
    }


    html << "<h2>Assembled sequence of segment " << segmentId << " at assembly stage " << assemblyStage << "</h2>";

    vector<Base> sequence;
    edge.getSequence(sequence);

    if(not sequenceBeginIsPresent) {
        sequenceBegin = 0;
    }
    if(not sequenceEndIsPresent) {
        sequenceEnd = sequence.size();
    }

    if(sequenceBegin >= sequence.size()) {
        sequenceBegin = sequence.size() - 1;
    }
    if(sequenceEnd > sequence.size()) {
        sequenceEnd = sequence.size();
    }
    if(sequenceEnd < sequenceBegin) {
        sequenceEnd = sequenceBegin;
    }

    html << "<div style='font-family:monospace'>";
    html << ">" << segmentName << "-" << sequenceBegin << "-" << sequenceEnd <<
        ", length " << sequenceEnd - sequenceBegin << "<br>";
    copy(sequence.begin() + sequenceBegin, sequence.begin() + sequenceEnd,
        ostream_iterator<Base>(html));
    html << "</div>";

    // Also write the sequence to LocalAssembly.fasta.
    ofstream fasta("Segment.fasta");
    fasta << ">" << segmentName << "-" << sequenceBegin << "-" << sequenceEnd <<
        " length " << sequenceEnd - sequenceBegin << endl;
    copy(sequence.begin() + sequenceBegin, sequence.begin() + sequenceEnd,
        ostream_iterator<Base>(fasta));

    html << "<p>The sequence was stored in Segment.fasta.";

}



void Assembler::exploreSegmentSteps(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string displaySteps = "none";
    HttpServer::getParameterValue(request, "displaySteps", displaySteps);

    string stepBeginString;
    HttpServer::getParameterValue(request, "stepBegin", stepBeginString);

    string stepEndString;
    HttpServer::getParameterValue(request, "stepEnd", stepEndString);

    string firstStepsCountString = "5";
    HttpServer::getParameterValue(request, "firstStepsCount", firstStepsCountString);

    string lastStepsCountString = "5";
    HttpServer::getParameterValue(request, "lastStepsCount", lastStepsCountString);

    string showSequenceDetailsString;
    const bool showSequenceDetails = getParameterValue(request, "showSequenceDetails", showSequenceDetailsString);



    // Start the form.
    html << "<h2>Assembly graph segment</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=8>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << " size=8>";



    // Options to control which segment steps are shown.
    html <<
        "<tr>"
        "<th class=left>Show segment steps"
        "<td class=left>"

        "<input type=radio required name=displaySteps value='none'" <<
        (displaySteps == "none" ? " checked=on" : "") << "> None"

        "<br><input type=radio required name=displaySteps value='all'" <<
        (displaySteps == "all" ? " checked=on" : "") << "> All"

        "<br><input type=radio required name=displaySteps value='range'" <<
        (displaySteps == "range" ? " checked=on" : "") << "> Steps in position range "
        "<input type=text name=stepBegin size=8 style='text-align:center' value='" << stepBeginString << "'> to "
        "<input type=text name=stepEnd size=8 style='text-align:center' value='" << stepEndString << "'>"

        "<br><input type=radio required name=displaySteps value='first'" <<
        (displaySteps == "first" ? " checked=on" : "") << "> First "
        "<input type=text name=firstStepsCount size=8 style='text-align:center' value='" << firstStepsCountString << "'>"
        " steps"

        "<br><input type=radio required name=displaySteps value='last'" <<
        (displaySteps == "last" ? " checked=on" : "") << "> Last "
        "<input type=text name=lastStepsCount size=8 style='text-align:center' value='" << lastStepsCountString << "'>"
        " steps"

        "<br><input type=checkbox name=showSequenceDetails" << (showSequenceDetails ? " checked" : "") <<
        "> Also show sequence details for these steps";
        ;



    // End the form.
    html <<
        "</table>"
        "<p><input type=submit value='Get segment information'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    const AssemblyGraph::vertex_descriptor v0 = source(e, assemblyGraph);
    const AssemblyGraph::vertex_descriptor v1 = target(e, assemblyGraph);
    const AnchorId anchorId0 = assemblyGraph[v0].anchorId;
    const AnchorId anchorId1 = assemblyGraph[v1].anchorId;



    html << "<h2>Segment " << segmentId << " at assembly stage " << assemblyStage << "</h2>";

    // Summary table.
    html <<
        "<table>"
        "<tr><th class=left>First anchor<td class = centered>" << anchorIdToString(anchorId0) <<
        "<tr><th class=left>Last anchor<td class = centered>" << anchorIdToString(anchorId1) <<
        "<tr><th class=left>Number of steps<td class = centered>" << edge.size();
    if(not edge.empty()) {
        html <<
            "<tr><th class=left>Average coverage<td class = centered>" << uint64_t(std::round(edge.averageCoverage()));
    }
    html <<
        "<tr><th class=left>Estimated length<td class = centered>" << edge.offset() <<
        "<tr><th class=left>Assembled<td class = centered>" << (edge.wasAssembled ? "Yes" : "No");
    if(edge.wasAssembled) {
        html <<
            "<tr><th class=left>Assembled length<td class = centered>" << edge.sequenceLength();

    }
    html << "</table>";

    html <<
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(anchorIdToString(anchorId0)) <<
        "'>See the first anchor in the local anchor graph</a>"
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(anchorIdToString(anchorId1)) <<
        "'>See the last anchor in the local anchor graph</a>"
        "<br><a href='exploreLocalAnchorGraph?anchorIdsString=" <<
        HttpServer::urlEncode(
            anchorIdToString(anchorId0) + " " +
            anchorIdToString(anchorId1)) <<
        "'>See the first and last anchor in the local anchor graph</a>";


    // Figure out the step position range to use.
    uint64_t stepBegin = invalid<uint64_t>;
    uint64_t stepEnd = invalid<uint64_t>;
    if(displaySteps == "all") {
        stepBegin = 0;
        stepEnd = edge.size();
    } else if(displaySteps == "range") {
        try {
            stepBegin = atoul(stepBeginString);
        } catch(std::exception& e) {
            throw runtime_error("Begin " + stepBeginString + " is not valid. Must be a number.");
        }
        try {
            stepEnd = atoul(stepEndString);
        } catch(std::exception& e) {
            throw runtime_error("End " + stepEndString + " is not valid. Must be a number.");
        }
        if(stepBegin >= edge.size()) {
            stepBegin = edge.size() - 1;
        }
        if(stepBegin > edge.size()) {
            stepBegin = edge.size();
        }
        if(stepEnd < stepBegin) {
            stepEnd = stepBegin + 1;
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "first") {
        stepBegin = 0;
        try {
            stepEnd = atoul(firstStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("First anchors count " + firstStepsCountString + " is not valid. Must be a number.");
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "last") {
        stepEnd = edge.size();
        uint64_t count = invalid<uint64_t>;
        try {
            count = atoul(lastStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("Last anchors count " + lastStepsCountString + " is not valid. Must be a number.");
        }
        if(count > edge.size()) {
            stepBegin = 0;
        } else {
            stepBegin = stepEnd - count;
        }
    }



    // Details table showing the requested steps.
    if((displaySteps != "none") and (not edge.empty())) {

        html <<
            "<br>"
            "<table>"
            "<tr><th>Step<th>AnchorIdA<th>AnchorIdB<th>Coverage<th>Estimated<br>Length";
        if(edge.wasAssembled) {
            html << "<th>Actual<br>Length";
            if(showSequenceDetails) {
                html <<
                    "<th>Sequence<br>begin"
                    "<th>Sequence<br>end"
                    "<th>Sequence";
            }
        }

        uint64_t sequencePosition = 0;
        if(showSequenceDetails) {
            for(uint64_t i=0; i<stepBegin; i++) {
                sequencePosition += edge[i].sequence.size();
            }
        }

        for(uint64_t stepId=stepBegin; stepId!=stepEnd; ++stepId) {
            const AssemblyGraphEdgeStep& step = edge[stepId];
            const string url = "exploreSegmentStep?assemblyStage=" +
                assemblyStage + "&segmentName=" + to_string(segmentId) + "&stepId=" + to_string(stepId);

            html <<
                "<tr>"
                "<td class=centered><a href='" << url << "'>" << stepId << "</a>"
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdA) <<
                "<td class=centered>" << anchorIdToString(step.anchorPair.anchorIdB) <<
                "<td class=centered>" << step.anchorPair.orientedReadIds.size() <<
                "<td class=centered>" << step.offset;
            if(edge.wasAssembled) {
                html << "<td class=centered>" << step.sequence.size();
                if(showSequenceDetails) {
                    html <<
                        "<td class=centered>" << sequencePosition <<
                        "<td class=centered>" << sequencePosition + step.sequence.size() <<
                        "<td style='font-family:monospace'>";
                    copy(step.sequence.begin(), step.sequence.end(), ostream_iterator<Base>(html));
                    sequencePosition += step.sequence.size();
                }
            }


        }

        html << "</table>";

        // Link to the local anchor graph showing these anchors.
        {
            const AssemblyGraphEdgeStep& firstStep = edge[stepBegin];
            const AnchorId anchorIdA = firstStep.anchorPair.anchorIdA;
            string urlAnchors = anchorIdToString(anchorIdA);
            for(uint64_t stepId=stepBegin; stepId!=stepEnd; ++stepId) {
                const AssemblyGraphEdgeStep& step = edge[stepId];
                const AnchorId anchorIdB = step.anchorPair.anchorIdB;
                urlAnchors += " ";
                urlAnchors += anchorIdToString(anchorIdB);
            }

            html << "<br><a href='exploreLocalAnchorGraph?anchorIdsString="<< HttpServer::urlEncode(urlAnchors) <<
            "'>See these anchors in the local anchor graph</a>";
        }
    }
}



void Assembler::exploreSegmentStepSupport(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string displaySteps = "none";
    HttpServer::getParameterValue(request, "displaySteps", displaySteps);

    string stepBeginString;
    HttpServer::getParameterValue(request, "stepBegin", stepBeginString);

    string stepEndString;
    HttpServer::getParameterValue(request, "stepEnd", stepEndString);

    string firstStepsCountString = "5";
    HttpServer::getParameterValue(request, "firstStepsCount", firstStepsCountString);

    string lastStepsCountString = "5";
    HttpServer::getParameterValue(request, "lastStepsCount", lastStepsCountString);

    double pixelsPerBase = 0.001;
    HttpServer::getParameterValue(request, "pixelsPerBase", pixelsPerBase);

    double pixelsPerStep = 20.;
    HttpServer::getParameterValue(request, "pixelsPerStep", pixelsPerStep);


    // Start the form.
    html << "<h2>Segment steps support</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=8>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << " size=8>";



    // Options to control which segment steps are shown.
    html <<
        "<tr>"
        "<th class=left>Show segment steps"
        "<td class=left>"

        "<input type=radio required name=displaySteps value='all'" <<
        (displaySteps == "all" ? " checked=on" : "") << "> All"

        "<br><input type=radio required name=displaySteps value='range'" <<
        (displaySteps == "range" ? " checked=on" : "") << "> Steps in position range "
        "<input type=text name=stepBegin size=8 style='text-align:center' value='" << stepBeginString << "'> to "
        "<input type=text name=stepEnd size=8 style='text-align:center' value='" << stepEndString << "'>"

        "<br><input type=radio required name=displaySteps value='first'" <<
        (displaySteps == "first" ? " checked=on" : "") << "> First "
        "<input type=text name=firstStepsCount size=8 style='text-align:center' value='" << firstStepsCountString << "'>"
        " steps"

        "<br><input type=radio required name=displaySteps value='last'" <<
        (displaySteps == "last" ? " checked=on" : "") << "> Last "
        "<input type=text name=lastStepsCount size=8 style='text-align:center' value='" << lastStepsCountString << "'>"
        " steps"

        "<tr><th>Pixels per base"
        "<td class=centered>"
        "<input type=text name=pixelsPerBase size=8 style='text-align:center' value='" << pixelsPerBase << "'>"

        "<tr><th>Pixels per step"
        "<td class=centered>"
        "<input type=text name=pixelsPerStep size=8 style='text-align:center' value='" << pixelsPerStep << "'>"
        ;



    // End the form.
    html <<
        "</table>"
        "<p><input type=submit value='Get segment support'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }
    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    // Figure out the step position range to use.
    uint64_t stepBegin = invalid<uint64_t>;
    uint64_t stepEnd = invalid<uint64_t>;
    if(displaySteps == "all") {
        stepBegin = 0;
        stepEnd = edge.size();
    } else if(displaySteps == "range") {
        try {
            stepBegin = atoul(stepBeginString);
        } catch(std::exception& e) {
            throw runtime_error("Begin " + stepBeginString + " is not valid. Must be a number.");
        }
        try {
            stepEnd = atoul(stepEndString);
        } catch(std::exception& e) {
            throw runtime_error("End " + stepEndString + " is not valid. Must be a number.");
        }
        if(stepBegin >= edge.size()) {
            stepBegin = edge.size() - 1;
        }
        if(stepBegin > edge.size()) {
            stepBegin = edge.size();
        }
        if(stepEnd < stepBegin) {
            stepEnd = stepBegin + 1;
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "first") {
        stepBegin = 0;
        try {
            stepEnd = atoul(firstStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("First anchors count " + firstStepsCountString + " is not valid. Must be a number.");
        }
        if(stepEnd > edge.size()) {
            stepEnd = edge.size();
        }
    } else if(displaySteps == "last") {
        stepEnd = edge.size();
        uint64_t count = invalid<uint64_t>;
        try {
            count = atoul(lastStepsCountString);
        } catch(std::exception& e) {
            throw runtime_error("Last anchors count " + lastStepsCountString + " is not valid. Must be a number.");
        }
        if(count > edge.size()) {
            stepBegin = 0;
        } else {
            stepBegin = stepEnd - count;
        }
    }
    const uint64_t stepCount = stepEnd - stepBegin;


    // Find the OrientedReadIds that appear in the AnchorPairs of these steps.
    vector<OrientedReadId> orientedReadIds;
    for(uint64_t i=stepBegin; i<stepEnd; ++i) {
        const AssemblyGraphEdgeStep& step = edge[i];
        const AnchorPair anchorPair = step.anchorPair;
        std::ranges::copy(step.anchorPair.orientedReadIds, back_inserter(orientedReadIds));
    }
    deduplicate(orientedReadIds);



    // Find the base position offsets of the two anchors of each step
    // relative to the base position of the first anchor of the first step.
    vector< pair<uint64_t, uint64_t> > offsetTable(stepCount);
    uint64_t offset = 0;
    for(uint64_t step=stepBegin; step!=stepEnd; step++) {
        const uint64_t length = (edge.wasAssembled ? edge[step].sequence.size() : edge[step].offset);
        pair<uint64_t, uint64_t>& offsets = offsetTable[step - stepBegin];
        offsets.first = offset;
        offset += length;
        offsets.second = offset;
    }



    // Svg constants.
    // Avoid complications with aspect ratio adjustments in the browser.

    // Horizontal direction.
    const double maxOffsetBases = double(offsetTable.back().second);
    const double basesPerStep = pixelsPerStep / pixelsPerBase;
    const double borderPixels = 5.;
    const double borderBases = borderPixels / pixelsPerBase;
    const double svgWidthBases = maxOffsetBases + double(stepCount) * basesPerStep + 2. * borderBases;
    const double svgWidthPixels = svgWidthBases * pixelsPerBase;

    // Vertical direction.
    const double svgHeightPixels = 10.;
    const double svgHeightBases = svgHeightPixels / pixelsPerBase;
    const double yBases = svgHeightBases / 2.;

    // Horizontal and vertical direction.
    const double dotRadiusBases = yBases;



    // Begin the table.
    html <<"<br><table>";



    // Write a table row for the segment.
    html <<
        "<tr><th>Segment " << segmentId <<
        "<td class=centered style='vertical-align:middle'>"
        "<svg width='" << svgWidthPixels << "' height='" << svgHeightPixels << "'"
        " viewbox='0 0 " << svgWidthBases << " " << svgHeightBases << "'"
        " style='background-color:#f0f0f0'"
        ">";

    // Write a rectangle for each step.
    for(uint64_t step=stepBegin; step<stepEnd; step++) {
        const pair<uint64_t, uint64_t>& offsets = offsetTable[step - stepBegin];
        const double xBases = borderBases + double(offsets.first) + double(step - stepBegin) * basesPerStep;
        const double widthBases = double(offsets.second - offsets.first) + basesPerStep;
        const string color = "LightBlue";
        html <<
            "<g><title>Step " << step << "</title>"
            "<rect x='" << xBases << "' y='0' width='" << widthBases << "' height='" << svgHeightBases <<
            "' style='fill:" << color << "' /></g>";
    }

    // Write a dot for each Anchor.
    for(uint64_t step=stepBegin; step<=stepEnd; step++) {
        uint64_t offset;
        AnchorId anchorId;
        if(step == stepEnd) {
            offset = offsetTable[step - 1 - stepBegin].second;
            anchorId = edge[step - 1].anchorPair.anchorIdB;
        } else {
            offset = offsetTable[step - stepBegin].first;
            anchorId = edge[step].anchorPair.anchorIdA;
        }
        const double xBases = borderBases + double(offset) + double(step - stepBegin) * basesPerStep;
        html <<
            "<g><title>" << anchorIdToString(anchorId) << "</title>"
            "<circle cx='" << xBases << "' cy='" << yBases << "' r=" << dotRadiusBases << " /></g>";
    }

    html << "</svg>";



    // Write a table row for each OrientedReadId.
    for(const OrientedReadId orientedReadId: orientedReadIds) {
        html << "<tr><th>" << orientedReadId <<
            "<td class=centered style='vertical-align:middle'>"
            "<svg width='" << svgWidthPixels << "' height='" << svgHeightPixels << "'"
            " viewbox='0 0 " << svgWidthBases << " " << svgHeightBases << "'"
            " style='background-color:#f0f0f0'"
            ">";

        // Write a rectangle for each step.
        for(uint64_t step=stepBegin; step<stepEnd; step++) {
            const AnchorPair& anchorPair = edge[step].anchorPair;
            const bool containsOrientedRead = std::ranges::binary_search(anchorPair.orientedReadIds, orientedReadId);
            if(containsOrientedRead) {
                const pair<uint64_t, uint64_t>& offsets = offsetTable[step - stepBegin];
                const double xBases = borderBases + double(offsets.first) + double(step - stepBegin) * basesPerStep;
                const double widthBases = double(offsets.second - offsets.first) + basesPerStep;
                const string color = "LightPink";
                html <<
                    "<g><title>Step " << step << "</title>"
                    "<rect x='" << xBases << "' y='0' width='" << widthBases << "' height='" << svgHeightBases <<
                    "' style='fill:" << color << "' /></g>";
            }
        }

        // Write a dot for each Anchor.
            for(uint64_t step=stepBegin; step<=stepEnd; step++) {
                uint64_t offset;
                AnchorId anchorId;
                if(step == stepEnd) {
                    offset = offsetTable[step - 1 - stepBegin].second;
                    anchorId = edge[step - 1].anchorPair.anchorIdB;
                } else {
                    offset = offsetTable[step - stepBegin].first;
                    anchorId = edge[step].anchorPair.anchorIdA;
                }
                const string color = (anchors().anchorContains(anchorId, orientedReadId) ? "Black" : "LightGrey");
                const double xBases = borderBases + double(offset) + double(step - stepBegin) * basesPerStep;
                html <<
                    "<g><title>" << anchorIdToString(anchorId) << "</title>"
                    "<circle cx='" << xBases << "' cy='" << yBases << "' r=" << dotRadiusBases <<
                    " fill=\"" << color << "\""
                    " /></g>";
            }

        html << "</svg>";
    }

    html << "</svg>";



    // End the table.
    html << "</table>";



    // Use SegmentStepSupport instead.
    vector<SegmentStepSupport> stepSupports;
    SegmentStepSupport::get(assemblyGraph, e, uint32_t(stepBegin), uint32_t(stepEnd), stepSupports);
    std::ranges::sort(stepSupports, std::ranges::less(),
        [](const SegmentStepSupport& s) {return std::tie(s.orientedReadId, s.stepId);}
        );

    html << "<br><h3>Details</h3>";
    SegmentStepSupport::writeHtml(html, assemblyGraph, stepSupports);

}



void Assembler::exploreSegmentStep(
    const vector<string>& request,
    ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);

    string segmentName;
    HttpServer::getParameterValue(request, "segmentName", segmentName);

    string stepIdString;
    const bool stepIdStringIsPresent = HttpServer::getParameterValue(request, "stepId", stepIdString);
    boost::trim(stepIdString);

    int localAssemblyVersion = 7;
    getParameterValue(request, "localAssemblyVersion", localAssemblyVersion);



    // LocalAssembly7::Options.
    LocalAssembly7::Options localAssembly7Options;
    getParameterValue(request, "aExtend", localAssembly7Options.aExtend);
    getParameterValue(request, "bExtend", localAssembly7Options.bExtend);

    string methodString = "Adaptive";
    getParameterValue(request, "method", methodString);
    localAssembly7Options.setMethod(methodString);

    getParameterValue(request, "commonCoverageThreshold", localAssembly7Options.commonCoverageThreshold);

    string disallowFastPathString;
    localAssembly7Options.allowFastPath = not HttpServer::getParameterValue(request,
        "disallowFastPath", disallowFastPathString);

    getParameterValue(request, "fastPathFractionThreshold", localAssembly7Options.fastPathFractionThreshold);
    getParameterValue(request, "maxAbpoaLength", localAssembly7Options.maxAbpoaLength);



    // Start the form.
    html << "<h2>Local assembly</h2><form><table>";

    html <<
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "<tr>"
        "<th class=left>Segment name"
        "<td class=centered><input type=text name=segmentName style='text-align:center' required";
    if(not segmentName.empty()) {
        html << " value='" << segmentName + "'";
    }
    html << ">";

    html <<
        "<tr>"
        "<th class=left>Segment step id"
        "<td class=centered><input type=text name=stepId style='text-align:center' required";
    if(stepIdStringIsPresent) {
        html << " value='" << stepIdString + "'";
    }
    html << ">";

    html <<
        "<tr><th class=left>Local assembly version<td class=centered>"
        "<input type=radio name=localAssemblyVersion value=6" <<
        (localAssemblyVersion == 6 ? " checked=on" : "") << "> 6"
        "<br><input type=radio name=localAssemblyVersion value=7" <<
        (localAssemblyVersion == 7 ? " checked=on" : "") << "> 7";

    html <<
        "<tr><th class=left>aExtend<td class=centered>"
        "<input type=text name=aExtend style='text-align:center' "
        "value='" << localAssembly7Options.aExtend << "'>";

    html <<
        "<tr><th class=left>bExtend<td class=centered>"
        "<input type=text name=bExtend style='text-align:center' "
        "value='" << localAssembly7Options.bExtend << "'>";

    html <<
        "<tr><th class=left>Method<td class=left>"
        "<input type=radio name=method value=Adaptive" <<
        (localAssembly7Options.method == LocalAssembly7::Method::Adaptive ? " checked=on" : "") << "> Adaptive"
        "<br><input type=radio name=method value=Abpoa" <<
        (localAssembly7Options.method == LocalAssembly7::Method::Abpoa ? " checked=on" : "") << "> Abpoa"
        "<br><input type=radio name=method value=Poasta" <<
        (localAssembly7Options.method == LocalAssembly7::Method::Poasta ? " checked=on" : "") << "> Poasta"
        "<br><input type=radio name=method value=TheseusOnly" <<
        (localAssembly7Options.method == LocalAssembly7::Method::TheseusOnly ? " checked=on" : "") <<
        "> Theseus, using only oriented reads on both anchors"
        "<br><input type=radio name=method value=TheseusAll" <<
        (localAssembly7Options.method == LocalAssembly7::Method::TheseusAll ? " checked=on" : "") <<
        "> Theseus, using all oriented reads on one or both anchors."
        "<br><input type=radio name=method value=DeBruijn" <<
        (localAssembly7Options.method == LocalAssembly7::Method::DeBruijn ? " checked=on" : "") << "> De Bruijn"
        ;

    html <<
        "<tr><th class=left>commonCoverageThreshold<td class=centered>"
        "<input type=text name=commonCoverageThreshold style='text-align:center' "
        "value='" << localAssembly7Options.commonCoverageThreshold << "'>";

    html <<
        "<tr><th class=left>Forbid fast path"
        "<td class=centered>"
        "<input type=checkbox name=disallowFastPath" <<
            (localAssembly7Options.allowFastPath ? "" : " checked") <<
            ">";

    html <<
        "<tr><th class=left>fastPathFractionThreshold<td class=centered>"
        "<input type=text name=fastPathFractionThreshold style='text-align:center' "
        "value='" << localAssembly7Options.fastPathFractionThreshold << "'>";

    html <<
        "<tr><th class=left>maxAbpoaLength<td class=centered>"
        "<input type=text name=maxAbpoaLength style='text-align:center' "
        "value='" << localAssembly7Options.maxAbpoaLength << "'>";



    // End the form.
    html <<
        "</table>"
        "<br><input type=submit value='Run the local assembly'>"
        "</form>";

    if(segmentName.empty()) {
        return;
    }
    if(not stepIdStringIsPresent) {
        return;
    }

    uint64_t stepId;
    try {
        stepId = atoul(stepIdString);
    } catch(std::exception& e) {
        throw runtime_error("Step id " + stepIdString + " is not valid. Must be a number.");
    }

    uint64_t segmentId = invalid<uint64_t>;
    try {
        segmentId = std::stol(segmentName);
    } catch(exception&) {
    }
    if(segmentId == invalid<uint64_t>) {
        html << "Segment name must be a number.";
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    const AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);

    // Find the AssemblyGraphEdge corresponding to the requested segment.
    auto it = assemblyGraph.edgeMap.find(segmentId);
    if(it == assemblyGraph.edgeMap.end()) {
        html << "<p>Assembly graph at stage " << assemblyStage <<
            " does not have segment " << segmentId;
        return;
    }

    const AssemblyGraph::edge_descriptor e = it->second;
    const AssemblyGraphEdge& edge = assemblyGraph[e];

    if(stepId >= edge.size()) {
        html << "<p>Step " << stepId << " is not valid for this segment, which has " <<
            edge.size() << " steps.";
        return;
    }



    // Write the AnchorPair to html.
    html << "<h2>AnchorPair for step " << stepId << " of segment " << segmentId << " at assembly stage " <<
        assemblyStage << "</h2>";
    edge[stepId].anchorPair.writeAllHtml(html, anchors());



    // Write the local assembly to html.
    html << "<h2>Local assembly for step " << stepId << " of segment " << segmentId << " at assembly stage " <<
        assemblyStage << "</h2>";



    // Gather OrientedReadIds from the previous and next step.
    vector<OrientedReadId> additionalOrientedReadIds;
    if(stepId > 0) {
        std::ranges::copy(edge[stepId - 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    if(stepId < edge.size() - 1) {
        std::ranges::copy(edge[stepId + 1].anchorPair.orientedReadIds, back_inserter(additionalOrientedReadIds));
    }
    deduplicate(additionalOrientedReadIds);



    // Combine the Oriented reads in the AnchorPair
    // and the additional OrientedReadIds.
    vector<OrientedReadId> orientedReadIds = additionalOrientedReadIds;
    const AnchorPair& anchorPair = edge[stepId].anchorPair;
    std::ranges::copy(anchorPair.orientedReadIds, back_inserter(orientedReadIds));
    deduplicate(orientedReadIds);



    // Do the local assembly.
    switch(localAssemblyVersion) {
    case 6:
        {
            LocalAssembly6 localAssembly(
                anchors(),
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                html,
                orientedReadIds);
            return;
        }
    case 7:
        {
            LocalAssembly7 localAssembly(
                localAssembly7Options,
                anchors(),
                anchorPair.anchorIdA,
                anchorPair.anchorIdB,
                html,
                orientedReadIds);
            return;
        }
    default:
        SHASTA2_ASSERT(0);
    }

}



AssemblyGraphPostprocessor& Assembler::getAssemblyGraph(
    const string& assemblyStage,
    const Options& options)
{
    auto it = assemblyGraphTable.find(assemblyStage);
    if(it == assemblyGraphTable.end()) {
        shared_ptr<AssemblyGraphPostprocessor> p =
            make_shared<AssemblyGraphPostprocessor>(anchors(), journeys(), options, assemblyStage);
        tie(it, ignore) = assemblyGraphTable.insert(make_pair(assemblyStage, p));
    }
    return *(it->second);
}



void Assembler::exploreTangleMatrix(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);
    boost::trim(assemblyStage);

    string entrancesString;
    HttpServer::getParameterValue(request, "entrances", entrancesString);
    boost::trim(entrancesString);

    string exitsString;
    HttpServer::getParameterValue(request, "exits", exitsString);
    boost::trim(exitsString);

    double epsilon = httpServerData.options->detangleEpsilon;
    HttpServer::getParameterValue(request, "epsilon", epsilon);



    // Start the form.
    html << "<h2>Assembly graph tangle matrix (new)</h2><form>";

    html <<
        "<table>"
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=30>";

    html <<
        "<tr title='Enter assembly segment names separated by commas or spaces.'>"
        "<th class=left>Entrances"
        "<td class=centered>"
        "<input type=text name=entrances style='text-align:center'";
    if(not entrancesString.empty()) {
        html << " value='" << entrancesString << "'";
    }
    html << " size=30>";

    html <<
        "<tr title='Enter assembly segment names separated by commas or spaces.'>"
        "<th class=left>Exits"
        "<td class=centered>"
        "<input type=text name=exits style='text-align:center'";
    if(not exitsString.empty()) {
        html << " value='" << exitsString << "'";
    }
    html << " size=30>";

    html << "<tr><th class=left>Epsilon for G-test evaluation"
        "<td class=centered>"
        "<input type=text name=epsilon style='text-align:center' value='" << epsilon << "'>";

    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Compute tangle matrix'>"
        "</form>";

    if(entrancesString.empty() or exitsString.empty()) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);



    // Find AssemblyGraph edges corresponding to the entrances.
    vector<AssemblyGraph::edge_descriptor> entrances;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(entrancesString, boost::char_separator<char>(", "));
        for(const string& vertexIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(vertexIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << vertexIdString << ". Must be a number.";
                return;
            }

            // Find the AssemblyGraphEdge corresponding to the requested segment.
            auto it = assemblyGraph.edgeMap.find(segmentId);
            if(it == assemblyGraph.edgeMap.end()) {
                html << "<p>Assembly graph at stage " << assemblyStage <<
                    " does not have segment " << segmentId;
                return;
            }

            entrances.push_back(it->second);
        }

    }
    std::ranges::sort(entrances, assemblyGraph.orderById);



    // Find AssemblyGraph edge corresponding to the exits.
    vector<AssemblyGraph::edge_descriptor> exits;
    {
        boost::tokenizer< boost::char_separator<char> > tokenizer(exitsString, boost::char_separator<char>(", "));
        for(const string& edgeIdString: tokenizer) {
            uint64_t segmentId = invalid<uint64_t>;
            try {
                segmentId = std::stol(edgeIdString);
            } catch(exception&) {
            }
            if(segmentId == invalid<uint64_t>) {
                html << "Invalid segment " << edgeIdString << ". Must be a number.";
                return;
            }

            // Find the AssemblyGraphEdge corresponding to the requested segment.
            auto it = assemblyGraph.edgeMap.find(segmentId);
            if(it == assemblyGraph.edgeMap.end()) {
                html << "<p>Assembly graph at stage " << assemblyStage <<
                    " does not have segment " << segmentId;
                return;
            }

            exits.push_back(it->second);
        }
    }
    std::ranges::sort(exits, assemblyGraph.orderById);

    // Compute the tangle matrix.
    const TangleMatrix1 tangleMatrix(assemblyGraph, entrances, exits, html);
    GTest gTest(tangleMatrix.tangleMatrix, epsilon, false, false);
    gTest.writeHtml(html);



    // Create a RestrictedAnchorGraph for each element of the top hypothesis
    // that is set to 1.
    if(gTest.hypotheses.empty()) {
        return;
    }
    const auto& bestConnectivityMatrix = gTest.hypotheses.front().connectivityMatrix;
    for(uint64_t iEntrance=0; iEntrance<entrances.size(); iEntrance++) {
        for(uint64_t iExit=0; iExit<exits.size(); iExit++) {
            if(bestConnectivityMatrix[iEntrance][iExit]) {

                const AssemblyGraph::edge_descriptor eEntrance = entrances[iEntrance];
                const AssemblyGraph::edge_descriptor eExit = exits[iExit];

                const AssemblyGraphEdge& entranceEdge = assemblyGraph[eEntrance];
                const AssemblyGraphEdge& exitEdge = assemblyGraph[eExit];

                const AnchorId entranceAnchorId = entranceEdge.back().anchorPair.anchorIdB;
                const AnchorId exitAnchorId = exitEdge.front().anchorPair.anchorIdA;

                if(entranceAnchorId == exitAnchorId) {
                    html << "<br>The two anchors are coincident.";
                } else {

                    html << "<h4>RestrictedAnchorGraph to connect entrance " <<
                        entranceEdge.id <<
                        " with exit " << exitEdge.id << "</h4>"
                        "Last AnchorId on entrance is " << anchorIdToString(entranceAnchorId) <<
                        "<br>First AnchorId on exit is " << anchorIdToString(exitAnchorId);


                    RestrictedAnchorGraph restrictedAnchorGraph(
                        anchors(), journeys(), tangleMatrix, iEntrance, iExit, html);

                    html << "<br>The final RestrictedAnchorGraph has " << num_vertices(restrictedAnchorGraph) <<
                        " vertices and " << num_edges(restrictedAnchorGraph) << " edges ";

                    html << "<br>";
                    restrictedAnchorGraph.writeOrientedReadsInVertices(html);

                    // Find the optimal path in the RestrictedAnchorGraph.
                    vector<RestrictedAnchorGraph::edge_descriptor> longestPath;
                    restrictedAnchorGraph.findOptimalPath(entranceAnchorId, exitAnchorId, longestPath);

                    // Write it out in Graphviz format.
                    restrictedAnchorGraph.writeHtml(html, {entranceAnchorId, exitAnchorId});
                }

            }
        }

    }
}



void Assembler::exploreSegmentPair(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string assemblyStage;
    HttpServer::getParameterValue(request, "assemblyStage", assemblyStage);
    boost::trim(assemblyStage);

    uint64_t segmentId0 = invalid<uint64_t>;
    HttpServer::getParameterValue(request, "segmentId0", segmentId0);

    uint64_t segmentId1 = invalid<uint64_t>;
    HttpServer::getParameterValue(request, "segmentId1", segmentId1);



    // Start the form.
    html << "<form>";

    html <<
        "<table>"
        "<tr>"
        "<th class=left>Assembly stage"
        "<td class=centered><input type=text name=assemblyStage style='text-align:center' required";
    if(not assemblyStage.empty()) {
        html << " value='" << assemblyStage + "'";
    }
    html << " size=10>";

    html <<
        "<tr>"
        "<th class=left>SegmentId0"
        "<td class=centered>"
        "<input type=text name=segmentId0 style='text-align:center'";
    if(segmentId0 != invalid<uint64_t>) {
        html << " value='" << segmentId0 << "'";
    }
    html << " size=10>";

    html <<
        "<tr>"
        "<th class=left>SegmentId1"
        "<td class=centered>"
        "<input type=text name=segmentId1 style='text-align:center'";
    if(segmentId1 != invalid<uint64_t>) {
        html << " value='" << segmentId1 << "'";
    }
    html << " size=10>";


    // End the form.
    html <<
        "</table>"
        "<input type=submit value='Go'>"
        "</form>";

    if((segmentId0 == invalid<uint64_t>) or (segmentId1 == invalid<uint64_t>)) {
        return;
    }

    // Get the AssemblyGraph for this assembly stage.
    AssemblyGraphPostprocessor& assemblyGraph = getAssemblyGraph(
        assemblyStage,
        *httpServerData.options);



    // Find the edges corresponding to these ids.
    auto it0 = assemblyGraph.edgeMap.find(segmentId0);
    if(it0 == assemblyGraph.edgeMap.end()) {
        html << "<br>Assembly stage " << assemblyStage <<
            " does not have segment " << segmentId0;
        return;
    }
    const AssemblyGraph::edge_descriptor e0 = it0->second;
    // const AssemblyGraphEdge& edge0 = assemblyGraph[e0];

    auto it1 = assemblyGraph.edgeMap.find(segmentId1);
    if(it1 == assemblyGraph.edgeMap.end()) {
        html << "<br>Assembly stage " << assemblyStage <<
            " does not have segment " << segmentId1;
        return;
    }
    const AssemblyGraph::edge_descriptor e1 = it1->second;




    html << "<h2>Segment pair " << segmentId0 << " " << segmentId1 << "</h2>";

    const uint32_t representativeRegionStepCount =
        uint32_t(assemblyGraph.options.representativeRegionStepCount);
    SegmentStepSupport::analyzeSegmentPair(html, assemblyGraph, e0, e1, representativeRegionStepCount);
}



void Assembler::exploreSimilarSequences(const vector<string>& request, ostream& html)
{
    // Get the options from the request.
    string sequence0String;
    HttpServer::getParameterValue(request, "sequence0", sequence0String);
    boost::trim(sequence0String);

    string sequence1String;
    HttpServer::getParameterValue(request, "sequence1", sequence1String);
    boost::trim(sequence1String);

    // Write the form.
    html <<
        "<h2>Similar sequences</h2><form>"
        "<table>"
        "<tr>"
        "<th class=left>First sequence"
        "<td class=centered><input type=text name=sequence0 style='text-align:center;font-family:monospace' required size=100"
        " value='" << sequence0String << "'>"
        "<tr>"
        "<th class=left>Second sequence"
        "<td class=centered><input type=text name=sequence1 style='text-align:center;font-family:monospace' required size=100"
        " value='" << sequence1String << "'>"
        "</table>"
        "<input type=submit value='Analyze'>"
        "</form>";

    if(sequence0String.empty() or sequence1String.empty()) {
        return;
    }

    // Fill in the Base sequences.
    vector<Base> sequence0;
    for(const char c: sequence0String) {
        sequence0.push_back(Base::fromCharacter(c));
    }
    vector<Base> sequence1;
    for(const char c: sequence1String) {
        sequence1.push_back(Base::fromCharacter(c));
    }

    // EXPOSE WHEN CODE STABILIZES.
    const vector<uint64_t> minRepeatCount = {0, 2, 2, 2, 2, 2, 2};
    const bool areSimilar = areSimilarSequences(sequence0, sequence1, minRepeatCount, html);
    if(areSimilar) {
        html << "<p>These sequence are similar and their differences are likely caused by sequencing errors.";
    }
}

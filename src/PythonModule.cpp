#ifdef SHASTA_PYTHON_API

// Shasta.
#include "Assembler.hpp"
#include "AssemblerOptions.hpp"
#include "AssemblyGraphPostprocessor.hpp"
#include "Base.hpp"
#include "deduplicate.hpp"
#include "Detangler.hpp"
#include "diploidBayesianPhase.hpp"
#include "extractKmer128.hpp"
#include "globalMsa.hpp"
#include "LongBaseSequence.hpp"
#include "mappedCopy.hpp"
#include "MultithreadedObject.hpp"
#include "PermutationDetangler.hpp"
#include "performanceLog.hpp"
#include "ShortBaseSequence.hpp"
#include "splitRange.hpp"
#include "testSpoa.hpp"
#include "testSubsetGraph.hpp"
#include "TrivialDetangler.hpp"
using namespace shasta;

// Pybind11
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
using namespace pybind11;



PYBIND11_MODULE(shasta2, shasta2Module)
{
    // Class AssemblerOptions::AssemblyGraphOptions.
    class_<AssemblerOptions::AssemblyGraphOptions>(shasta2Module, "AssemblyGraphOptions")
        .def_readonly("transitiveReductionThreshold", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionThreshold)
        .def_readonly("transitiveReductionA", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionA)
        .def_readonly("transitiveReductionB", &AssemblerOptions::AssemblyGraphOptions::transitiveReductionB)
        .def_readonly("minCommonCoverage", &AssemblerOptions::AssemblyGraphOptions::minCommonCoverage)
        ;


    // Class AssemblerOptions.
    class_<AssemblerOptions>(shasta2Module, "AssemblerOptions")

        // Constructor from the name of a configuration file.
        .def(pybind11::init<const string&>(),
            arg("configurationFileName") = "shasta2.conf")
        .def_readonly("k", &AssemblerOptions::k)
        .def_readonly("markerDensity", &AssemblerOptions::markerDensity)
        .def_readonly("minAnchorCoverage", &AssemblerOptions::minAnchorCoverage)
        .def_readonly("maxAnchorCoverage", &AssemblerOptions::maxAnchorCoverage)
        .def_readonly("minAnchorGraphEdgeCoverage", &AssemblerOptions::minAnchorGraphEdgeCoverage)
        .def_readonly("assemblyGraphOptions", &AssemblerOptions::assemblyGraphOptions)
        ;

    // Class Assembler.
    class_<Assembler>(shasta2Module, "Assembler")

        // Constructor.
        .def(pybind11::init<const string&>(),
            "Assembler constructor.",
            arg("largeDataFileNamePrefix") = "Data/")

        // Reads
        .def("histogramReadLength",
            &Assembler::histogramReadLength,
            "Create a histogram of read length and write it to a csv file.",
            arg("fileName") = "ReadLengthHistogram.csv")

        // K-mer checker.
        .def("createKmerChecker",
            &Assembler::createKmerChecker,
            arg("k"),
            arg("markerDensity"),
            arg("threadCount") = 0)
        .def("accessKmerChecker",
            &Assembler::accessKmerChecker)

         // Markers.
        .def("accessMarkers",
            &Assembler::accessMarkers)
        .def("createMarkers",
            &Assembler::createMarkers,
            arg("threadCount") = 0)

         // Marker k-mers.
        .def("createMarkerKmers",
            &Assembler::createMarkerKmers,
            arg("threadCount") = 0)
        .def("accessMarkerKmers",
            &Assembler::accessMarkerKmers)

        // Anchors.
       .def("createAnchors",
           &Assembler::createAnchors,
           arg("minAnchorCoverage"),
           arg("maxAnchorCoverage"),
           arg("maxHomopolymerLength"),
           arg("threadCount") = 0)
       .def("accessAnchors",
           &Assembler::accessAnchors,
           arg("writeAccess") = false)

       // Journeys.
      .def("createJourneys",
          &Assembler::createJourneys,
          arg("threadCount") = 0)
      .def("accessJourneys",
          &Assembler::accessJourneys)

      // AssemblyGraph.
     .def("createAssemblyGraph",
         &Assembler::createAssemblyGraph,
         arg("assemblerOptions"),
         arg("threadCount") = 0)
     .def("getAssemblyGraph",
         &Assembler::getAssemblyGraph)
    ;



    // AssemblyGraph and related access functions.
    class_<AssemblyGraph::vertex_descriptor>(shasta2Module, "AssemblyGraphVertexDescriptor")
        ;
    class_<AssemblyGraph::edge_descriptor>(shasta2Module, "AssemblyGraphEdgeDescriptor")
        ;
    class_<AssemblyGraphVertex>(shasta2Module, "AssemblyGraphVertex")
        .def_readwrite("anchorId", &AssemblyGraphVertex::anchorId)
        ;
    class_<AssemblyGraphEdge>(shasta2Module, "AssemblyGraphEdge")
        .def_readwrite("id", &AssemblyGraphEdge::id)
        ;
    class_<AssemblyGraphPostprocessor>(shasta2Module, "AssemblyGraph")
        .def("getVertexDescriptor", &AssemblyGraphPostprocessor::getVertexDescriptor)
        .def("getVertex", &AssemblyGraphPostprocessor::getVertex, return_value_policy::reference)
        .def("getEdgeDescriptor", &AssemblyGraphPostprocessor::getEdgeDescriptor)
        .def("getEdge", &AssemblyGraphPostprocessor::getEdge, return_value_policy::reference)
        .def("transitiveReduction", &AssemblyGraphPostprocessor::transitiveReduction)
        .def("detangleVertices",
            (
                void (AssemblyGraphPostprocessor::*)
                (TrivialDetangler&)
            )
            &AssemblyGraphPostprocessor::detangleVertices)
        .def("detangleEdges",
            (
                void (AssemblyGraphPostprocessor::*)
                (TrivialDetangler&)
            )
            &AssemblyGraphPostprocessor::detangleEdges)
        .def("detangleVertices",
            (
                void (AssemblyGraphPostprocessor::*)
                (PermutationDetangler&)
            )
            &AssemblyGraphPostprocessor::detangleVertices)
        .def("detangleEdges",
            (
                void (AssemblyGraphPostprocessor::*)
                (PermutationDetangler&)
            )
            &AssemblyGraphPostprocessor::detangleEdges)
        ;
    class_<TrivialDetangler>(shasta2Module, "TrivialDetangler")
        .def(init<uint64_t>())
        ;
    class_<PermutationDetangler>(shasta2Module, "PermutationDetangler")
        .def(init<uint64_t>())
        ;


    // Non-member functions exposed to Python.
    shasta2Module.def("openPerformanceLog",
        openPerformanceLog
        );
    shasta2Module.def("testMultithreadedObject",
        testMultithreadedObject
        );
    shasta2Module.def("testMemoryMappedVector",
        testMemoryMappedVector
        );
    shasta2Module.def("testBase",
        testBase
        );
    shasta2Module.def("testShortBaseSequence",
        testShortBaseSequence
        );
    shasta2Module.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shasta2Module.def("testExtractKmer128",
        testExtractKmer128
        );
    shasta2Module.def("testSplitRange",
        testSplitRange
        );
    shasta2Module.def("testSpoa",
        testSpoa
        );
    shasta2Module.def("testDeduplicateAndCount",
        testDeduplicateAndCount
        );
    shasta2Module.def("testBitReversal",
        testBitReversal
        );
    shasta2Module.def("mappedCopy",
        mappedCopy
        );
    shasta2Module.def("testLongBaseSequence",
        testLongBaseSequence
        );
    shasta2Module.def("testDiploidBayesianPhase",
        testDiploidBayesianPhase
        );
    shasta2Module.def("testSubsetGraph",
        testSubsetGraph
        );
    shasta2Module.def("globalMsaPython",
        globalMsaPython
        );
}

#endif

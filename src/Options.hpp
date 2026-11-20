#pragma once

// Shasta.
#include "cstdint.hpp"
#include "iosfwd.hpp"
#include "string.hpp"
#include "vector.hpp"

// CLI11. Requires package libcli11-dev.
#include "CLI/CLI.hpp"

namespace shasta2 {
    class Options;
}


class shasta2::Options : public CLI::App {
public:

    // Constructor from command line options.
    Options(int argumentCount, char** arguments);

    // Constructor from a configuration file.
    // Used by the Python API.
    Options(const string& configurationFileName = "");

    // Construct with just default values. 
    // Used by createBashCompletionScript.
    Options(int);

    // Write a configuration file.
    void write(ostream&) const;

    bool isHelp = false;

    // Default values of vector options.
    const vector<string> inputFileNamesDefault;
    const vector<uint64_t> maxAnchorRepeatLengthDefault = {6, 4, 4, 4, 4};
    const vector<uint64_t> minAnchorDistinctSubkmerCountDefault = {4, 12, 24};

	// Options with a simple type are defined in OptionsDefine.hpp
	#define SHASTA2_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
		type name = defaultValue;
    #define SHASTA2_VECTOR_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
        vector<type> name = defaultValue;
    #define SHASTA2_BOOL_OPTION_DEFINE(name, optionName, defaultValue, description) \
        bool name = defaultValue;
	#include "OptionsDefine.hpp"
    #undef SHASTA2_OPTION_DEFINE
    #undef SHASTA2_VECTOR_OPTION_DEFINE
    #undef SHASTA2_BOOL_OPTION_DEFINE

    // Return the actual number of threads to be used.
    uint64_t actualThreadCount() const;

private:
    void addOptions();
    void validate();
    void load(istream&);
};



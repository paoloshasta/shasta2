// Shasta.
#include "Options.hpp"
#include "configurationTable.hpp"
using namespace shasta2;

// Standard library.
#include "array.hpp"
#include "fstream.hpp"
#include "iostream.hpp"
#include "stdexcept.hpp"
#include <thread>



Options::Options(int argc, char** argv) :
    CLI::App("Shasta2. Under development.")
{
    allow_config_extras(false);
    get_formatter()->column_width(10);

    addOptions();

    // First get the options from the arguments.
    try {
        parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        exit(e);
        isHelp = true;
        return;
    }

    if((command != "assemble") and (command != "explore")) {
        return;
    }

    // When initializing from command line and --command is assemble or explore, --config is mandatory.
    if(config.empty()) {
        throw runtime_error("Option --config is missing. "
            "It is mandatory when running an assembly "
            "and should specify a built-in assembly configuration "
            "or a configuration file.");
    }

    // Now load from the specified configuration.
    // Values that were already specified on the command line are not overwritten.
    const auto it = configurationTable.find(config);
    if(it == configurationTable.end()) {

        // From configuration file.
        ifstream configurationFile(config);
        if(not configurationFile) {
            cout << "Option --config does not specify a built-in assembly configuration "
                "and is not a valid configuration file." << endl;
            cout << "Built-in assembly configurations are:" << endl;
            for(const auto&[name, _]: configurationTable) {
                cout << name << endl;
            }
            throw runtime_error("Invalid --config option. See above for a list of built-in assembly configurations.");

        }
        load(configurationFile);

    } else {

        // From the configuration table.
        const string& configString = it->second;
        std::istringstream s(configString);
        load(s);
    }

    validate();
}



void Options::addOptions()
{
	#define SHASTA2_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
		add_option(optionName, name, description)->capture_default_str();
    #define SHASTA2_VECTOR_OPTION_DEFINE(type, name, optionName, defaultValue, description) \
        add_option(optionName, name, description)->capture_default_str();
    #define SHASTA2_BOOL_OPTION_DEFINE(name, optionName, defaultValue, description) \
        add_flag(optionName, name, description)->capture_default_str();

	#include "OptionsDefine.hpp"

	#undef SHASTA2_OPTION_DEFINE
    #undef SHASTA2_VECTOR_OPTION_DEFINE
    #undef SHASTA2_BOOL_OPTION_DEFINE
}



// Constructor from a configuration file.
Options::Options(const string& configurationFileName)
{

    allow_config_extras(false);
    addOptions();

    config = configurationFileName;
    ifstream s(config);
    if(not s) {
        throw runtime_error("Could not open configuration file " + config);
    }
    load(s); 

    validate();
}



// Construct with just default values. 
// Used by createBashCompletionScript.
Options::Options(int)
{
    addOptions();
}



void Options::load(istream& s) 
{
    parse_from_stream(s);
}



void Options::write(ostream& s) const
{
    s << config_to_str(true,true);
}



void Options::validate()
{
    // Adjust the numbers of threads, if necessary.
    if(threadCount == 0) {
        threadCount = std::thread::hardware_concurrency();
    }
}



// Return the actual number of threads to be used.
uint64_t Options::actualThreadCount() const
{
    if(threadCount == 0) {
        return std::thread::hardware_concurrency();
    } else {
        return threadCount;
    }
}

// Shasta.
#include "Options.hpp"
using namespace shasta2;

// Standard library.
#include "array.hpp"
#include "stdexcept.hpp"
#include <thread>



Options::Options(int argc, char** argv) :
    CLI::App("Shasta2. Under development.")
{
    allow_config_extras(false);
    set_config("--config", "", "Specify a configuration file.");

    get_formatter()->column_width(10);

    addOptions();

    try {
        parse(argc, argv);
    } catch(const CLI::ParseError& e) {
         exit(e);
         if(e.get_name() == "CallForHelp") {
             isHelp = true;
         } else {
             throw runtime_error("Error parsing options.");
         }
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
Options::Options(const string& fileName)
{

    allow_config_extras(false);
    set_config("--config", "", "Configuration file.");
    get_formatter()->column_width(20);

    addOptions();

    if(not fileName.empty()) {
        // Construct arguments "--config fileName".
        int argc = 3;
        const string name = "shasta2";
        const string keyword = "--config";
        char* arg0 = const_cast<char*>(name.c_str());
        char* arg1 = const_cast<char*>(keyword.c_str());
        char* arg2 = const_cast<char*>(fileName.c_str());
        const array<char*, 3> argv = {arg0, arg1, arg2};

        try {
            parse(argc, &argv.front());
        } catch(const CLI::ParseError& e) {
             exit(e);
             throw runtime_error("Error parsing options.");
        }
    }

    validate();
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


// Shasta.
#include "Assembler.hpp"
#include "Options.hpp"
#include "filesystem.hpp"
#include "HomopolymerModelTable.hpp"
#include "html.hpp"
#include "htmlLog.hpp"
#include "memoryInformation.hpp"
#include "performanceLog.hpp"
#include "ReadSummary.hpp"
#include "Tee.hpp"
#include "timestamp.hpp"
using namespace shasta2;

// Boost libraries.
#include <boost/program_options.hpp>
#include  <boost/chrono/process_cpu_clocks.hpp>

//  Linux.
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>


// Standard library.
#include "chrono.hpp"
#include <filesystem>
#include "iostream.hpp"
#include "iterator.hpp"
#include "stdexcept.hpp"



namespace shasta2 {
    namespace main {

        void main(int argumentCount, char** arguments);

        void setupRunDirectory(
            const string& memoryMode,
            const string& memoryBacking,
            size_t& pageSize,
            string& dataDirectory
            );

        void setupHugePages();
        void segmentFaultHandler(int);
        void setupSegmentFaultHandler();
        void finalMemoryCleanup(const Options&);

        // Functions that implement --command keywords
        void assemble(const Options&, int argumentCount, char** arguments);
        void saveBinaryData(const Options&);
        void cleanupBinaryData(const Options&);
        void explore(const Options&);
        void listCommands();
        void listHomopolymerModels();
        void createBashCompletionScript();

        const std::set<string> commands = {
            "assemble",
            "saveBinaryData",
            "cleanupBinaryData",
            "explore",
            "listCommands",
            "listHomopolymerModels",
            "createBashCompletionScript"};

        void writeHtmlLogInitialOutput(
            int argumentCount,
            char** arguments,
            uint64_t threadCount);
        void writeHtmlLogFinalOutput(double elapsedTime, double cpuTime, uint64_t threadCount);

    }

    // This is used to duplicate cout output to stdout.log.
    Tee tee;
    ofstream shastaLog;
}



int main(int argumentCount, char** arguments)
{
    try {
        shasta2::main::main(argumentCount, arguments);
    }

    catch(const boost::program_options::error_with_option_name& e) {
        cout << "Invalid option: " << e.what() << endl;
        return 1;
    }

    catch (const runtime_error& e) {
        cout << timestamp << e.what() << endl;
        return 2;
    }

    catch (const std::bad_alloc& e) {
        cout << timestamp << e.what() << endl;
        cout << "Memory allocation failure." << endl;
        cout << "This assembly requires more memory than available." << endl;
        cout << "Rerun on a larger machine." << endl;
        return 2;
    }

    catch (const exception& e) {
        cout << timestamp << e.what() << endl;
        return 3;
    }

    catch (...) {
        cout << timestamp << "Terminated after catching a non-standard exception." << endl;
        return 4;
    }

    return 0;
}



void shasta2::main::segmentFaultHandler(int)
{
    char message[] = "\nA segment fault occurred. Please report it by filing an "
        "issue on the Shasta2 repository and attaching the entire log output. "
        "To file an issue, point your browser to https://github.com/paoloshasta/shasta2/issues\n";
    ::write(fileno(stderr), message, sizeof(message));
    ::_exit(1);
}



void shasta2::main::setupSegmentFaultHandler()
{
    struct sigaction action;
    ::memset(&action, 0, sizeof(action));
    action.sa_handler = &segmentFaultHandler;
    sigaction(SIGSEGV, &action, 0);
}


#ifndef SHASTA2_STATIC_LIBRARY
void shasta2::main::main(int argumentCount, char** arguments)
{
    setupSegmentFaultHandler();

    // Parse command line options and the configuration file, if one was specified.
    Options options(argumentCount, arguments);
    if(options.isHelp) {
        return;
    }


    // Check that we have a valid command.
    auto it = commands.find(options.command);
    if(it ==commands.end()) {
        const string message = "Invalid command " + options.command;
        listCommands();
        throw runtime_error(message);
    }



    // Execute the requested command.
    if(options.command == "assemble") {
        assemble(options, argumentCount, arguments);
        return;
    } else if(options.command == "cleanupBinaryData") {
        cleanupBinaryData(options);
        return;
    } else if(options.command == "saveBinaryData") {
        saveBinaryData(options);
        return;
    } else if(options.command == "explore") {
        explore(options);
        return;
    } else if(options.command == "listCommands") {
        listCommands();
        return;
    } else if(options.command == "listHomopolymerModels") {
        listHomopolymerModels();
        return;
    } else if(options.command == "createBashCompletionScript") {
        createBashCompletionScript();
        return;
    }

    // We already checked for a valid command above, so if we get here
    // the above logic is missing code for one of the valid commands.
    SHASTA2_ASSERT(0);

}
#endif



// Implementation of --command assemble.
void shasta2::main::assemble(
    const Options& options,
    int argumentCount, char** arguments)
{
    SHASTA2_ASSERT(options.command == "assemble");

    const auto steadyClock0 = std::chrono::steady_clock::now();
    const auto userClock0 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock0 = boost::chrono::process_system_cpu_clock::now();


    // Various checks for option validity.

    // Check --k.
    if(options.k > 126) {
        throw runtime_error("Invalid value specified for --k. "
            "Must be no more than 126");
    }
    if((options.k % 2) == 1) {
        throw runtime_error("Invalid value specified for --k. Must be even.");
    }

    // Check that we have at least one input file.
    if(options.inputFileNames.empty()) {
        throw runtime_error("Specify at least one input file "
            "using command line option --input.");
    }



    // Find absolute paths of the input files.
    // We will use them below after changing directory to the output directory.
    vector<string> inputFileAbsolutePaths;
    for(const string& inputFileName: options.inputFileNames) {
        if(!std::filesystem::exists(inputFileName)) {
            throw runtime_error("Input file not found: " + inputFileName);
        }
        if(!std::filesystem::is_regular_file(inputFileName)) {
            throw runtime_error("Input file is not a regular file: " + inputFileName);
        }
        inputFileAbsolutePaths.push_back(filesystem::getAbsolutePath(inputFileName));
    }
    string externalAnchorsNameAbsolutePath;
    if(not options.externalAnchorsName.empty()) {
        externalAnchorsNameAbsolutePath = filesystem::getAbsolutePath(options.externalAnchorsName);
    }
    string externalAnchorGraphNameAbsolutePath;
    if(not options.externalAnchorGraphName.empty()) {
        externalAnchorGraphNameAbsolutePath = filesystem::getAbsolutePath(options.externalAnchorGraphName);
    }

    // Check --homopolymer-model. It must be in the homopolymerModelTable,
    // or else it must be an absolute path.
    if(not options.homopolymerModelName.empty()) {
        if(not options.useMsaRepair) {
            throw runtime_error("--homopolymer-model requires --use-msa-repair.");
        }
        if(not homopolymerModelTable.contains(options.homopolymerModelName)) {
            if(options.homopolymerModelName[0] != '/') {
                cout << "Invalid homopolymer model. Option --homopolymer-model must specify "
                    "a built-in homopolymer model or an absolute path: " +
                    options.homopolymerModelName << endl;
                listHomopolymerModels();
                throw runtime_error("Invalid homopolymer model. "
                        "Option --homopolymer-model must specify "
                        "a built-in homopolymer model or an absolute path: " + options.homopolymerModelName);
            }
            if(!std::filesystem::is_regular_file(options.homopolymerModelName)) {
                throw runtime_error("Homopolymer model not found: " + options.homopolymerModelName);
            }
        }
    }



    // Create the assembly directory. If it exists, stop.
    bool exists = std::filesystem::exists(options.assemblyDirectory);
    if (exists) {
        throw runtime_error(
            options.assemblyDirectory +
            " already exists. Remove it first \n"
            "or use --assemblyDirectory to specify a different assembly directory."
        );
    } else {
        SHASTA2_ASSERT(std::filesystem::create_directory(options.assemblyDirectory));
    }

    // Make the assembly directory current.
    std::filesystem::current_path(options.assemblyDirectory);

    // Write a configuration to the assembly directory.
    {
        ofstream configurationFile("shasta2.conf");
        options.write(configurationFile);
    }

    // Open the performance log.
    openPerformanceLog("performance.log");
    performanceLog << timestamp << "Assembly begins." << endl;

    // Open the html log.
    openHtmlLog("AssemblyLog.html");
    writeHtmlBegin(htmlLog, "Shasta2 assembly");
    writeHtmlLogInitialOutput(argumentCount, arguments, options.threadCount);

    // Open stdout.log and "tee" (duplicate) stdout to it.
    shastaLog.open("stdout.log");
    tee.duplicate(cout, shastaLog);

    // Echo out the command line options.
    cout << timestamp << "Assembly begins.\nCommand line:" << endl;
    for(int i=0; i<argumentCount; i++) {
        cout << arguments[i] << " ";
    }
    cout << endl;

    // Set up the run directory as required by the memoryMode and memoryBacking options.
    size_t pageSize = 0;
    string dataDirectory;
    setupRunDirectory(
        options.memoryMode,
        options.memoryBacking,
        pageSize,
        dataDirectory);


    {
        // Create the Assembler.
        Assembler assembler(dataDirectory, pageSize);

        // Run the assembly.
        assembler.assemble(
            options,
            inputFileAbsolutePaths,
            externalAnchorsNameAbsolutePath,
            externalAnchorGraphNameAbsolutePath);
    }

    // Final memory cleanup must happen after the Assembler is destroyed.
    finalMemoryCleanup(options);

    // Write performance information.
    const auto steadyClock1 = std::chrono::steady_clock::now();
    const auto userClock1 = boost::chrono::process_user_cpu_clock::now();
    const auto systemClock1 = boost::chrono::process_system_cpu_clock::now();
    const double elapsedTime = 1.e-9 * double((
        std::chrono::duration_cast<std::chrono::nanoseconds>(steadyClock1 - steadyClock0)).count());
    const double userCpuTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(userClock1 - userClock0)).count());
    const double systemCpuTime = 1.e-9 * double((
        boost::chrono::duration_cast<boost::chrono::nanoseconds>(systemClock1 - systemClock0)).count());
    const double cpuTime = userCpuTime + systemCpuTime;
    writeHtmlLogFinalOutput(elapsedTime, cpuTime, options.threadCount);
    writeHtmlEnd(htmlLog);

    cout << timestamp << "Assembly ends." << endl;
    performanceLog << timestamp << "Assembly ends." << endl;
}



// Set up the run directory as required by the memoryMode and memoryBacking options.
void shasta2::main::setupRunDirectory(
    const string& memoryMode,
    const string& memoryBacking,
    size_t& pageSize,
    string& dataDirectory
    )
{

    if(memoryMode == "anonymous") {

        if(memoryBacking == "disk") {

            // This combination is meaningless.
            throw runtime_error("\"--memoryMode anonymous\" is not allowed in combination "
                "with \"--memoryBacking disk\".");

        } else if(memoryBacking == "4K") {

            // Anonymous memory on 4KB pages.
            // This combination is the default.
            // It does not require root privilege.
            dataDirectory = "";
            pageSize = 4096;

        } else if(memoryBacking == "2M") {

            // Anonymous memory on 2MB pages.
            // This may require root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            // Root privilege is not required if 2M pages have already
            // been set up as required.
            setupHugePages();
            pageSize = 2 * 1024 * 1024;

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else if(memoryMode == "filesystem") {

        if(memoryBacking == "disk") {

            // Binary files on disk.
            // This does not require root privilege.
            SHASTA2_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;

        } else if(memoryBacking == "4K") {

            // Binary files on the tmpfs filesystem
            // (filesystem in memory backed by 4K pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            SHASTA2_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 4096;
            const string command = "sudo mount -t tmpfs -o size=0 tmpfs Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else if(memoryBacking == "2M") {

            // Binary files on the hugetlbfs filesystem
            // (filesystem in memory backed by 2M pages).
            // This requires root privilege, which is obtained using sudo
            // and may result in a password prompting depending on sudo set up.
            setupHugePages();
            SHASTA2_ASSERT(std::filesystem::create_directory("Data"));
            dataDirectory = "Data/";
            pageSize = 2 * 1024 * 1024;
            const uid_t userId = ::getuid();
            const gid_t groupId = ::getgid();
            const string command = "sudo mount -t hugetlbfs -o pagesize=2M"
                ",uid=" + to_string(userId) +
                ",gid=" + to_string(groupId) +
                " none Data";
            const int errorCode = ::system(command.c_str());
            if(errorCode != 0) {
                throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
                    " running command: " + command);
            }

        } else {
            throw runtime_error("Invalid value specified for --memoryBacking: " + memoryBacking +
                "\nValid values are: disk, 4K, 2M.");
        }

    } else {
        throw runtime_error("Invalid value specified for --memoryMode: " + memoryMode +
            "\nValid values are: anonymous, filesystem.");
    }
}



// This function sets nr_overcommit_hugepages for 2MB pages
// to a little below total memory.
// If the setting needs to be modified, it acquires
// root privilege via sudo. This may result in the
// user having to enter a password.
void shasta2::main::setupHugePages()
{

    // Get the total memory size.
    const uint64_t totalMemoryBytes = sysconf(_SC_PAGESIZE) * sysconf(_SC_PHYS_PAGES);

    // Figure out how much memory we want to allow for 2MB pages.
    const uint64_t MB = 1024 * 1024;
    const uint64_t GB = MB * 1024;
    const uint64_t maximumHugePageMemoryBytes = totalMemoryBytes - 8 * GB;
    const uint64_t maximumHugePageMemoryHugePages = maximumHugePageMemoryBytes / (2 * MB);

    // Check what we have it set to.
    const string fileName = "/sys/kernel/mm/hugepages/hugepages-2048kB/nr_overcommit_hugepages";
    ifstream file(fileName);
    if(!file) {
        throw runtime_error("Error opening " + fileName + " for read.");
    }
    uint64_t currentValue = 0;
    file >> currentValue;
    file.close();

    // If it's set to at least what we want, don't do anything.
    // When this happens, root access is not required.
    if(currentValue >= maximumHugePageMemoryHugePages) {
        return;
    }

    // Use sudo to set.
    const string command =
        "sudo sh -c \"echo " +
        to_string(maximumHugePageMemoryHugePages) +
        " > " + fileName + "\"";
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command: " + command);
    }

}



// Implementation of --command saveBinaryData.
// This copies Data to DataOnDisk.
void shasta2::main::saveBinaryData(
    const Options& options)
{
    SHASTA2_ASSERT(options.command == "saveBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        options.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        throw runtime_error(dataDirectory + " does not exist, nothing done.");
    }

    // Check that the DataOnDisk directory does not exist.
    const string dataOnDiskDirectory =
        options.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        throw runtime_error(dataOnDiskDirectory + " already exists, nothing done.");
    }

    // Copy Data to DataOnDisk.
    const string command = "cp -rp " + dataDirectory + " " + dataOnDiskDirectory;
    const int errorCode = ::system(command.c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " running command:\n" + command);
    }
    cout << "Binary data successfully saved." << endl;
}



// Implementation of --command cleanupBinaryData.
void shasta2::main::cleanupBinaryData(
    const Options& options)
{
    SHASTA2_ASSERT(options.command == "cleanupBinaryData");

    // Locate the Data directory.
    const string dataDirectory =
        options.assemblyDirectory + "/Data";
    if(!std::filesystem::exists(dataDirectory)) {
        cout << dataDirectory << " does not exist, nothing done." << endl;
        return;
    }

    // Unmount it and remove it.
    ::system(("sudo umount " + dataDirectory).c_str());
    const int errorCode = ::system(string("rm -rf " + dataDirectory).c_str());
    if(errorCode != 0) {
        throw runtime_error("Error " + to_string(errorCode) + ": " + strerror(errorCode) +
            " removing " + dataDirectory);
    }
    cout << "Cleanup of " << dataDirectory << " successful." << endl;

    // If the DataOnDisk directory exists, create a symbolic link
    // Data->DataOnDisk.
    const string dataOnDiskDirectory =
        options.assemblyDirectory + "/DataOnDisk";
    if(std::filesystem::exists(dataOnDiskDirectory)) {
        std::filesystem::current_path(options.assemblyDirectory);
        const string command = "ln -s DataOnDisk Data";
        ::system(command.c_str());
    }

}

// Implementation of --command explore.
void shasta2::main::explore(
    const Options& options)
{

    // Go to the assembly directory.
    std::filesystem::current_path(options.assemblyDirectory);

    // Check that we have the binary data.
    if(!std::filesystem::exists("Data")) {
        throw runtime_error("Binary directory \"Data\" not available "
        " in assembly directory " +
        options.assemblyDirectory +
        ". Use \"--memoryMode filesystem\", possibly followed by "
        "\"--command saveBinaryData\" and \"--command cleanupBinaryData\" "
        "if you want to make sure the binary data are persistently available on disk."
        );
        return;
    }

    // Create the Assembler.
    Assembler assembler("Data/");
    assembler.fillServerFunctionTable();


    // Access all available binary data.
    assembler.httpServerData.options = &options;
    assembler.accessAllSoft();


    // Start the http server.
    bool localOnly;
    bool sameUserOnly;
    if(options.exploreAccess == "user") {
        localOnly = true;
        sameUserOnly = true;
    } else if(options.exploreAccess == "local") {
        localOnly = true;
        sameUserOnly = false;
    } else if (options.exploreAccess == "unrestricted"){
        localOnly = false;
        sameUserOnly = false;
    } else {
        throw runtime_error("Invalid value specified for --exploreAccess. "
            "Only use this option if you understand its security implications."
        );
    }
    assembler.explore(
        options.port,
        localOnly,
        sameUserOnly);
}



void shasta2::main::listCommands()
{
    cout << "Valid commands are:" << endl;
    for(const string& command: commands) {
        cout << command << endl;
    }
}



void shasta2::main::listHomopolymerModels()
{
    cout << "Valid built-in homopoymer models are:" << endl;
    for(const auto&[name, ignore]: homopolymerModelTable) {
        cout << name << endl;
    }
}



void shasta2::main::createBashCompletionScript()
{
    const string fileName = "shasta2Completion.sh";
    ofstream file(fileName);

    file << "#!/bin/bash\n";
    file << "complete -o default -W \"\\\n";

    // Options.
    const Options options;
    for(const CLI::Option* option: options.get_options()) {
        for(const string& name: option->get_lnames()) {
            file << "--" << name << " \\\n";
        }
    }

    // Commands.
    for(const auto& command: commands) {
        file << command << " \\\n";
    }

    // Other keywords.
    file << "filesystem anonymous \\\n";
    file << "disk 4K 2M \\\n";
    file << "user local unrestricted \\\n";

    // Finish the "complete" command.
    file << "\" shasta2\n";

    cout << "Use \"source shasta2Completion.sh.\" for automatic completion "
        "of shasta2 commands in the bash shell." << endl;
}



void shasta2::main::writeHtmlLogInitialOutput(
    int argumentCount,
    char** arguments,
    uint64_t threadCount)
{
    const uint64_t virtualProcessorCount = std::thread::hardware_concurrency();
    if(threadCount == 0) {
        threadCount = virtualProcessorCount;
    }
    const uint64_t memoryBytes = getTotalPhysicalMemory();
    const double memoryGiB = double(memoryBytes) / (1024. * 1024. * 1024);

    htmlLog <<
        "<h1>Shasta2 assembly</h1>\n"
        "<table>\n"
        "<tr><th class=left>Begin time<td>" << timestamp;

    htmlLog << "<tr><th class=left>Command<td>";
    for(int i=0; i<argumentCount; i++) {
        htmlLog << arguments[i];
        if(i != argumentCount - 1) {
            htmlLog << " ";
        }
    }
    htmlLog << "\n";

    htmlLog <<
        "<tr><th class=left>Assembly directory<td>" <<
        std::filesystem::current_path().string() << "\n" <<
        "<tr><th class=left>Number of threads<td>" << threadCount << "\n" <<
        "<tr><th class=left>Number of virtual processors<td>" << virtualProcessorCount << "\n" <<
        "<tr><th class=left>Total usable physical memory<td>" << memoryBytes <<
        " bytes, " << memoryGiB << " GiB\n";

    htmlLog << "</table>\n" << flush;
}



void shasta2::main::writeHtmlLogFinalOutput(
    double elapsedTime,
    double cpuTime,
    uint64_t threadCount)
{
    const uint64_t virtualProcessorCount = std::thread::hardware_concurrency();
    if(threadCount == 0) {
        threadCount = virtualProcessorCount;
    }

    const uint64_t memoryBytes = getPeakMemoryUsage();
    const double memoryGiB = double(memoryBytes) / (1024. * 1024. * 1024);

    const double averageConcurrency =
        cpuTime / elapsedTime;

    const auto oldPrecision = htmlLog.precision(2);
    const auto oldFlags = htmlLog.setf(std::ios_base::fixed, std::ios_base::floatfield);

    htmlLog <<
        "<h2>Performance summary</h2>\n"
        "<table>\n"
        "<tr><th class=left>End time<td>" << timestamp <<
        "<tr><th class=left>Elapsed time<td>" << elapsedTime << " s = " <<
        elapsedTime / 60. << " m = " << elapsedTime / 3600. << " h\n"
        "<tr><th class=left>CPU time<td>" << cpuTime << " s = " <<
        cpuTime / 60. << " m = " << cpuTime / 3600. << " h\n"
        "<tr><th class=left>Average concurrency<td>" << averageConcurrency << "\n" <<

        "<tr><th class=left>Average CPU utilization<br>(relative to " <<
        threadCount << " threads)<td>" <<
        averageConcurrency/double(threadCount) << "\n" <<

        "<tr><th class=left>Average CPU utilization<br>(relative to " <<
        virtualProcessorCount << " virtual processors)<td>" <<
        averageConcurrency/double(virtualProcessorCount) << "\n" <<

        "<tr><th class=left>Peak virtual memory<td>" << memoryBytes << " bytes = " <<
        double(memoryBytes) / 1.e9 << " GB = " <<
        memoryGiB << " GiB\n" <<
        "</table>\n";

    htmlLog.precision(oldPrecision);
    htmlLog.flags(oldFlags);
}



// If --memory-mode filesystem and --keep-binary-data
// is not used, get rid of the binary data directory.
// This must be called after the Assembler destructor.
void shasta2::main::finalMemoryCleanup(const Options& options)
{
    // If --memory-mode anonymous, the process cleanup will take care of it.
    if(options.memoryMode == "anonymous") {
        return;
    }

    // --memory-mode filesystem
    SHASTA2_ASSERT(options.memoryMode == "filesystem");

    // If --keep-binary-data, don't do anything.
    if(options.keepBinaryData) {
        return;
    }

    // Remove all files in the Data directory.
    vector<string> fileNames;
    for(const auto& directoryEntry : std::filesystem::directory_iterator("Data")) {
        fileNames.push_back(string(directoryEntry.path()));
    }
    for(const string& fileName: fileNames) {
        try {
            std::filesystem::remove(fileName);
        } catch(...) {
            // Ignore.
        }
    }

    // Remove the Data directory, if --memory-backing disk.
    // Don't do anything if --memory-backing 2M, to avoid complications.
    if(options.memoryBacking == "disk") {
        std::filesystem::remove_all("Data");
    }
}

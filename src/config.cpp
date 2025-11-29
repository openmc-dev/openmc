#include "openmc/config.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include <dirent.h>

namespace openmc {

std::string list_files(const std::string& path) {
    std::string result = "";
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir(path.c_str())) != NULL) {
        while ((ent = readdir(dir)) != NULL) {
            if (std::string(ent->d_name) != "." && std::string(ent->d_name) != "..") {
                result += ent->d_name;
                result += " ";
            }
        }
        closedir(dir);
    }
    return result;
}

void print_config_usage() {
    std::cout << "Usage: openmc config [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -p, --prefix         Print the installation prefix" << std::endl;
    std::cout << "  -i, --include-dir    Print the include directory" << std::endl;
    std::cout << "  -l, --lib-dir        Print the library directory" << std::endl;
    std::cout << "  -c, --cmake-dir      Print the CMake configuration directory" << std::endl;
    std::cout << "  -L, --libs           Print the all libraries in the library directory" << std::endl;
    std::cout << "  -e, --extra-lib-dir  Print the extra library directory" << std::endl;
    std::cout << "  -E, --extra-libs     Print the extra libraries" << std::endl;
}

std::string exec(const std::string& cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd.c_str(), "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    // Remove trailing newline character if present
    if (!result.empty() && result.back() == '\n') {
        result.pop_back();
    }
    return result;
}

int openmc_config_main(int argc, char** argv) {
    std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty()) {
        print_config_usage();
        return 1;
    }

    // Default values
    std::string install_prefix_val;
    std::string include_dir_val;
    std::string lib_dir_val;
    std::string cmake_dir_val;
    std::string libs_val;
    std::string extra_libs_val;
    std::string extra_lib_dir_val;

#ifdef SKBUILD
    try {
        install_prefix_val = exec("python3 -c \"import sys; print(sys.prefix)\"");
        include_dir_val = exec("python3 -c \"import openmc.paths; print(''.join(openmc.paths.include_path))\"");
        lib_dir_val = exec("python3 -c \"import openmc.paths; print(''.join(openmc.paths.lib_path))\"");
        cmake_dir_val = exec("python3 -c \"import openmc.paths; print(openmc.paths.cmake_path)\" ");
        libs_val = exec("python3 -c \"import openmc.paths; print(' '.join(openmc.paths.lib))\"");
        extra_libs_val = exec("python3 -c \"import openmc.paths; print(' '.join(openmc.paths.extra_lib))\"");
        extra_lib_dir_val = exec("python3 -c \"import openmc.paths; print(openmc.paths.extra_lib_path)\"");
    } catch (const std::runtime_error& e) {
        std::cerr << "Fatal error: Could not get configuration from openmc.paths via Python: " << e.what() << std::endl;
        return 1;
    }
#else
    // Hardcoded values for CMake-only build
    install_prefix_val = CMAKE_INSTALL_PREFIX;
    include_dir_val = std::string(CMAKE_INSTALL_PREFIX) + "/" + CMAKE_INSTALL_INCLUDEDIR;
    lib_dir_val = std::string(CMAKE_INSTALL_PREFIX) + "/" + CMAKE_INSTALL_LIBDIR;
    cmake_dir_val = std::string(CMAKE_INSTALL_PREFIX) + "/" + CMAKE_INSTALL_LIBDIR + "/cmake/OpenMC";
    libs_val = list_files(lib_dir_val);
    extra_libs_val = "";
    extra_lib_dir_val = "";
#endif

    for (const auto& arg : args) {
        if (arg == "--prefix" || arg == "-p") {
            std::cout << install_prefix_val << std::endl;
            return 0;
        } else if (arg == "--include-dir" || arg == "-i") {
            std::cout << include_dir_val << std::endl;
            return 0;
        } else if (arg == "--lib-dir" || arg == "-l") {
            std::cout << lib_dir_val << std::endl;
            return 0;
        } else if (arg == "--cmake-dir" || arg == "-c") {
            std::cout << cmake_dir_val << std::endl;
            return 0;
        } else if (arg == "--libs" || arg == "-L") {
            std::cout << libs_val << std::endl;
            return 0;
        } else if (arg == "--extra-lib-dir" || arg == "-e") {
            std::cout << extra_lib_dir_val << std::endl;
            return 0;
        } else if (arg == "--extra-libs" || arg == "-E") {
            std::cout << extra_libs_val << std::endl;
            return 0;
        }
    }

    print_config_usage();
    return 1;
}

} // namespace openmc

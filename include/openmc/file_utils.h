#ifndef OPENMC_FILE_UTILS_H
#define OPENMC_FILE_UTILS_H

#include <string>

namespace openmc {

// NOTE: This is a thin wrapper over std::filesystem because we
// pass strings around a lot. Objects like settings::path_input
// are extern std::string to play with other libraries and languages

//! Determine if a path is a directory
//! \param[in] path Path to check
//! \return Whether the path is a directory
bool dir_exists(const std::string& path);

//! Determine if a file exists
//! \param[in] filename Path to file
//! \return Whether file exists
bool file_exists(const std::string& filename);

//! Determine directory containing given file
//! \param[in] filename Path to file
//! \return Name of directory containing file excluding the final directory
//! separator
std::string dir_name(const std::string& filename);

// Gets the file extension of whatever string is passed in. This is defined as
// a sequence of strictly alphanumeric characters which follow the last period,
// i.e. at least one alphabet character is present, and zero or more numbers.
// If such a sequence of characters is not found, an empty string is returned.
std::string get_file_extension(const std::string& filename);

} // namespace openmc

#endif // OPENMC_FILE_UTILS_H

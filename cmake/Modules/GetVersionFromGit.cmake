# GetVersionFromGit.cmake
# Standalone script to retrieve versioning information from Git or .git_archival.txt.
# Customizable for any project by setting variables before including this file.

# Configurable variables:
#   - VERSION_PREFIX: Prefix for version tags (default: "v").
#   - VERSION_SUFFIX: Suffix for version tags (default: "[~+-]([a-zA-Z0-9]+)").
#   - VERSION_REGEX: Regex to extract version (default: "(?[0-9]+\\.[0-9]+\\.[0-9]+)").
#   - ARCHIVAL_FILE: Path to .git_archival.txt (default: "${CMAKE_SOURCE_DIR}/.git_archival.txt").
#   - DESCRIBE_NAME_KEY: Key for describe name in .git_archival.txt (default: "describe-name: ").
#   - COMMIT_HASH_KEY: Key for commit hash in .git_archival.txt (default: "commit: ").

# Default Format Example:
#   1.2.3 v1.2.3 v1.2.3-rc1

set(VERSION_PREFIX "v" CACHE STRING "Prefix used in version tags")
set(VERSION_SUFFIX "[~+-]([a-zA-Z0-9]+)" CACHE STRING "Suffix used in version tags")
set(VERSION_REGEX "?([0-9]+\\.[0-9]+\\.[0-9]+)" CACHE STRING "Regex for extracting version")
set(ARCHIVAL_FILE "${CMAKE_SOURCE_DIR}/.git_archival.txt" CACHE STRING "Path to .git_archival.txt")
set(DESCRIBE_NAME_KEY "describe-name: " CACHE STRING "Key for describe name in .git_archival.txt")
set(COMMIT_HASH_KEY "commit: " CACHE STRING "Key for commit hash in .git_archival.txt")


# Combine prefix and regex
set(VERSION_REGEX_WITH_PREFIX "^${VERSION_PREFIX}${VERSION_REGEX}")

# Ensure Git is available
find_package(Git REQUIRED)

# Attempt to retrieve version from Git
if(EXISTS "${CMAKE_SOURCE_DIR}/.git" AND GIT_FOUND)
    message(STATUS "Using git describe for versioning")

    # Extract the version string
    execute_process(
        COMMAND git describe --tags --dirty
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE VERSION_STRING
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    # Extract the commit hash
    execute_process(
        COMMAND git rev-parse HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
    message(STATUS "Using archival file for versioning: ${ARCHIVAL_FILE}")
    if(EXISTS "${ARCHIVAL_FILE}")
        file(READ "${ARCHIVAL_FILE}" ARCHIVAL_CONTENT)

        # Extract the describe-name line
        string(REGEX MATCH "${DESCRIBE_NAME_KEY}([^\\n]+)" VERSION_STRING "${ARCHIVAL_CONTENT}")
        if(VERSION_STRING MATCHES "${DESCRIBE_NAME_KEY}(.*)")
            set(VERSION_STRING "${CMAKE_MATCH_1}")
        else()
            message(FATAL_ERROR "Could not extract version from ${ARCHIVAL_FILE}")
        endif()

        # Extract the commit hash
        string(REGEX MATCH "${COMMIT_HASH_KEY}([a-f0-9]+)" COMMIT_HASH "${ARCHIVAL_CONTENT}")
        if(COMMIT_HASH MATCHES "${COMMIT_HASH_KEY}([a-f0-9]+)")
            set(COMMIT_HASH "${CMAKE_MATCH_1}")
        else()
            message(FATAL_ERROR "Could not extract commit hash from ${ARCHIVAL_FILE}")
        endif()
    else()
        message(FATAL_ERROR "Neither git describe nor ${ARCHIVAL_FILE} is available for versioning.")
    endif()
endif()

# Ensure version string format
if(VERSION_STRING MATCHES "${VERSION_REGEX_WITH_PREFIX}")
    set(VERSION_NO_SUFFIX "${CMAKE_MATCH_1}")
else()
    message(FATAL_ERROR "Invalid version format: Missing base version in ${VERSION_STRING}")
endif()

# Check for development state
if(VERSION_STRING MATCHES "-([0-9]+)-g([0-9a-f]+)")
    set(DEV_STATE "true")
    set(COMMIT_COUNT "${CMAKE_MATCH_1}")
    string(REGEX REPLACE "-([0-9]+)-g([0-9a-f]+)" "" VERSION_WITHOUT_META "${VERSION_STRING}")
else()
    set(DEV_STATE "false")
    set(VERSION_WITHOUT_META "${VERSION_STRING}")
endif()

# Split and set version components
string(REPLACE "." ";" VERSION_LIST "${VERSION_NO_SUFFIX}")
list(GET VERSION_LIST 0 VERSION_MAJOR)
list(GET VERSION_LIST 1 VERSION_MINOR)
list(GET VERSION_LIST 2 VERSION_PATCH)

# Increment patch number for dev versions
if(DEV_STATE)
    math(EXPR VERSION_PATCH "${VERSION_PATCH} + 1")
endif()

# Export variables
set(OPENMC_VERSION_MAJOR "${VERSION_MAJOR}")
set(OPENMC_VERSION_MINOR "${VERSION_MINOR}")
set(OPENMC_VERSION_PATCH "${VERSION_PATCH}")
set(OPENMC_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
set(OPENMC_COMMIT_HASH "${COMMIT_HASH}")
set(OPENMC_DEV_STATE "${DEV_STATE}")
set(OPENMC_COMMIT_COUNT "${COMMIT_COUNT}")

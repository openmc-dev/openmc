# Function to generate and install a Python script
# This function creates a Python script with a specified name, writes content to it, 
# and installs it to a specified location.
#
# Arguments:
#   SCRIPT_NAME - The name of the script to be generated (e.g., 'my_script').
#
function(generate_and_install_python_script SCRIPT_NAME)
    if(SKBUILD)
        # Strip any leading/trailing whitespace from the script name
        string(STRIP "${SCRIPT_NAME}" CLEAN_SCRIPT_NAME)

        # Define the output path for the generated script using the cleaned name
        set(GENERATED_SCRIPT "${CMAKE_BINARY_DIR}/${CLEAN_SCRIPT_NAME}")

        # Generate the script content directly
        file(WRITE "${GENERATED_SCRIPT}" "#!/usr/bin/env python3\n\n")
        file(APPEND "${GENERATED_SCRIPT}" "import os\n")
        file(APPEND "${GENERATED_SCRIPT}" "import sys\n")
        file(APPEND "${GENERATED_SCRIPT}" "import sysconfig\n\n")
        file(APPEND "${GENERATED_SCRIPT}" "if __name__ == '__main__':\n")
        file(APPEND "${GENERATED_SCRIPT}" "    os.execv(\n")
        file(APPEND "${GENERATED_SCRIPT}" "        os.path.join(sysconfig.get_path('platlib'), '${SKBUILD_PROJECT_NAME}', '${CMAKE_INSTALL_BINDIR}', '${CLEAN_SCRIPT_NAME}'),\n")
        file(APPEND "${GENERATED_SCRIPT}" "        sys.argv,\n")
        file(APPEND "${GENERATED_SCRIPT}" "    )\n")

        # Install the generated script
        install(
            PROGRAMS "${GENERATED_SCRIPT}"
            DESTINATION "${SKBUILD_SCRIPTS_DIR}"
        )
    endif()
endfunction()

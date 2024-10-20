# Function to create and install the binary scripts
function(generate_and_install_scripts SCRIPT_LIST)
    foreach(SCRIPT_NAME ${SCRIPT_LIST})
        # Define the output path for each generated script
        set(GENERATED_SCRIPT "${CMAKE_BINARY_DIR}/${SCRIPT_NAME}")

        # Generate the script content directly
        file(WRITE "${GENERATED_SCRIPT}" "#!/usr/bin/env python3\n\n")
        file(APPEND "${GENERATED_SCRIPT}" "import os\n")
        file(APPEND "${GENERATED_SCRIPT}" "import sys\n")
        file(APPEND "${GENERATED_SCRIPT}" "import sysconfig\n\n")
        file(APPEND "${GENERATED_SCRIPT}" "if __name__ == '__main__':\n")
        file(APPEND "${GENERATED_SCRIPT}" "    os.execv(\n")
        file(APPEND "${GENERATED_SCRIPT}" "        os.path.join(sysconfig.get_path('platlib'), '${SKBUILD_PROJECT_NAME}', '${CMAKE_INSTALL_BINDIR}', '${SCRIPT_NAME}'),\n")
        file(APPEND "${GENERATED_SCRIPT}" "        sys.argv,\n")
        file(APPEND "${GENERATED_SCRIPT}" "    )\n")

        # Install the generated script
        install(
            PROGRAMS "${GENERATED_SCRIPT}"
            DESTINATION "${SKBUILD_SCRIPTS_DIR}"
        )
    endforeach()
endfunction()

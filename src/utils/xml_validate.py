#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import glob
import lxml.etree as etree
from subprocess import call
from optparse import OptionParser

# Command line parsing
parser = OptionParser()
parser.add_option('-r', '--relaxng-path', dest='relaxng',
                  help="Path to RelaxNG files.")
parser.add_option('-i', '--input-path', dest='inputs', default=os.getcwd(),
                  help="Path to OpenMC input files." )
(options, args) = parser.parse_args()

# Colored output
if sys.stdout.isatty():
    OK = '\033[92m'
    FAIL = '\033[91m'
    NOT_FOUND = '\033[93m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
else:
    OK = ''
    FAIL = ''
    ENDC = ''
    BOLD = ''
    NOT_FOUND = ''

# Get absolute paths
if options.relaxng is not None:
    relaxng_path = os.path.abspath(options.relaxng)
if options.inputs is not None:
    inputs_path = os.path.abspath(options.inputs)

# Search for relaxng path if not set
if options.relaxng is None:
    xml_validate_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    if "bin" in xml_validate_path:
        relaxng_path = os.path.join(xml_validate_path, "..", "share", "relaxng")
    elif os.path.join("src", "utils") in xml_validate_path:
        relaxng_path = os.path.join(xml_validate_path, "..", "relaxng")
    else:
        raise Exception("Set RelaxNG path with -r command line option.")
if not os.path.exists(relaxng_path):
    raise Exception("RelaxNG path: {0} does not exist, set with -r "
                    "command line option.".format(relaxng_path))

# Make sure there are .rng files in RelaxNG path
rng_files = glob.glob(os.path.join(relaxng_path, "*.rng"))
if len(rng_files) == 0:
    raise Exception("No .rng files found in RelaxNG "
                    "path: {0}.".format(relaxng_path))

# Get list of xml input files
xml_files = glob.glob(os.path.join(inputs_path, "*.xml"))
if len(xml_files) == 0:
    raise Exception("No .xml files found at input path: {0}"
                    ".".format(inputs_path))

# Begin loop around input files
for xml_file in xml_files:

    text = "Validating {0}".format(os.path.basename(xml_file))
    print(text + '.'*(30 - len(text)), end="")

    # Validate the XML file
    try:
        xml_tree = etree.parse(xml_file)
    except etree.XMLSyntaxError as e:
        print(BOLD + FAIL + '[XML ERROR]' + ENDC)
        print("    {0}".format(e))
        continue

    # Get xml_filename prefix
    xml_prefix = os.path.basename(xml_file)
    xml_prefix = xml_prefix.split(".")[0]

    # Search for rng file
    rng_file = os.path.join(relaxng_path, xml_prefix + ".rng")
    if rng_file in rng_files:

        # read in RelaxNG
        relaxng_doc = etree.parse(rng_file)
        relaxng = etree.RelaxNG(relaxng_doc)

        # validate xml file again RelaxNG
        try:
            relaxng.assertValid(xml_tree)
            print(BOLD + OK + '[VALID]' + ENDC)
        except (etree.DocumentInvalid, TypeError) as e:
            print(BOLD + FAIL + '[NOT VALID]' + ENDC)
            print("    {0}".format(e))

    # RNG file does not exist
    else:
        print(BOLD + NOT_FOUND + '[NO RELAXNG FOUND]' + ENDC)

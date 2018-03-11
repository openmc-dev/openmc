#!/usr/bin/env python3

import glob
from string import whitespace
from sys import exit
from textwrap import TextWrapper


################################################################################
# Fortran reading/parsing backend.
################################################################################


class LineOfCode(object):
    """Contains and provides basic info about a line of Fortran 90 code."""
    def __init__(self, number, content, is_continuation=False,
                 string_terminator=None):
        # Initialize member variables.
        self.number = number
        self.content = content
        self.is_continuation = is_continuation
        self.is_continued = False
        assert (string_terminator is None or
                string_terminator == "'" or
                string_terminator == '"')
        self.initial_string_terminator = string_terminator
        self.final_string_terminator = None
        self.is_comment_only = False

        # Parse the string.
        self.parse()

    def __repr__(self):
        rep = 'LineOfCode: line number = {0:d}\n'.format(self.number)
        rep += '\tis_continuation = ' + str(self.is_continuation) + '\n'
        rep += '\tis_continued = ' + str(self.is_continued) + '\n'
        rep += ('\tinitial_string_terminator = '
                + str(self.initial_string_terminator) + '\n')
        rep += ('\tfinal_string_terminator = '
                + str(self.final_string_terminator) + '\n')
        rep += '\tis_comment_only = ' + str(self.is_comment_only) + '\n'
        rep += '\tContent:\n' + self.content
        return rep

    def __str__(self):
        return repr(self)

    def parse(self):
        # Initialize the string variables.
        if self.initial_string_terminator == "'":
            in_a_single_quote = True
            in_a_double_quote = False
        elif self.initial_string_terminator == '"':
            in_a_single_quote = False
            in_a_double_quote = True
        else:
            in_a_single_quote = False
            in_a_double_quote = False

        # Initialize other variables.
        in_a_comment = False
        ends_with_ampersand = False
        has_real_content = False

        # Parse through the line.
        for char in self.content:
            # Check for the start or end of strings.
            if char == "'" and in_a_single_quote:
                in_a_single_quote = False
            elif char == "'" and not in_a_double_quote:
                in_a_single_quote = True
            elif char == '"' and in_a_double_quote:
                in_a_double_quote = False
            elif char == '"' and not in_a_single_quote:
                in_a_double_quote = True

            # Check for the start of a comment.
            if (char == "!" and not in_a_single_quote
                 and not in_a_double_quote):
                in_a_comment = True
                break # Don't care about anything in a comment.

            # Check for a continuation ampersand.
            if char == "&":
                ends_with_ampersand = True
            elif all([char != w for w in whitespace]):
                ends_with_ampersand = False

            # Check to see if there is any real content (non-whitespace, not in
            # a comment)
            if (not in_a_comment) and all([char != w for w in whitespace]):
                has_real_content = True

            # Make sure nothing went terribly wrong.
            assert not (in_a_single_quote and in_a_double_quote)
            assert not (in_a_single_quote and in_a_comment)
            assert not (in_a_double_quote and in_a_comment)

        # Is this line continued onto the next?
        if ends_with_ampersand:
            self.is_continued = True

        # Is a multiline string continued on the next line?
        if in_a_single_quote:
            assert self.is_continued
            self.final_string_terminator = "'"
        if in_a_double_quote:
            assert self.is_continued
            self.final_string_terminator = '"'

        # Is this line a comment-only line?
        if (not has_real_content) and in_a_comment:
            self.is_comment_only = True

    def get_indent(self):
        """Return the number of indentation spaces prefixing this line."""
        return len(self.content) - len(self.content.lstrip(' '))

    def contains_whitespace_only(self):
        """Return True/False if all characters in the line are whitespace."""
        if len(self.content.strip('\n')) == 0: return False
        is_ws = [ any([char == ws for ws in whitespace])
                  for char in self.content.strip('\n') ]
        return all(is_ws)

    def has_trailing_whitespace(self):
        """Return True/False if this line ends with a whitespace character."""
        stripped = self.content.strip('\n')
        if len(stripped) == 0: return False
        return any([stripped[-1] == ws for ws in whitespace])

    def contains_tab(self):
        """Return True/False if a tab character appears in the line."""
        return '\t' in self.content


def read_lines_of_code(fname):
    line_num = 0
    cont = False
    str_term = None
    with open(fname) as fin:
        for line in fin:
            loc = LineOfCode(line_num, line, is_continuation=cont,
                             string_terminator=str_term)
            cont = loc.is_continued
            str_term = loc.final_string_terminator
            line_num += 1
            yield loc


################################################################################
# Error checking.
################################################################################


def print_error(fname, line_number, err_msg):
    header = "Error in file {0}, line {1}:".format(fname, line_number + 1)

    tw = TextWrapper(width=80, initial_indent=' '*4, subsequent_indent=' '*4)
    body = '\n'.join(tw.wrap(err_msg))

    print(header + '\n' + body + '\n')


def check_source(fname):
    """Make sure the given Fortran source file meets OpenMC standards.

    If errors are found, messages will be printed to the screen describing the
    error.  The function will return True if no errors were found or False
    otherwise.
    """
    good_code = True
    base_indent = 0

    for loc in read_lines_of_code(fname):
        # Check for extra whitespace errors.
        if loc.contains_whitespace_only():
            good_code = False
            print_error(fname, loc.number, "Line contains whitespace "
                 "characters but no content.  Please remove whitespace.")
        elif loc.has_trailing_whitespace():
            good_code = False
            print_error(fname, loc.number, "Line has trailing whitespace after"
                 " the content.  Please remove trailing whitespace.")
        if loc.contains_tab():
            good_code = False
            print_error(fname, loc.number, "Line contains a tab character.  "
                 "Please replace with single whitespace characters.")

        # Check indentation.
        current_indent = loc.get_indent()
        if ((not loc.is_continuation) and (not loc.is_comment_only) and
              (not loc.contains_whitespace_only()) and current_indent % 2 != 0):
            good_code = False
            print_error(fname, loc.number, "Line is indented an odd number of "
                 "spaces.  All non-continuation lines should be indented an "
                 "even number of spaces.")
        if loc.is_continuation and current_indent < base_indent + 5:
            good_code = False
            print_error(fname, loc.number, "Continuation lines must be "
                 "indented by at least 5 spaces, but this line is indented {0}"
                 " spaces.".format(current_indent - base_indent))

        # Set base indentation for next lines.
        if not loc.is_continuation:
            base_indent = current_indent

    return good_code


################################################################################
# Main loop.
################################################################################


if __name__ == '__main__':
    # Get a list of the F90 source files.
    f90_files = glob.glob('../src/*.F90')

    # Make sure we found the source files.
    assert len(f90_files) != 0, 'No .F90 source files found.'

    # Make sure all the F90 source files meet our standards.
    good_code = [check_source(fname) for fname in f90_files]
    if not all(good_code):
        print("ERROR: The Fortran source code does not meet OpenMC's standards")
        exit(-1)
    else:
        print("SUCCESS! Looks like the Fortran source meets our standards")
        exit(0)

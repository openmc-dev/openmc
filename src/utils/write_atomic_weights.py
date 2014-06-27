#!/usr/bin/env python2
"""Write expand_natural_element subroutine using NIST atomic weight data.

In order to use this program, follow these steps:
*  Make a data file from the NIST Atomic Weights and Isotopic Compositions
   database <http://nist.gov/pml/data/comp.cfm>.  Check the 'Linearized ASCII
   Output' and 'All isotopes' options.  Please add access date and reference
   information to the top of file.  The data file should look like the
   NIST_weights.txt file in the openmc /src/utils directory.
*  Run this program with the sytax 'write_atomic_weights.py DATA TEMPLATE OUT'.
   TEMPLATE should be an input_xml.F90 file with the 'Start of the natural
   element cases' and 'End of the natural element cases' flags.  The OUT file
   will be written exactly like the TEMPLATE file but with the inclusion of the
   atomic weight data. Note that TEMPLATE and OUT can be the same file; in that
   case, the template will be overwritten.  A typical invokation would be:
   ./write_atomic_weights.py NIST_weights.txt ../input_xml.F90 ../input_xml.F90

Usage information can be obtained by running 'write_atomic_weights.py --help':

    usage: write_atomic_weights.py [-h] DATA TEMPLATE OUT

    Rewrite an input_xml.F90 file with atomic weight data included in the
    expand_natural_element subroutine.

    positional arguments:
      DATA        Text file containing NIST atomic weight data.
      TEMPLATE    input_xml.F90 file with 'Start of', 'End of...' flags.
      OUT         Name of .F90 file that will be written.

    optional arguments:
      -h, --help  show this help message and exit

"""

from __future__ import print_function

import argparse


class Element():
    """Structure for containing atomic data."""
    def __init__(self, symbol, z, stand_atom_weight):
        self.symbol = symbol
        self.z = z
        self.sam = stand_atom_weight
        self.a_list = []
        self.ram_list = []
        self.frac_list = []
    def add_isotope(self, a, rel_atom_mass, iso_frac):
        self.a_list.append(a)
        self.ram_list.append(rel_atom_mass)
        self.frac_list.append(iso_frac)


def _parse_args():
    # Create argument parser.
    parser = argparse.ArgumentParser(description=''.join(['Rewrite an ',
        'input_xml.F90 file with atomic weight data included in the ',
        'expand_natural_element subroutine.']))
    parser.add_argument('data', metavar='DATA', type=str,
                        help='Text file containing NIST atomic weight data.')
    parser.add_argument('template', metavar='TEMPLATE', type=str,
                  help="input_xml.F90 file with 'Start of', 'End of...' flags.")
    parser.add_argument('output', metavar='OUT', type=str,
                        help='Name of .F90 file that will be written.')

    # Parse and return commandline arguments.
    return parser.parse_args()


def read_NIST_file(fname):
    """Return dataset from NIST linearized ASCII atomic weights database."""
    fin = open(fname)

    line = ' '
    elements = []
    while line != '':
        # Read lines until the start of a data entry is reached.
        line = fin.readline()
        if line[:13] != 'Atomic Number': continue

        # Extract the atomic number from this line.
        z = int(line[16:].strip())

        # Read the next line and extract the elemental symbol.
        line = fin.readline()
        assert line[:13] == 'Atomic Symbol'
        symbol = line[16:].strip()

        # Read the next line and extract the mass number.
        line = fin.readline()
        assert line[:11] == 'Mass Number'
        a = int(line[14:].strip())

        # Read the next line and extract the relative atomic mass if given.
        line = fin.readline()
        assert line[:20] == 'Relative Atomic Mass'
        if len(line) < 25:
            # No relative atomic mass given; isotope is synthetic or trace.
            # Move on to the next data entry.
            continue
        ram = float(line[23:].split('(')[0]) # Number between index 23 and '('

        # Read the next line and extract the isotopic fraction if given.
        line = fin.readline()
        assert line[:20] == 'Isotopic Composition'
        if len(line) < 25:
            # No fraction given; isotope is synthetic or trace.
            # Move on to the next data entry.
            continue
        frac = float(line[23:].split('(')[0]) # Number between index 23 and '('

        # Read the next line and extract the standard atomic weight.
        line = fin.readline()
        assert line[:22] == 'Standard Atomic Weight'
        if len(line) < 27:
            # No standard weight given; isotope is synthetic or trace.
            # Move on to the next data entry.
            continue
        sam = float(line[25:].split('(')[0]) # Number between index 25 and '('

        # All relevant data has now been read; add it to the elements array.
        # Check to see if this element is already in the list.
        if len(elements) > 0:
            if elements[-1].z == z:
                # Element is already present; add this isotope.
                elements[-1].add_isotope(a, ram, frac)
                # Move on to the next data entry.
                continue
        # Element is not already present; add it.
        elements.append(Element(symbol, z, sam))
        elements[-1].add_isotope(a, ram, frac)

    fin.close()
    return elements


def print_elements(elements):
    """Print atomic data from a list of elements."""
    for el in elements:
        print('sym={0}, z={1}, sam={2} :'.format(el.symbol, el.z, el.sam))
        for i in range(len(el.a_list)):
            print('    a={0}, ram={1}, frac={2}'.format(el.a_list[i],
                                               el.ram_list[i], el.frac_list[i]))

def write_code(elements, input_fname, output_fname):
    """Write code using input template and atomic data to a file."""
    # Read the template.  The code will be contained in a list where every line
    # is a seperate list element.
    fin = open(input_fname, mode='r')
    code = fin.readlines()
    fin.close()

    # Find and remove the signal phrase.
    start_ind = code.index('    ! Start of the natural element cases.\n')
    end_ind = code.index('    ! End of the natural element cases.\n')
    del code[start_ind+1:end_ind]
    ind = start_ind + 1

    # Insert atomic data inside case statements.
    for elem in elements:
        # Insert case statements ex: "case('h')", "case('he')", "case('li')"
        code.insert(ind, "    case ('{0}')\n".format(elem.symbol.lower()))
        ind += 1

        # Handle elements with incomplete cross section evaluations.
        if elem.symbol == 'C':
            st = '      ! No evaluations split up Carbon into isotopes yet\n'
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, line)
                ind += 1

        elif elem.symbol == 'O':
            st = '      if (default_expand == JEFF_32) then\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = ''.join(['      elseif (default_expand >= JENDL_32 .and.',
                          ' default_expand <= JENDL_40) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(16, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            ram = ((elem.ram_list[0]*elem.frac_list[0]
                    + elem.ram_list[2]*elem.frac_list[2])
                   /(elem.frac_list[0] + elem.frac_list[2]))
            el.add_isotope(16, ram, elem.frac_list[0] + elem.frac_list[2])
            el.add_isotope(17, elem.ram_list[1], elem.frac_list[1])
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'V':
            st = ''.join(['      if (default_expand == ENDF_BVII0 .or. ',
                          'default_expand == JEFF_311 &\n'])
            code.insert(ind, st)
            ind += 1
            st = '           .or. default_expand == JEFF_32 .or. &\n'
            code.insert(ind, st)
            ind += 1
            st = ''.join(['           (default_expand >= JENDL_32 .and. ',
                          'default_expand <= JENDL_33)) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Zn':
            st = ''.join(['      if (default_expand == ENDF_BVII0 .or. ',
                          'default_expand == &\n'])
            code.insert(ind, st)
            ind += 1
            st = '           JEFF_311 .or. default_expand == JEFF_312) then\n'
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Ga':
            st = ''.join(['      if (default_expand == JEFF_311 .or. ',
                          'default_expand == JEFF_312) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Ta':
            st = '      if (default_expand == ENDF_BVII0 .or. &\n'
            code.insert(ind, st)
            ind += 1
            st = ''.join(['           (default_expand >= JEFF_311 .and. ',
                          'default_expand <= JEFF_312) .or. &\n'])
            code.insert(ind, st)
            ind += 1
            st = ''.join(['           (default_expand >= JENDL_32 .and. ',
                          'default_expand <= JENDL_40)) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(181, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'W':
            st = ''.join(['      if (default_expand == ENDF_BVII0 .or. ',
                          'default_expand == JEFF_311 &\n'])
            code.insert(ind, st)
            ind += 1
            st = '           .or. default_expand == JEFF_312 .or. &\n'
            code.insert(ind, st)
            ind += 1
            st = ''.join(['           (default_expand >= JENDL_32 .and. ',
                          'default_expand <= JENDL_33)) then\n'])
            code.insert(ind, st)
            ind += 1
            st = '        ! Combine W-180 with W-182\n'
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            ram = ((elem.ram_list[0]*elem.frac_list[0]
                    + elem.ram_list[1]*elem.frac_list[1])
                   /(elem.frac_list[0] + elem.frac_list[1]))
            el.add_isotope(182, ram, elem.frac_list[0] + elem.frac_list[1])
            for i in range(len(elem.a_list))[2:]:
                el.add_isotope(elem.a_list[i], elem.ram_list[i],
                               elem.frac_list[i])
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Os':
            st = ''.join(['      if (default_expand == JEFF_311 .or. ',
                          'default_expand == JEFF_312) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Pt':
            st = ''.join(['      if (default_expand == JEFF_311 .or. ',
                          'default_expand == JEFF_312) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        elif elem.symbol == 'Tl':
            st = ''.join(['      if (default_expand == JEFF_311 .or. ',
                          'default_expand == JEFF_312) then\n'])
            code.insert(ind, st)
            ind += 1
            el = Element(elem.symbol, elem.z, elem.sam)
            el.add_isotope(0, elem.sam, 1.0)
            for line in code_lines(el):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      else\n'
            code.insert(ind, st)
            ind += 1
            for line in code_lines(elem):
                code.insert(ind, ''.join(['  ', line]))
                ind += 1

            st = '      end if\n'
            code.insert(ind, st)
            ind += 1

        else:
            for line in code_lines(elem):
                code.insert(ind, line)
                ind += 1

        code.insert(ind, '\n')
        ind += 1

    # Write the output file.
    fout = open(output_fname, mode='w')
    fout.write(''.join(code))
    fout.close()


def code_lines(el):
    """Return a list of lines of code that divide element el into isotopes."""
    line1 = "      call list_names % append('{0}{1:03d}.' // xs)\n"
    line2 = "      if (density > 0.0) then\n"
    line3 = "        call list_density % append(density * {0:.14f}_8)\n"
    line4 = "      else\n"
    line5 = "        call list_density % append(density * {0:.14f}_8)\n"
    line6 = "      end if\n"
    lines = []
    for i in range(len(el.a_list)):
        lines.append(line1.format(el.z, el.a_list[i]))
    lines.append(line2)
    for i in range(len(el.a_list)):
        lines.append(line3.format(el.frac_list[i]))
    lines.append(line4)
    for i in range(len(el.a_list)):
        lines.append(line5.format(el.ram_list[i] / el.sam * el.frac_list[i]))
    lines.append(line6)
    return lines


if __name__ == '__main__':
    args = _parse_args()
    elements = read_NIST_file(args.data)
    write_code(elements, args.template, args.output)

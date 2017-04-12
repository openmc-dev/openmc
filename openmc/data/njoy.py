#!/usr/bin/env python

import argparse
from collections import namedtuple
from io import StringIO
import os
import shutil
from subprocess import Popen, PIPE, STDOUT
import sys
import tempfile

from .endf import Evaluation

_NJOY_2012 = False

# For a given MAT number, give a name for the ACE table and a list of ZAID
# identifiers
ThermalTuple = namedtuple('ThermalTuple', ['name', 'zaids', 'nmix'])
_THERMAL_DATA = {
    1: ThermalTuple('hh2o', [1001], 1),
    2: ThermalTuple('parah', [1001], 1),
    3: ThermalTuple('orthoh', [1001], 1),
    5: ThermalTuple('hyh2', [1001], 1),
    7: ThermalTuple('hzrh', [1001], 1),
    8: ThermalTuple('hcah2', [1001], 1),
    10: ThermalTuple('hice', [1001], 1),
    11: ThermalTuple('dd2o', [1002], 1),
    12: ThermalTuple('parad', [1002], 1),
    13: ThermalTuple('orthod', [1002], 1),
    26: ThermalTuple('be', [4009], 1),
    27: ThermalTuple('bebeo', [4009], 1),
    31: ThermalTuple('graph', [6000, 6012, 6013], 1),
    33: ThermalTuple('lch4', [1001], 1),
    34: ThermalTuple('sch4', [1001], 1),
    37: ThermalTuple('hch2', [1001], 1),
    39: ThermalTuple('lucite', [1001], 1),
    40: ThermalTuple('benz', [1001, 6000, 6012], 2),
    41: ThermalTuple('od2o', [8016, 8017, 8018], 1),
    43: ThermalTuple('sisic', [14028, 14029, 14030], 1),
    44: ThermalTuple('csic', [6000, 6012, 6013], 1),
    46: ThermalTuple('obeo', [8016, 8017, 8018], 1),
    47: ThermalTuple('sio2-a', [8016, 8017, 8018, 14028, 14029, 14030], 3),
    48: ThermalTuple('uuo2', [92238], 1),
    49: ThermalTuple('sio2-b', [8016, 8017, 8018, 14028, 14029, 14030], 3),
    50: ThermalTuple('oice', [8016, 8017, 8018], 1),
    52: ThermalTuple('mg24', [12024], 1),
    53: ThermalTuple('al27', [13027], 1),
    55: ThermalTuple('yyh2', [39089], 1),
    56: ThermalTuple('fe56', [26056], 1),
    58: ThermalTuple('zrzrh', [40000, 40090, 40091, 40092, 40094, 40096], 1),
    59: ThermalTuple('cacah2', [20040, 20042, 20043, 20044, 20046, 20048], 1),
    75: ThermalTuple('ouo2', [8016, 8017, 8018], 1),
}


_PENDF_TEMPLATE = """
reconr / %%%%%%%%%%%%%%%%%%% Reconstruct XS for neutrons %%%%%%%%%%%%%%%%%%%%%%%
20 22
'{library} PENDF for {zsymam}'/
{mat} 2/
0.001 0.0 0.003/ err tempr errmax
'{library}: {zsymam}'/
'Processed by NJOY'/
0/
stop
"""

_ACE_TEMPLATE = """
reconr / %%%%%%%%%%%%%%%%%%% Reconstruct XS for neutrons %%%%%%%%%%%%%%%%%%%%%%%
20 21
'{library} PENDF for {zsymam}'/
{mat} 2/
0.001 0.0 0.003/ err tempr errmax
'{library}: {zsymam}'/
'Processed by NJOY'/
0/
broadr / %%%%%%%%%%%%%%%%%%%%%%% Doppler broaden XS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
20 21 22
{mat} {num_temp} 0 0 0. /
0.001 -1.0e6 0.003 /
{temps}
0/
heatr / %%%%%%%%%%%%%%%%%%%%%%%%% Add heating kerma %%%%%%%%%%%%%%%%%%%%%%%%%%%%
20 22 23 /
{mat} {num_partials} /
{heatr_partials} /
purr / %%%%%%%%%%%%%%%%%%%%%%%% Add probability tables %%%%%%%%%%%%%%%%%%%%%%%%%
20 23 24
{mat} {num_temp} 5 20 64 /
{temps}
1.e10 1.e4 1.e3 1.e2 1.e1
0/
"""

_ACE_TEMPLATE_ACER = """acer / %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make ACE file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
20 24 0 {nace} {ndir}
1 0 1 .{ext} /
'{library}: {zsymam} at {temperature}'/
{mat} {temperature}
1 1/
/
"""

_ACE_THERMAL_TEMPLATE = """
reconr / %%%%%%%%%%%%%%%%%%% Reconstruct XS for neutrons %%%%%%%%%%%%%%%%%%%%%%%
20 22
'{library} PENDF for {zsymam}'/
{mat} 2/
0.001 0. 0.001/ err tempr errmax
'{library}: PENDF for {zsymam}'/
'Processed by NJOY'/
0/
broadr / %%%%%%%%%%%%%%%%%%%%%%% Doppler broaden XS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
20 22 23
{mat} 1 0 0 0./
0.001 2.0e+6 0.001/ errthn thnmax errmax
{temperature}
0/
thermr / %%%%%%%%%%%%%%%% Add thermal scattering data (free gas) %%%%%%%%%%%%%%%
0 23 62
0 {mat} 12 1 1 0 {iform} 1 221 1/
{temperature}
0.001 4.0
thermr / %%%%%%%%%%%%%%%% Add thermal scattering data (bound) %%%%%%%%%%%%%%%%%%
60 62 27
{mat_thermal} {mat} 16 1 {inelastic} {elastic} {iform} {natom} 222 1/
{temperature}
0.001 4.0
acer / %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Make ACE file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
20 27 0 28 29
2 0 1 .01/
'{library}: {zsymam_thermal} processed by NJOY'/
{mat} {temperature} '{data.name}' /
{zaids} /
222 64 {mt_elastic} {elastic_type} {data.nmix} {energy_max} 2/
stop
"""

def run(commands, tapein, tapeout, stdout=False, njoy_exec='njoy'):
    # Create temporary directory -- it would be preferable to use
    # TemporaryDirectory(), but it is only available in Python 3.2
    tmpdir = tempfile.mkdtemp()
    try:
        # Copy evaluations to appropriates 'tapes'
        for tape_num, filename in tapein.items():
            tmpfilename = os.path.join(tmpdir, 'tape{}'.format(tape_num))
            shutil.copy(filename, tmpfilename)

        # Start up NJOY process
        njoy = Popen([njoy_exec], cwd=tmpdir, stdin=PIPE, stdout=PIPE,
                   stderr=STDOUT, universal_newlines=True)

        njoy.stdin.write(commands)
        njoy.stdin.flush()
        while True:
            # If process is finished, break loop
            line = njoy.stdout.readline()
            if not line and njoy.poll() is not None:
                break

            if stdout:
                # If user requested output, print to screen
                print(line, end='')

        # Copy output files back to original directory
        for tape_num, filename in tapeout.items():
            tmpfilename = os.path.join(tmpdir, 'tape{}'.format(tape_num))
            if os.path.isfile(tmpfilename):
                shutil.move(tmpfilename, filename)
    finally:
        shutil.rmtree(tmpdir)

    return njoy.returncode


def make_pendf(filename, pendf='pendf', stdout=False):
    """Generate ACE file from an ENDF file

    Parameters
    ----------
    filename : str
        Path to ENDF file
    pendf : str, optional
        Path of pointwise ENDF file to write
    stdout : bool
        Whether to display NJOY standard output

    Returns
    -------
    int
        Return code of NJOY process

    """
    ev = Evaluation(filename)
    mat = ev.material
    zsymam = ev.target['zsymam']

    # Determine name of library
    library = '{}-{}.{}'.format(*ev.info['library'])

    commands = _PENDF_TEMPLATE.format(**locals())
    tapein = {20: filename}
    tapeout = {22: pendf}
    return run(commands, tapein, tapeout, stdout)


def make_ace(filename, temperatures=None, ace='ace', xsdir='xsdir',
             keep_pendf=False, **kwargs):
    """Generate ACE file from an ENDF file

    Parameters
    ----------
    filename : str
        Path to ENDF file
    temperatures : iterable of float, optional
        Temperatures in Kelvin to produce ACE files at
    ace : str, optional
        Path of ACE file to write
    xsdir : str, optional
        Path of xsdir file to write
    keep_pendf : bool
        Whether to retain PENDF file (written to 'pendf')
    **kwargs
        Keyword arguments passed to :func:`openmc.data.njoy.run`

    Returns
    -------
    int
        Return code of NJOY process

    """
    ev = Evaluation(filename)
    mat = ev.material
    zsymam = ev.target['zsymam']

    # Determine name of library
    library = '{}-{}.{}'.format(*ev.info['library'])

    # Determine which partial heating values are needed to self-shield heating
    # for ptables
    partials = [302, 402]
    if ev.target['fissionable']:
        partials.append(318)
    heatr_partials = ' '.join(map(str, partials))
    num_partials = len(partials)

    if temperatures is None:
        temperatures = [293.15]
    num_temp = len(temperatures)
    temps = ' '.join(str(i) for i in temperatures)

    commands = _ACE_TEMPLATE.format(**locals())
    tapein = {20: filename}
    tapeout = {}
    if keep_pendf:
        tapeout[21] = 'pendf'
    fname = '{}_{:.1f}'
    for i, temperature in enumerate(temperatures):
        # Extend input with an ACER run for each temperature
        nace = 25 + 2*i
        ndir = 25 + 2*i + 1
        ext = '{:02}'.format(i + 1)
        commands += _ACE_TEMPLATE_ACER.format(**locals())

        # Indicate tapes to save for each ACER run
        tapeout[nace] = fname.format(ace, temperature)
        tapeout[ndir] = fname.format(xsdir, temperature)
    commands += 'stop\n'
    retcode = run(commands, tapein, tapeout, **kwargs)

    if retcode == 0:
        with open(ace, 'w') as ace_file, open(xsdir, 'w') as xsdir_file:
            for temperature in temperatures:
                # Get contents of ACE file
                text = open(fname.format(ace, temperature), 'r').read()

                # If the target is metastable, make sure that ZAID in the ACE file reflects
                # this by adding 400
                if ev.target['isomeric_state'] > 0:
                    mass_first_digit = int(text[3])
                    if mass_first_digit <= 2:
                        text = text[:3] + str(mass_first_digit + 4) + text[4:]

                # Concatenate into destination ACE file
                ace_file.write(text)

                # Concatenate into destination xsdir file
                text = open(fname.format(xsdir, temperature), 'r').read()
                xsdir_file.write(text)

        # Remove ACE/xsdir files for each temperature
        for temperature in temperatures:
            os.remove(fname.format(ace, temperature))
            os.remove(fname.format(xsdir, temperature))

    return retcode


def make_ace_thermal(filename, filename_thermal, temperature, output=None):
    ev = Evaluation(filename)
    mat = ev.material
    zsymam = ev.target['zsymam']

    ev_thermal = Evaluation(filename_thermal)
    mat_thermal = ev_thermal.material
    zsymam_thermal = ev_thermal.target['zsymam']

    data = _THERMAL_DATA[mat_thermal]
    zaids = ' '.join(str(zaid) for zaid in data.zaids[:3])
    energy_max = ev_thermal.info['energy_max']

    # Determine name of library
    library = '{}-{}.{}'.format(*ev_thermal.info['library'])

    # Determine if thermal elastic is present
    if (7, 2) in ev_thermal.section:
        elastic = 1
        mt_elastic = 223

        # Determine whether elastic is incoherent (0) or coherent (1)
        file_obj = StringIO(ev_thermal.section[7, 2])
        elastic_type = openmc.data.endf.get_head_record(file_obj)[2]
    else:
        elastic = 0
        mt_elastic = 0
        elastic_type = 0

    # Determine number of principal atoms
    file_obj = StringIO(ev_thermal.section[7, 4])
    items = openmc.data.endf.get_head_record(file_obj)
    items, values = openmc.data.endf.get_list_record(file_obj)
    natom = int(values[5])

    # A few parameters are different in NJOY 99 and NJOY 2012
    if _NJOY_2012:
        iform = 0
        inelastic = 2
    else:
        iform = ''
        inelastic = 4

    with tempfile.TemporaryDirectory() as tmpdirname:
        # Copy evaluation to tape20
        shutil.copy(filename, os.path.join(tmpdirname, 'tape20'))
        shutil.copy(filename_thermal, os.path.join(tmpdirname, 'tape60'))

        commands = _ACE_THERMAL_TEMPLATE.format(**locals())
        njoy = run(['njoy'], cwd=tmpdirname, input=commands, stdout=PIPE,
                   stderr=STDOUT, universal_newlines=True)
        if njoy.returncode != 0:
            print('Failed to produce ACE for {}'.format(filename_thermal))

        if output is None:
            output = data.name + '.ace'
        output_xsdir = data.name + '.xsdir'
        output_njoy = data.name + '.output'

        ace_file = os.path.join(tmpdirname, 'tape28')
        xsdir_file = os.path.join(tmpdirname, 'tape29')
        output_file = os.path.join(tmpdirname, 'output')
        if os.path.exists(ace_file):
            shutil.move(ace_file, output)
        if os.path.exists(xsdir_file):
            shutil.move(xsdir_file, output_xsdir)
        if os.path.exists(output_file):
            shutil.move(output_file, output_njoy)

    # Write NJOY output
    open(data.name + '.njo', 'w').write(njoy.stdout)

    return njoy

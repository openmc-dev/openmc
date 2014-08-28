#!/usr/bin/env python2

from __future__ import division, print_function

import sys
import itertools
import re
import warnings

from statepoint import StatePoint

alphanum = re.compile(r"[\W_]+")

err = False

################################################################################
def parse_options():
  """Process command line arguments"""



  def tallies_callback(option, opt, value, parser):
    """Option parser function for list of tallies"""
    global err
    try:
      setattr(parser.values, option.dest, [int(v) for v in value.split(',')])
    except:
      p.print_help()
      err = True

  def scores_callback(option, opt, value, parser):
    """Option parser function for list of scores"""
    global err
    try:
      scores = {}
      entries = value.split(',')
      for e in entries:
        tally,score = [int(i) for i in e.split('.')]
        if not tally in scores: scores[tally] = []
        scores[tally].append(score)
      setattr(parser.values, option.dest, scores)
    except:
      p.print_help()
      err = True

  def filters_callback(option, opt, value, parser):
    """Option parser function for list of filters"""
    global err
    try:
      filters = {}
      entries = value.split(',')
      for e in entries:
        tally,filter_,bin = [i for i in e.split('.')]
        tally,bin = int(tally),int(bin)
        if not tally in filters: filters[tally] = {}
        if not filter_ in filters[tally]: filters[tally][filter_] = []
        filters[tally][filter_].append(bin)
      setattr(parser.values, option.dest, filters)
    except:
      p.print_help()
      err = True

  from optparse import OptionParser
  usage = r"""%prog [options] <statepoint_file>

The default is to process all tallies and all scores into one file. Subsets
can be chosen using the options.  For example, to only process tallies 2 and 4
with all scores on tally 2 and only scores 1 and 3 on tally 4:

%prog -t 2,4 -s 4.1,4.3 <statepoint_file>

Likewise if you have additional filters on a tally you can specify a subset of
bins for each filter for that tally. For example to process all tallies and
scores, but only energyin bin #1 in tally 2:

%prog -f 2.energyin.1 <statepoint_file>

You can list the available tallies, scores, and filters with the -l option:

%prog -l <statepoint_file>"""
  p = OptionParser(usage=usage)
  p.add_option('-t', '--tallies', dest='tallies', type='string', default=None,
                action='callback', callback=tallies_callback,
                help='List of tally indices to process, separated by commas.' \
                     ' Default is to process all tallies.')
  p.add_option('-s', '--scores', dest='scores', type='string', default=None,
               action='callback', callback=scores_callback,
               help='List of score indices to process, separated by commas, ' \
                 'specified as {tallyid}.{scoreid}.' \
                 ' Default is to process all scores in each tally.')
  p.add_option('-f', '--filters', dest='filters', type='string', default=None,
               action='callback', callback=filters_callback,
               help='List of filter bins to process, separated by commas, ' \
                 'specified as {tallyid}.{filter}.{binid}. ' \
                 'Default is to process all filter combinaiton for each score.')
  p.add_option('-l', '--list', dest='list', action='store_true',
               help='List the tally and score indices available in the file.')
  p.add_option('-o', '--output', action='store', dest='output',
               default='tally', help='path to output SILO file.')
  p.add_option('-e', '--error', dest='valerr', default=False,
               action='store_true',
               help='Flag to extract errors instead of values.')
  p.add_option('-v', '--vtk', action='store_true', dest='vtk',
               default=False, help='Flag to convert to VTK instead of SILO.')
  parsed = p.parse_args()

  if not parsed[1]:
    p.print_help()
    return parsed, err

  if parsed[0].valerr:
    parsed[0].valerr = 1
  else:
    parsed[0].valerr = 0

  return parsed, err

################################################################################
def main(file_, o):
  """Main program"""

  sp = StatePoint(file_)
  sp.read_results()

  validate_options(sp, o)

  if o.list:
    print_available(sp)
    return

  if o.vtk:
    if not o.output[-4:] == ".vtm": o.output += ".vtm"
  else:
    if not o.output[-5:] == ".silo": o.output += ".silo"

  if o.vtk:
    try:
      import vtk
    except:
      print('The vtk python bindings do not appear to be installed properly.\n'
            'On Ubuntu: sudo apt-get install python-vtk\n'
            'See: http://www.vtk.org/')
      return
  else:
    try:
      import silomesh
    except:
      print('The silomesh package does not appear to be installed properly.\n'
            'See: https://github.com/nhorelik/silomesh/')
      return

  if o.vtk:
    blocks = vtk.vtkMultiBlockDataSet()
    blocks.SetNumberOfBlocks(5)
    block_idx = 0
  else:
    silomesh.init_silo(o.output)

  # Tally loop #################################################################
  for tally in sp.tallies:

    # skip non-mesh tallies or non-user-specified tallies
    if o.tallies and not tally.id in o.tallies: continue
    if not 'mesh' in tally.filters: continue

    print("Processing Tally {}...".format(tally.id))

    # extract filter options and mesh parameters for this tally
    filtercombos = get_filter_combos(tally)
    meshparms = get_mesh_parms(sp, tally)
    nx,ny,nz = meshparms[:3]
    ll = meshparms[3:6]
    ur = meshparms[6:9]

    if o.vtk:
      ww = [(u-l)/n for u,l,n in zip(ur,ll,(nx,ny,nz))]
      grid = grid = vtk.vtkImageData()
      grid.SetDimensions(nx+1,ny+1,nz+1)
      grid.SetOrigin(*ll)
      grid.SetSpacing(*ww)
    else:
      silomesh.init_mesh('Tally_{}'.format(tally.id), *meshparms)

    # Score loop ###############################################################
    for sid,score in enumerate(tally.scores):

      # skip non-user-specified scrores for this tally
      if o.scores and tally.id in o.scores and not sid in o.scores[tally.id]:
        continue

      # Filter loop ############################################################
      for filterspec in filtercombos:

        # skip non-user-specified filter bins
        skip = False
        if o.filters and tally.id in o.filters:
          for filter_,bin in filterspec[1:]:
            if filter_ in o.filters[tally.id] and \
               not bin in o.filters[tally.id][filter_]:
              skip = True
              break
        if skip: continue

        # find and sanitize the variable name for this score
        varname = get_sanitized_filterspec_name(tally, score, filterspec)
        if o.vtk:
          vtkdata = vtk.vtkDoubleArray()
          vtkdata.SetName(varname)
          dataforvtk = {}
        else:
          silomesh.init_var(varname)

        lbl = "\t Score {}.{} {}:\t\t{}".format(tally.id, sid+1, score, varname)

        # Mesh fill loop #######################################################
        for x in range(1,nx+1):
          sys.stdout.write(lbl+" {0}%\r".format(int(x/nx*100)))
          sys.stdout.flush()
          for y in range(1,ny+1):
            for z in range(1,nz+1):
              filterspec[0][1] = (x,y,z)
              val = sp.get_value(tally.id-1, filterspec, sid)[o.valerr]
              if o.vtk:
                # vtk cells go z, y, x, so we store it now and enter it later
                i = (z-1)*nx*ny + (y-1)*nx + x-1
                dataforvtk[i] = float(val)
              else:
                silomesh.set_value(float(val), x, y, z)

        # end mesh fill loop
        print()
        if o.vtk:
          for i in range(nx*ny*nz):
            vtkdata.InsertNextValue(dataforvtk[i])
          grid.GetCellData().AddArray(vtkdata)
          del vtkdata

        else:
          silomesh.finalize_var()

      # end filter loop

    # end score loop
    if o.vtk:
      blocks.SetBlock(block_idx, grid)
      block_idx += 1
    else:
      silomesh.finalize_mesh()

  # end tally loop
  if o.vtk:
    writer = vtk.vtkXMLMultiBlockDataWriter()
    writer.SetFileName(o.output)
    writer.SetInput(blocks)
    writer.Write()
  else:
    silomesh.finalize_silo()

################################################################################
def get_sanitized_filterspec_name(tally, score, filterspec):
  """Returns a name fit for silo vars for a given filterspec, tally and score"""

  comboname = "_"+" ".join(["{}_{}".format(filter_, bin)
                        for filter_, bin in filterspec[1:]])
  if len(filterspec[1:]) == 0: comboname = ''
  varname = 'Tally_{}_{}{}'.format(tally.id, score, comboname)
  varname = alphanum.sub('_', varname)
  return varname

################################################################################
def get_filter_combos(tally):
  """Returns a list of all filter spec combinations, excluding meshes

  Each combo has the mesh spec as the first element, to be set later.
  These filter specs correspond with the second argument to StatePoint.get_value
  """

  specs = []

  if len(tally.filters) == 1:
    return [[['mesh', [1, 1, 1]]]]

  filters = list(tally.filters.keys())
  filters.pop(filters.index('mesh'))
  nbins = [tally.filters[f].length for f in filters]

  combos = [ [b] for b in range(nbins[0])]
  for i,b in enumerate(nbins[1:]):
    prod = list(itertools.product(combos, range(b)))
    if i == 0:
      combos = prod
    else:
      combos = [[v for v in p[0]] + [p[1]] for p in prod]

  for c in combos:
    spec = [['mesh', [1, 1, 1]]]
    for i,bin in enumerate(c):
      spec.append((filters[i], bin))
    specs.append(spec)

  return specs

################################################################################
def get_mesh_parms(sp, tally):
  meshid = tally.filters['mesh'].bins[0]
  for i,m in enumerate(sp.meshes):
    if m.id == meshid:
      mesh = m
  return mesh.dimension + mesh.lower_left + mesh.upper_right

################################################################################
def print_available(sp):
  """Prints available tallies/scores in a statepoint"""

  print("Available tally and score indices:")
  for tally in sp.tallies:
    mesh = ""
    if not 'mesh' in tally.filters: mesh = "(no mesh)"
    print("\tTally {} {}".format(tally.id, mesh))
    scores = ["{}.{}: {}".format(tally.id, sid, score)
              for sid, score in enumerate(tally.scores)]
    for score in scores:
      print("\t\tScore {}".format(score))
      for filter_ in tally.filters:
        if filter_ == 'mesh': continue
        for bin in range(tally.filters[filter_].length):
          print("\t\t\tFilters: {}.{}.{}".format(tally.id, filter_, bin))

################################################################################
def validate_options(sp,o):
  """Validates specified tally/score options for the current statepoint"""

  available_tallies = [t.id for t in sp.tallies]
  if o.tallies:
    for otally in o.tallies:
      if not otally in available_tallies:
        warnings.warn('Tally {} not in statepoint file'.format(otally))
        continue
      else:
        for tally in sp.tallies:
          if tally.id == otally: break
      if not 'mesh' in tally.filters:
        warnings.warn('Tally {} contains no mesh'.format(otally))
      if o.scores and otally in o.scores.keys():
        for oscore in o.scores[otally]:
          if oscore > len(tally.scores):
            warnings.warn('No score {} in tally {}'.format(oscore, otally))

  if o.scores:
    for otally in o.scores.keys():
      if not otally in available_tallies:
        warnings.warn('Tally {} not in statepoint file'.format(otally))
        continue
      if o.tallies and not otally in o.tallies:
        warnings.warn(
          'Skipping scores for tally {}, excluded by tally list'.format(otally))
        continue

  if o.filters:
    for otally in o.filters.keys():
      if not otally in available_tallies:
        warnings.warn('Tally {} not in statepoint file'.format(otally))
        continue
      if o.tallies and not otally in o.tallies:
        warnings.warn(
          'Skipping filters for tally {}, excluded by tally list'.format(otally))
        continue
      for tally in sp.tallies:
        if tally.id == otally: break
      for filter_ in o.filters[otally]:
        if filter_ == 'mesh':
          warnings.warn('Cannot specify mesh filter bins')
          continue
        if not filter_ in tally.filters.keys():
          warnings.warn(
                  'Tally {} does not contain filter {}'.format(otally, filter_))
          continue
        for bin in o.filters[otally][filter_]:
          if bin >= tally.filters[filter_].length:
            warnings.warn(
                 'No bin {} in tally {} filter {}'.format(bin, otally, filter_))

################################################################################
# monkeypatch to suppress the source echo produced by warnings
def formatwarning(message, category, filename, lineno, line):
  return "{}:{}: {}: {}\n".format(filename, lineno, category.__name__, message)
warnings.formatwarning = formatwarning

################################################################################
if __name__ == '__main__':
    (options, args), err = parse_options()
    if args and not err:
      main(args[0],options)

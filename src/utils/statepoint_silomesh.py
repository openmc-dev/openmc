#!/usr/bin/env python 

import warnings

from statepoint import StatePoint

import silomesh # https://github.com/nhorelik/silomesh/

def parse_options():
  """Process command line arguments"""

    def tallies_callback(option, opt, value, parser):
      """Option parser function for list of tallies"""
      try:
        setattr(parser.values, option.dest, [int(v) for v in value.split(',')])
      except:
        p.print_help()

    def scores_callback(option, opt, value, parser):
      """Option parser function for list of scores"""
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

    from optparse import OptionParser
    usage = r"""%prog [options] <statepoint_file>

The default is to process all tallies and all scores into one silo file. Subsets
can be chosen using the options.  For example, to only process tallies 2 and 4
with all scores on tally 2 and only scores 1 and 3 on tally 4:

%prog -t 2,4 -f 4.1,4.3 <statepoint_file>

You can list the available tallies and scores with the -l option:

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
    p.add_option('-l', '--list', dest='list', action='store_true',
              help='List the tally and score indices available in the file.')
    parsed = p.parse_args()
    
    if not parsed[1]: p.print_help()

    return parsed

################################################################################
def main(file_,o):
  """Main program"""
  
  sp = StatePoint(file_)
  sp.read_results()

  validate_options(sp,o)
  
  if o.list:
    print_available(sp)
    return

  for tally in sp.tallies:
    if o.tallies and not tally.id in o.tallies: continue
    if not 'mesh' in tally.filters: continue
    nx,ny,nz = get_mesh_dims(sp, tally)
    

  print o.tallies
  print o.scores

################################################################################

def get_mesh_dims(sp, tally):
  meshid = tally.filters['mesh'].bins[0]
  for i,m in enumerate(sp.meshes):
    if m.id == meshid:
      mesh = m
  return mesh.dimension

def print_available(sp):
  """Prints available tallies/scores in a statepoint"""
  
  print "Available tally and score indices:"
  for tally in sp.tallies:
    mesh = ""
    if not 'mesh' in tally.filters: mesh = "(no mesh)"
    print "\tTally {} {}".format(tally.id,mesh)
    scores = ["{}.{}: {}".format(tally.id, sid+1, score) 
              for sid,score in enumerate(tally.scores)]
    for score in scores:
      print "\t\tScore {}".format(score)

def validate_options(sp,o):
  """Validates specified tally/score options for the current statepoint"""
  
  if o.tallies:
    available_tallies = [t.id for t in sp.tallies]
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
            warnings.warn('No score {} in tally {}'.format(oscore,otally))
              
  if o.scores:
    for otally in o.scores.keys():
      if not otally in available_tallies:
        warnings.warn('Tally {} not in statepoint file'.format(otally))
        continue
  
# monkeypatch to suppress the source echo produced by warnings
def formatwarning(message, category, filename, lineno, line):
  return "%s:%s: %s: %s\n" % (filename, lineno, category.__name__, message)
warnings.formatwarning = formatwarning

if __name__ == '__main__':

    (options, args) = parse_options()
    
    if args:
      main(args[0],options)

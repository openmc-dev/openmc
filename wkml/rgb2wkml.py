#!/usr/bin/python

# A script to read in the X11 color file, and create a fortran module exporting the defined colors.
# This will fail on any file with more than (maxContLines-1)^2 colors.

from string import Template

maxContLines=39 # in f95

header = """module m_wkml_color_def

  implicit none
  private

  type color
    integer       :: r, g, b
    character(20) :: name
  end type color

  public :: color
  public :: colorarray
"""
c = Template("color($r, $g, $b, '$name')")

fIn = open("rgb.txt")

fIn.readline() # throw away the first line, which is the CVS tag.

colors = []
for line in fIn:
    values = line.split()
    (r, g, b) = [int(x) for x in values[:3]]
    name = ' '.join(values[3:])
    if len(name)>20: continue
    name = name+(20-len(name))*' '
    colors.append({'r':r, 'g':g, 'b':b, 'name':name})
    
numArraysNeeded = len(colors)/(maxContLines) + 1

print header

ii = 0
for i in range(numArraysNeeded-1):
    print "  type(color), parameter :: colorArray%d(%d)" % (i, maxContLines) , 
    print "= (/ &" 
    for j in range(maxContLines-1):
        print "    "+c.substitute(colors[ii])+', &'
        ii += 1
    print "    "+c.substitute(colors[ii])+'/)'
    ii += 1
    print

print
print "  type(color), parameter :: colorArray%d(%d)" % (numArraysNeeded-1, len(colors)-ii) , 
print "= (/ &"
for j in range(len(colors)-ii-1):
    print "    "+c.substitute(colors[ii])+', &'
    ii += 1
print "    "+c.substitute(colors[ii])+'/)'
print

print "  type(color), parameter :: colorArray(%d)" % len(colors), 
print "= (/ &"
for i in range(numArraysNeeded-1):
    print "    colorArray%d" % i+", &"
print "    colorArray%d" % (numArraysNeeded-1), 
print "/)"

print "end module m_wkml_color_def"





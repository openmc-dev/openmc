#!/usr/bin/env python

# Package to read OpenMC input XML files for use in Python scripts
# Right now functionality for tallies.xml is all that is supported.

# Only limited xml checking is performed as it is assumed this is to be used
# for post-processing, and thus the inputs have been successfully run.

from lxml import objectify
import os.path

class tally:
    def __init__(self,node,default_id):
        # Need to scan both the contained objects and attributes for info
        # Will first scan attrib, then objects

        # id:
        self.id = node.attrib.get("id") # Will be None if it doesnt exist
        if self.id == None: # Scan the objects.
            self.id = node.find("id")
            if self.id == None: # It wasnt found, i'll have to use default_id
                self.id = default_id
            else: # It was found, but convert to int.
                self.id = int(self.id)
        else: # Convert self.id to a python int
            self.id = int(self.id['id']) # Getting the dictionary's value as int
        
        # label:
        self.label = node.attrib.get("label") # Will be None if it doesnt exist
        if self.label == None: # Scan the objects.
            self.label = node.find("label")
            if self.label == None: # It wasnt found, use the id
                self.label = str(self.id)
            else:
                self.label = self.label.text
        else: # Convert self.label to a python str
            self.label = self.label['label'].text
            
        # estimator:
        self.estimator = node.attrib.get("estimator") # Will be None if it doesnt exist
        if self.estimator == None: # Scan the objects.
            self.estimator = node.find("estimator")
            if self.estimator == None: # It wasnt found, Set to Default
                # to denote that we are not over-riding the score's default.
                self.estimator = 'Default'
            else:
                self.estimator = self.estimator.text
        else: # Convert self.estimator to a python str
            self.estimator = self.estimator['estimator'].text
            
        # scores:
        scoresList = node.attrib.get("scores") # Will be None if it doesnt exist
        if scoresList == None: # Scan the objects.
            scoresList = node.find("scores")
            if scoresList == None: # It wasnt found, Odd, print error, store
                # as "None"
                scoresList = "None"
                print "Error: Empty <scores> in tally ", self.id
            else:
                scoresList = scoresList.text
        else: # Convert self.estimator to a python str
            scoresList = scoresList['scores'].text
        # Now split scoresList and store as self.scores            
        self.scores = scoresList.split()
        
        # nuclides:
        nuclidesList = node.attrib.get("nuclides") # Will be None if it doesnt exist
        if nuclidesList == None: # Scan the objects.
            nuclidesList = node.find("nuclides")
            if nuclidesList == None: # It wasnt found, store the default, 
                nuclidesList = "total"
            else:
                nuclidesList = nuclidesList.text
        else: # Convert self.estimator to a python str
            nuclidesList = nuclidesList['nuclides'].text
        # Now split scoresList and store as self.scores            
        self.nuclides = nuclidesList.split()
        
        # Get filter. Filters are assumed to NOT be present in the attribute.
        filterNode = node.find("filters")
        if filterNode == None: # Uhoh. Error,
            print "Error: Empty <filters> in tally ", self.id
        else:
            # Pass to the filters constructor.
            self.filters = tally_filters(filterNode)
        
            
class tally_filters:
    def __init__(self,node):
        self.cell = 0
        
class tally_mesh:
    def __init__(self,node, default_id):
        self.id = 0

# The talliesXML class holds all the information in tallies.xml and provides 
# interfaces for getting useful information out of them

class talliesXML(object):
    def __init__(self,filename):
        # Initialize tallies and meshes to a list.
        self.tallies = []
        self.meshes = []
        
        # Check if filename exists
        if os.path.isfile(filename):
            # Get the data from tallies.xml            
            tree = objectify.parse(filename)
            root = tree.getroot()
            
            # Find and Get value of assume_seperate
            separate_val = root.find("assume_separate")
            # If the user did not specify assume_separate, then the default
            # shall be used (False)
            # Assume_separate can possibly exist as an attrib of the root.
            # Therefore we must check that true.
            # Other items, due to their length, likely will not exist as an 
            # attrib.
            if separate_val == None:
                # Now check the attrib
                separate_val = root.attrib.get("assume_separate")
                if separate_val == None: # Then we still found nothing
                    self.assume_separate = False
                else: # got it, it is a tag.
                    # Convert to a string, and for easy comparison capitalize
                    # and strip the white space
                    separate_val = separate_val.text.upper().lstrip().rstrip()
                    if separate_val == 'NO':
                        self.assume_separate = False
                    elif separate_val == 'YES':
                        self.assume_separate = True
            else: 
                # Convert to a string, and for easy comparison capitalize it
                # and strip the white space
                separate_val = separate_val.text.upper()
                separate_val = separate_val.lstrip().rstrip()
                if separate_val == 'NO':
                    self.assume_separate = False
                elif separate_val == 'YES':
                    self.assume_separate = True
            
            # Perform tally processing.
            i = 0
            for tallyNode in root.findall("tally"):
                self.tallies.append(tally(tallyNode, i))
                i = i + 1
                
            # Perform mesh processing.
            i = 0
            for meshNode in root.findall("mesh"):
                self.meshes.append(tally_mesh(meshNode, i))
                i = i + 1
            
        else:
            print "Error, ", filename, " does not exist. " \
                "Returning null talliesXML object."
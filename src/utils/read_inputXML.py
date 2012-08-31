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
        else: # Convert scoresList to a python str
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
        else: # Convert nuclidesList to a python str
            nuclidesList = nuclidesList['nuclides'].text
        # Now split nuclidesList and store as self.scores            
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
        # cells:
        cellList = node.attrib.get("cell") # Will be None if it doesnt exist
        if cellList == None: # Scan the objects.
            cellList = node.find("cell")
            if cellList != None: # It wasnt found, 
                cellList = cellList.text
        else: # Convert cellList to a python string, for later conversion to int
            cellList = cellList['cell'].text
        # Now split cellList and store as self.cell            
        if cellList != None:        
            self.cell = cellList.split()
            # convert to int
            for i in range(len(self.cell)):
                self.cell[i] = int(self.cell[i])
        else:
            self.cell = [None]
            
        # surfaces:
        surfList = node.attrib.get("surf") # Will be None if it doesnt exist
        if surfList  == None: # Scan the objects.
            surfList  = node.find("surf")
            if surfList  != None: # It wasnt found, 
                surfList  = surfList .text
        else: # Convert surfList to a python string, for later conversion to int
            surfList  = surfList ['surf'].text
        # Now split surfList and store as self.surf            
        if surfList  != None:        
            self.surf = surfList .split()
            # convert to int
            for i in range(len(self.surf)):
                self.surf[i] = int(self.surf[i])
        else:
            self.surf = [None]
            
        # universe:
        universeList = node.attrib.get("universe") # Will be None if it doesnt exist
        if universeList  == None: # Scan the objects.
            universeList  = node.find("universe")
            if universeList  != None: # It wasnt found, 
                universeList  = universeList .text
        else: # Convert universeList to a python string, for later conversion to int
            universeList  = universeList ['universe'].text
        # Now split universeList and store as self.universe            
        if universeList  != None:        
            self.universe = universeList .split()
            # convert to int
            for i in range(len(self.universe)):
                self.universe[i] = int(self.universe[i])
        else:
            self.universe = [None]
            
        # material:
        materialList = node.attrib.get("material") # Will be None if it doesnt exist
        if materialList  == None: # Scan the objects.
            materialList  = node.find("material")
            if materialList  != None: # It wasnt found, 
                materialList  = materialList .text
        else: # Convert materialList to a python string, for later conversion to int
            materialList  = materialList ['material'].text
        # Now split materialList and store as self.material            
        if materialList  != None:        
            self.material = materialList .split()
            # convert to int
            for i in range(len(self.material)):
                self.material[i] = int(self.material[i])
        else:
            self.material = [None]
            
        # mesh:
        self.mesh = node.attrib.get("mesh") # Will be None if it doesnt exist
        if self.mesh == None: # Scan the objects.
            self.mesh = node.find("mesh")
            if self.mesh != None: # It was found, but convert to int.
                self.mesh = int(self.mesh)
        else: # Convert self.mesh to a python int
            self.mesh = int(self.mesh['mesh']) # Getting the dictionary's value as int
            
        # cellborn:
        cellbornList = node.attrib.get("cellborn") # Will be None if it doesnt exist
        if cellbornList  == None: # Scan the objects.
            cellbornList  = node.find("cellborn")
            if cellbornList  != None: # It wasnt found, 
                cellbornList  = cellbornList .text
        else: # Convert cellbornList to a python string, for later conversion to int
            cellbornList  = cellbornList['cellborn'].text
        # Now split cellbornList and store as self.cellborn            
        if cellbornList  != None:        
            self.cellborn = cellbornList .split()
            # convert to int
            for i in range(len(self.cellborn)):
                self.cellborn[i] = int(self.cellborn[i])
        else:
            self.cellborn = [None]
        
        # energy:
        energyList = node.attrib.get("energy") # Will be None if it doesnt exist
        if energyList == None: # Scan the objects.
            energyList = node.find("energy")
            if energyList != None: # It wasnt found, 
                energyList = energyList.text
        else: # Convert energyList to a python string, for later conversion to int
            energyList = energyList['energy'].text
        # Now split energyList and store as self.energy            
        if energyList != None:        
            self.energy = energyList.split()
            # convert to double
            for i in range(len(self.energy)):
                self.energy[i] = float(self.energy[i])
        else:
            self.energy = [None]
            
        # energyout:
        energyoutList = node.attrib.get("energyout") # Will be None if it doesnt exist
        if energyoutList == None: # Scan the objects.
            energyoutList = node.find("energyout")
            if energyoutList != None: # It wasnt found, 
                energyoutList = energyoutList.text
        else: # Convert energyoutList to a python string, for later conversion to int
            energyoutList = energyoutList['energyout'].text
        # Now split energyoutList and store as self.energyout            
        if energyoutList != None:        
            self.energyout = energyoutList.split()
            # convert to double
            for i in range(len(self.energyout)):
                self.energyout[i] = float(self.energyout[i])
        else:
            self.energyout = [None]
        
class tally_mesh:
    def __init__(self, node, default_id):
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
            
        # type:
        self.type = node.attrib.get("type") # Will be None if it doesnt exist
        if self.type == None: # Scan the objects.
            self.type = node.find("type")
            if self.type != None: 
                self.type = self.type.text
        else: # Convert self.type to a python str
            self.type = self.type['type'].text
            
        # dimension:
        dimensionList = node.attrib.get("dimension") # Will be None if it doesnt exist
        if dimensionList == None: # Scan the objects.
            dimensionList = node.find("dimension")
            if dimensionList != None: # It wasnt found, 
                dimensionList = dimensionList.text
        else: # Convert dimensionList to a python string, for later conversion to int
            dimensionList = dimensionList['dimension'].text
        # Now split dimensionList and store as self.dimension            
        if dimensionList != None:        
            self.dimension = dimensionList.split()
            # convert to int
            for i in range(len(self.dimension)):
                self.dimension[i] = int(self.dimension[i])
        else:
            self.dimension = [None]
            
        # lower_left:
        lower_leftList = node.attrib.get("lower_left") # Will be None if it doesnt exist
        if lower_leftList == None: # Scan the objects.
            lower_leftList = node.find("lower_left")
            if lower_leftList != None: # It wasnt found, 
                lower_leftList = lower_leftList.text
        else: # Convert lower_leftList to a python string, for conversion to float
            lower_leftList = lower_leftList['lower_left'].text
        # Now split lower_leftList and store as self.lower_left            
        if lower_leftList != None:        
            self.lower_left = lower_leftList.split()
            # convert to float
            for i in range(len(self.lower_left)):
                self.lower_left[i] = float(self.lower_left[i])
        else:
            self.lower_left = [None]
            
        # upper_right:
        upper_rightList = node.attrib.get("upper_right") # Will be None if it doesnt exist
        if upper_rightList == None: # Scan the objects.
            upper_rightList = node.find("upper_right")
            if upper_rightList != None: # It wasnt found, 
                upper_rightList = upper_rightList.text
        else: # Convert upper_rightList to a python string, for conversion to float
            upper_rightList = upper_rightList['upper_right'].text
        # Now split upper_rightList and store as self.upper_right            
        if upper_rightList != None:        
            self.upper_right = upper_rightList.split()
            # convert to float
            for i in range(len(self.upper_right)):
                self.upper_right[i] = float(self.upper_right[i])
        else:
            self.upper_right = [None]
            
        # width:
        widthList = node.attrib.get("width") # Will be None if it doesnt exist
        if widthList == None: # Scan the objects.
            widthList = node.find("width")
            if widthList != None: # It wasnt found, 
                widthList = widthList.text
        else: # Convert widthList to a python string, for conversion to float
            widthList = widthList['width'].text
        # Now split widthList and store as self.width            
        if widthList != None:        
            self.width = widthList.split()
            # convert to float
            for i in range(len(self.width)):
                self.width[i] = float(self.width[i])
        else:
            self.width = [None]

# The talliesXML class holds all the information in tallies.xml and provides 
# interfaces for getting useful information out of them

class talliesXML:
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
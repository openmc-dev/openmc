import os
import copy
import openmc
import openmc.deplete
import openmc.data
import matplotlib.pyplot as plt

######################################################################################
#           Create a chain that includes spont. fission by reading from endf files.
######################################################################################

#initiate an instance of the Chain class
mychain = openmc.deplete.Chain()

#Find all relevant endf files
decay_files = []
decay_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/decay/"
for file in os.listdir(decay_folder):
    if file[-4:] == "endf":
        decay_files.append(decay_folder + file)


nfy_files = []
nfy_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/nfy/"
for file in os.listdir(nfy_folder):
    if file[-4:] == "endf":
        nfy_files.append(nfy_folder + file)


xs_files = []
xs_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/neutrons/"
for file in os.listdir(xs_folder):
    if file[-4:] == "endf":
        xs_files.append(xs_folder + file)

sfy_files = []
sfy_folder = "/home/shauksson/evaluations/ENDF-B-VIII.0/sfy/"
for file in os.listdir(sfy_folder):
    if file[-4:] == "endf":
        sfy_files.append(sfy_folder + file)

#Create the chain by reading from endf files
mychain = mychain.from_endf(decay_files,nfy_files,xs_files,sfy_files)


#Check whether information in Nuclide classes included in mychain agree with the attribute spont_fission_yields
agree = True
i = mychain.nuclide_dict['U238']
for product in mychain.spont_fission_yields[0]['U238'].keys():
    agree = agree and mychain.spont_fission_yields[0]['U238'][product] == mychain.nuclides[i].spont_yield_data[0][product]
print("Does information in Nuclide classes agree with spont_fission_yields?",agree)
print("")

#Example of contents of one instance of Nuclide in the chain.
print("Attributes of nuclide U238")
print("name",mychain.nuclides[i].name)
print("half_life",mychain.nuclides[i].half_life)
print("decay_energy",mychain.nuclides[i].decay_energy)
print("n_decay_modes",mychain.nuclides[i].n_decay_modes)
print("decay_modes",mychain.nuclides[i].decay_modes)
print("decay_modes, example",mychain.nuclides[i].decay_modes[0].type)
print("n_reaction_paths",mychain.nuclides[i].n_reaction_paths)
print("reactions",mychain.nuclides[i].reactions)
print("sources",mychain.nuclides[i].sources)
print("yield_data",mychain.nuclides[i].yield_data)
print("spont_yield_data",mychain.nuclides[i].spont_yield_data)
print("yield_energies",mychain.nuclides[i].yield_energies)
print("")



######################################################################################
#           Can read and write xml files with information about chains including spont. fission
######################################################################################

#Can print chain to xml file
mychain.export_to_xml("./chain.xml")

#Can then read the chain from a xml file. Gives correct spontaneous fission.
newchain = openmc.deplete.Chain.from_xml("./chain.xml")
new_U238index = newchain.nuclide_dict['U238']
print("Was spontanteous fission read correctly from xml file?",newchain.nuclides[new_U238index].spont_yield_data ==  mychain.nuclides[i].spont_yield_data)
print("")


######################################################################################
#           Can check consistency and reduce chains including spont. fission
######################################################################################


#Checks whether chain is consistent, including spontaneous fission
print("Chain is consistent: ",mychain.validate())
#Have changed spontaneous fission yields manually, so that spontaneous fission yields of U238 do not sum to 2.0.
wrongchain = openmc.deplete.Chain.from_xml("./chain_wrong.xml")
print("Expect to get false when validating incorrect chain:",wrongchain.validate(strict=False))


#Create reduced chain 
reduced_chain = mychain.reduce(["U238"],level=1)
#Check that spont. fission has been included correctly.
reduced_U238index = reduced_chain.nuclide_dict['U238']
print("sf is copied correctly by reduce: ",reduced_chain.nuclides[reduced_U238index].spont_yield_data == mychain.nuclides[i].spont_yield_data)



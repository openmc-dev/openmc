import os
import openmc
import openmc.deplete
import openmc.data
import matplotlib.pyplot as plt

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

#read from endf file
mychain = mychain.from_endf(decay_files,nfy_files,xs_files,sfy_files)


#try printing out the chain, especially spontaneous fission variables
#print(mychain)
#print(mychain.nuclides)
#print(mychain.reactions)
#print(mychain.nuclide_dict)
#print(mychain.stable_nuclides)
#print(mychain.unstable_nuclides)
print(mychain.fission_yields[0].keys())
print(mychain.spont_fission_yields[0].keys())
print(mychain.spont_fission_yields[0]['U238'])
print(sum(mychain.spont_fission_yields[0]['U238'].values()))


#Print chain to xml file
mychain.export_to_xml("./chain.xml")

#Read in the chain from the xml file and print to a new xml file.
newchain = openmc.deplete.Chain.from_xml("./chain.xml")
newchain.export_to_xml("./chain_new.xml")




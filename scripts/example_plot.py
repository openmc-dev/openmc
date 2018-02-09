"""An example file showing how to plot data from a simulation."""

import matplotlib.pyplot as plt

from opendeplete import read_results, \
                        evaluate_single_nuclide, \
                        evaluate_reaction_rate, \
                        evaluate_eigenvalue

# Set variables for where the data is, and what we want to read out.
result_folder = "test"

# Load data
results = read_results(result_folder + "/results.h5")

cell = "5"
nuc = "Gd157"
rxn = "(n,gamma)"

# Total number of nuclides
plt.figure()
# Pointwise data
x, y = evaluate_single_nuclide(results, cell, nuc)
plt.semilogy(x, y)

plt.xlabel("Time, s")
plt.ylabel("Total Number")
plt.savefig("number.pdf")

# Reaction rate
plt.figure()
x, y = evaluate_reaction_rate(results, cell, nuc, rxn)
plt.plot(x, y)
plt.xlabel("Time, s")
plt.ylabel("Reaction Rate, 1/s")

plt.savefig("rate.pdf")

# Eigenvalue
plt.figure()
x, y = evaluate_eigenvalue(results)
plt.plot(x, y)
plt.xlabel("Time, s")
plt.ylabel("Eigenvalue")

plt.savefig("eigvl.pdf")

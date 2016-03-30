import openmc


su = openmc.Summary('summary.h5')
sp = openmc.StatePoint('statepoint.1.h5')
sp.link_with_summary(su)
print(sp.tallies[1].get_pandas_dataframe(summary=su))

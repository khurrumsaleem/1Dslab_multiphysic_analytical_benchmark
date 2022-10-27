import openmc

sp = openmc.StatePoint('statepoint.500.h5')
flux_tally = sp.get_tally(scores=['flux'])
kappa_fission_tally = sp.get_tally(scores=['kappa-fission'])

print(flux_tally)
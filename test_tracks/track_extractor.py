import openmc
tracks = openmc.Tracks('tracks.h5')
num_tracks = len(tracks)

for track in tracks:
  # get number of tracks for each particle and use them below
  num_particle_tracks = len(track.particle_tracks)
  for i in range(0,num_particle_tracks):
    neutron = track.particle_tracks[i]
    print(neutron.states['u'])
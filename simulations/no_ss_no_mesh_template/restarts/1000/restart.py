import openmc
import openmc.lib

openmc.lib.init(output=True)
openmc.lib.import_properties("properties.h5")
openmc.lib.run(threads=16)
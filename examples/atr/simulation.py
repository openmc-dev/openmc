import openmc

settings = openmc.Settings()
settings.batches = 30
settings.inactive = 10
settings.particles = 4096
settings.export_to_xml()

openmc.run()
#!/usr/bin/env python

import openmc
import numpy as np
import matplotlib.pyplot as plt


######################    Material Definition   #################################### 

fuel = openmc.Material(name = 'uo2')
fuel.add_nuclide('U235', 3.8393e-5)
fuel.add_nuclide('U238', 1.8917e-2)
fuel.add_nuclide('O16', 4.1707e-2)
fuel.add_nuclide('Pu239', 6.5875e-4)
fuel.add_nuclide('Pu240', 4.2323e-5)
fuel.add_nuclide('Pu241', 7.0246e-6)
fuel.set_density('sum')
fuel.temperature = 300

zirconium = openmc.Material(name='clad')
zirconium.add_element('Zr', 4.2300e-2)
zirconium.set_density('sum')
zirconium.temperature = 300

water = openmc.Material(name='Borated')
water.add_nuclide('H1', 6.694e-2)
water.add_nuclide('O16', 3.347e-2)
water.add_nuclide('B10', 6.6262e-6)
water.add_nuclide('B11', 2.6839e-5)
water.add_s_alpha_beta('c_H_in_H2O')
water.set_density('sum')
water.temperature = 323.6

mat = openmc.Materials([fuel, water, zirconium])
mat.export_to_xml()

###########################   Geometry       #################################

fuel_or = openmc.ZCylinder(x0=0.0, y0=0.0, r=0.386)
clad_ir = openmc.ZCylinder(x0=0.0, y0=0.0, r=0.45)

min_z = openmc.ZPlane(z0=-137.413, boundary_type= 'vacuum')
max_z = openmc.ZPlane(z0=+137.413, boundary_type= 'vacuum')


pin_cell_universe = openmc.Universe(name='Fuel Pin')

fuel_re = openmc.Cell(name='fuel region', region = -fuel_or, fill = fuel)
pin_cell_universe.add_cell(fuel_re)

clad_re = openmc.Cell(name='clad region',region = +fuel_or & -clad_ir, fill = zirconium)
pin_cell_universe.add_cell(clad_re)

water_re = openmc.Cell(name='water region',region = +clad_ir, fill = water)
pin_cell_universe.add_cell(water_re)


root_cell = openmc.Cell(name='root cell')
root_cell.fill = pin_cell_universe

hexa = openmc.model.hexagonal_prism(edge_length=0.6375, orientation='y', origin=(0.0, 0.0), boundary_type='reflective', corner_radius=0.)
root_cell.region =  +min_z & -max_z & hexa

# Create root Universe

root_universe = openmc.Universe(universe_id=0, name='root universe')
root_universe.add_cell(root_cell)

# Create Geometry and set root Universe

geometry = openmc.Geometry(root_universe)

geometry.export_to_xml()


plot = openmc.Plot(plot_id=1)
plot.origin = [0, 0, 0]
plot.width  = [1.3, 1.3]
plot.pixels = [500, 500]
plot.color_by = 'material'
# Instantiate a Plots object and export to XML
plot_file = openmc.Plots([plot])
plot_file.export_to_xml()


settings_file = openmc.Settings()
settings_file.batches = 100
settings_file.inactive = 20
settings_file.particles = 1000
settings_file.export_to_xml()



# Instantiate an empty Tallies
tallies_file = openmc.Tallies()

# Create mesh tally to score flux and fission rate (defaults to analog estimator because of nu-scatter)
tally = openmc.Tally(name='fuel rxn rates')
tally.scores = ['nu-fission', 'nu-scatter', 'total']
tallies_file.append(tally)
# Export to "tallies.xml"
tallies_file.export_to_xml()

# Load the statepoint file
sp = openmc.StatePoint('statepoint.100.h5')

# Extract eigenvalue
k_openmc = sp.k_combined

print(k_openmc)

# Get tally from file
tally = sp.get_tally(name='fuel rxn rates')


# Extract tallies

nu_fission_rate = tally.get_slice(scores=['nu-fission'])
total_rate = tally.get_slice(scores=['total'])
nu_scatter_rate = tally.get_slice(scores=['nu-scatter'])

# calculate k-inf
kinf = nu_fission_rate / (total_rate - nu_scatter_rate)

print ("k-inf is {}".format(kinf.mean[0,0,0]))


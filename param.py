
# Internal Paraboloid surface is a*x**2 + b*y**2 - c*z = 0
# with 1/a = 15000**2, 1/b=15000^2, and 1/c=50000 (numbers are microns)
# The total thickness from inner to outer surface is "nominal_thickness"
# There are n_layer and the distance between layers is layer_thickness.
# Layers are normal to the surface along the gradient relative to the
# inner surface (only the inner surface is a true paraboloid):
# (2*a*x, 2*b*y, -c)
# Regions (or abstract cells) are space filling and defined
# by the location of the corner with the least z, least layer, least angle.
# I.e. not the center of the region.
# They are laid out end to end around circles of constant z within a layer.
# The number of regions around a circle is an integer and so the length
# of each region is as close to "nominal_region_length" as possible.
# Note that circles with greater z and farther from the inner surface have
# longer circumferences and generally, therefore, have more regions.
# Circles within a layer are separated by "layer_surface_circle_distance
# The inside layer circle closest to the tip has radius of 500 um.

# Connections along a circle are obvious. Regions of adjacent circles within
# a layer overlap in the lenght direction. Regions of adjacent layers
# overlap in both the length and adjacent circle dimension.

from common import pr

cell_length = 100.
cell_diameter = 30.

# Overall macroscopic shape
nominal_height = 50000.
nominal_base_radius = 15000.
nominal_thickness = 5000.
hole_radius = 500.

# Region discretization
cellbased = False
layer_surface_circle_distance = cell_diameter if cellbased else 500.
layer_thickness = cell_diameter if cellbased else 500.
nominal_region_length = cell_length if cellbased else 500.

# Simulate region
simulation_center = (15000., 0., 50000.) # (x,y,z) coordinate
simulation_region = 2000. # Only simulate cells closer to simulation_center

n_layer = int(nominal_thickness/layer_thickness)

abc = (1/nominal_base_radius**2, 1/nominal_base_radius**2, 1/nominal_height)

pr("layer_surface_circle_distance %g" % layer_surface_circle_distance)
pr("layer_thickness %g" % layer_thickness)
pr("nominal_region_length %g" % nominal_region_length)


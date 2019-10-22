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

cell_length = 100. #um
cell_diameter = 30. #um

# Overall macroscopic shape
nominal_height = 50000.
nominal_base_radius = 15000.
nominal_thickness = 30. #100. # three layers. Should be 5000.
hole_radius = 500.

# Total gap junction of standard size cell is 30nS
# Gap junction conductance is disributed evenly over the surface.
stdcellarea = 2*cell_diameter**2 + 4*cell_diameter*cell_length #um2
meang = 1500 # (nS) over standard cell area
purkinje_gap_factor = 20 # this much larger than gap between regular cells.

# see purkspec.py about list of "rectangles" defined as
# ((angle, n), (bottom, top))
purkinje_spec = [
  ((0, 3), (0, -1)), # three cells wide from bottom to top at angle 0
    ((0, 120), (1200, 1202)), # npts[0][1200] = 768
    ((0, 120), (800, 802)), # know that the npts[0][800] = 617
    ((0, 100), (500, 502)), # know that npts[0][500] = 476
  ((180, 3), (0, -300)), # vertical line at opposite side
    ((180, 200), (800, 802)), # know that the npts[0][800] = 617
    ((180, 100), (500, 502)), # know that npts[0][500] = 476
]

# Region discretization
cellbased = True
layer_surface_circle_distance = cell_diameter if cellbased else 500.
layer_thickness = cell_diameter if cellbased else 500.
nominal_region_length = cell_length if cellbased else 500.

n_layer = int(nominal_thickness/layer_thickness)

# Simulate region
simulation_center = (15000., 0., 50000.) # (x,y,z) coordinate
simulation_region = 100000. # Only simulate cells closer to simulation_center

simulation_layers = (0, 1) # just 0 layer
simulation_circles = (0, None) # all circles
simulation_angledeg = (0., 90.) # quarter of the circumference

abc = (1/nominal_base_radius**2, 1/nominal_base_radius**2, 1/nominal_height)

#pr("layer_surface_circle_distance %g" % layer_surface_circle_distance)
#pr("layer_thickness %g" % layer_thickness)
#pr("nominal_region_length %g" % nominal_region_length)

def print_param():
  import sys
  types = [type(True), type(1), type(1.0), type((1,2,3))]
  try:
    p = sys.modules['param']
    for name in dir(p):
      if '__' not in name:
        val = getattr(p, name)
        if type(val) in types or "purkinje_spec" in name:
          pr("%s = %s"%(name, str(val)))
  except:
    print ('Error in param.print_param')
    pass

if __name__ != '__main__':
  print_param()


# Internal Paraboloid surface is a*x**2 + b*y**2 - c*z = 0
# with 1/a = 15000**2, 1/b=15000^2, and 1/c=50000 (numbers are microns)
# The thickness is 0.5cm 
# There are 5 layers moving outward separated by 5000/4 um.
# Layers are normal to the surface along the gradient:
# (2*a*x, 2*b*y, -c)
# Cells are nominally L=100 and diam=30 and are arranged end to end
# around circles of constant z. The number of cells on an interior
# circle is the greatest integer such that L >=100. All layer circles
# associated with the same z, have the same number of cells (length
# becomes slightly larger as layer index increases).
# Circles on the internal surface are
# separated by 300 um normal to the gradient. The inside layer circle closest
# to the tip has radius of 500 um.

# Connections between layers and along a circle are obvious.

nominal_height = 50000.
nominal_base_radius = 15000.
thickness = 5000.
hole_radius = 500.
internal_surface_circle_distance = 300.
n_layer = 5
nominal_cell_length = 100.
cell_diameter = 30.

abc = (1/nominal_base_radius**2, 1/nominal_base_radius**2, 1/nominal_height)


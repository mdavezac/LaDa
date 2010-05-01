from boost.mpi import world
from lada.vff import Vff
from lada.crystal import Structure
vff = Vff()
vff.lattice.set_types = ("In", "Ga"), ("As",)
vff.lattice.scale = 6.5
vff.add_bond = "In", "As", (2.62332, 21.6739, -112.0, 150.0)
vff.add_bond = "Ga", "As", (2.44795, 32.1530, -105.0, 150.0)
vff.add_angle = "As", "Ga", "As", ("tet", -4.099, 9.3703)
vff.add_angle = "Ga", "As", "Ga", ("tet", -4.099, 9.3703)
vff.add_angle = "In", "As", "In", ("tet", -5.753, 5.7599)
vff.add_angle = "As", "In", "As", ("tet", -5.753, 5.7599)
vff.add_angle = "Ga", "As", "In", (-0.35016, -4.926, 7.5651)


structure = Structure()
structure.set_cell = (10.0, 0.5, 0.5),\
                     (0.00, 0.0, 0.5),\
                     (0.00, 0.5, 0.0)
structure.add_atoms = ((0.00, 0.00, 0.00), "Ga"),\
                      ((0.25, 0.25, 0.25), "As"),\
                      ((1.00, 0.00, 0.00), "Ga"),\
                      ((1.25, 0.25, 0.25), "As"),\
                      ((2.00, 0.00, 0.00), "In"),\
                      ((2.25, 0.25, 0.25), "As"),\
                      ((3.00, 0.00, 0.00), "In"),\
                      ((3.25, 0.25, 0.25), "As"),\
                      ((4.00, 0.00, 0.00), "Ga"),\
                      ((4.25, 0.25, 0.25), "As"),\
                      ((5.00, 0.00, 0.00), "In"),\
                      ((5.25, 0.25, 0.25), "As"),\
                      ((6.00, 0.00, 0.00), "In"),\
                      ((6.25, 0.25, 0.25), "As"),\
                      ((7.00, 0.00, 0.00), "Ga"),\
                      ((7.25, 0.25, 0.25), "As"),\
                      ((8.00, 0.00, 0.00), "Ga"),\
                      ((8.25, 0.25, 0.25), "As"),\
                      ((9.00, 0.00, 0.00), "Ga"),\
                      ((9.25, 0.25, 0.25), "As"), 
print vff
print structure

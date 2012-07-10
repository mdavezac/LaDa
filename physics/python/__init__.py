""" Physical quantities. 


    Holds things related to physical quantities, eg atomic units using the *quantities* package.
"""
__docformat__ = "restructuredtext en"

import quantities as pq

__all__ = [ 'a0', 'bhor_radius',\
            'h', 'planck', 'h_bar', 'reduced_planck', 'electronic_mass',\
            'Ry', 'rydberg', 'Kb', 'boltzmann' ]


h = pq.UnitQuantity("planck", 4.1356673310e-15 * pq.eV * pq.s)
""" Planck's constant. """
planck = h
""" Planck's constant. """
h_bar = pq.UnitQuantity("h_bar", h / 2e0 / pq.pi, symbol='h_bar')
""" Reduced Planck's constant. """
reduced_plank = h_bar
""" Reduced Planck's constant. """

Ry = pq.UnitQuantity('Rydberg', 0.5 * pq.hartree, symbol='Ry')
""" Rydberg energy units. """
rydberg = Ry
""" Rydberg energy units. """

a0 = pq.UnitQuantity('bhor_radius', 0.529177249 * pq.angstrom, symbol='a0')
""" Bhor radius, unit of length of atomic units. """
bhor_radius = a0
""" Bhor radius, unit of length of atomic units. """

electronic_mass = pq.UnitQuantity("electronic_mass", h_bar**2 / (2e0 * Ry * a0**2)  )
""" Mass of the electron at rest.

    The value is obtained from a formula. It comes close enough and makes the
    Rydberg units consistent.
"""


Kb = pq.UnitQuantity("boltzmann", 8.617 * pq.eV / pq.K, symbol='Kb')
""" Boltzmann's constant. """
boltzmann = Kb
""" Boltzmann's constant. """

reduced_reciprocal_au = pq.UnitQuantity("reduced_reciprocal_au",
                                        2e0*pq.pi/a0, symbol='2pi/a0')

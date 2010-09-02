""" Physical quantities. 


    Holds things related to physical quantities, eg atomic units using the *quantities* package.
"""
__docformat__ = "restructuredtext en"

from _physics import Z, Symbol, Charge, Mass 
import quantities as pq

__all__ = [ 'Z', 'Symbol', 'Charge', 'Mass', 'a0', 'bhor_radius',\
            'h', 'planck', 'h_bar', 'reduced_planck', 'electronic_mass',\
            'Ry', 'rydberg', 'Kb', 'boltzmann' ]

a0 = pq.UnitQuantity('bhor_radius', 0.5219177 * pq.angstrom, symbol='a0')
""" Bhor radius, unit of length of atomic units. """
bhor_radius = a0
""" Bhor radius, unit of length of atomic units. """

h = pq.UnitQuantity("planck", 4.1357 * pq.eV * pq.s)
""" Planck's constant. """
planck = h
""" Planck's constant. """
h_bar = pq.UnitQuantity("h_bar", h / 2e0 / pq.pi, symbol='h_bar')
""" Reduced Planck's constant. """
reduced_plank = h_bar
""" Reduced Planck's constant. """

electronic_mass = pq.UnitQuantity("electronic_mass", 9.1095e-31 * pq.kg)
""" Mass of the electron at rest. """

Ry = pq.UnitQuantity('Rydberg', 13.60569172 * pq.eV, symbol='Ry')
""" Rydberg energy units. """
rydberg = Ry
""" Rydberg energy units. """

Kb = pq.UnitQuantity("boltzmann", 8.617 * pq.eV / pq.K, symbol='Kb')
""" Boltzmann's constant. """
boltzmann = Kb
""" Boltzmann's constant. """

reduced_reciprocal_au = pq.UnitQuantity("reduced_reciprocal_au",
                                        2e0*pq.pi/a0, symbol='* 2pi/a0')

#!/usr/bin/env python

from re import compile
from lada import jobs
from lada.opt import read_input
from lada.opt.changedir import Changedir
from lada.crystal import A2BX4, fill_structure, specie_list
from lada.vasp import Vasp, Specie, specie, files
from lada.vasp.methods import RelaxCellShape

# names we need to create input.
input_dict = { "Vasp": Vasp, "U": specie.U, "nlep": specie.nlep, "RelaxCellShape": RelaxCellShape }
# reads input file
input = read_input("input.py", input_dict)
# list of all lattices in A2BX4.
lattices = [ A2BX4.S1(), A2BX4.S1I(), A2BX4.S2(), A2BX4.S3(), A2BX4.S3I(),
             A2BX4.b1(), A2BX4.b10(), A2BX4.b10I(), A2BX4.b11(), A2BX4.b12(),
             A2BX4.b15(), A2BX4.b16(), A2BX4.b18(), A2BX4.b19(), A2BX4.b1I(),
             A2BX4.b2(), A2BX4.b20(), A2BX4.b21(), A2BX4.b2I(), A2BX4.b33(),
             A2BX4.b36(), A2BX4.b37(), A2BX4.b38(), A2BX4.b4(), A2BX4.b4I(), 
             A2BX4.b5(), A2BX4.b5I(), A2BX4.b8(), A2BX4.b9(), A2BX4.b9I(), 
             A2BX4.d1(), A2BX4.d1I(), A2BX4.d3(), A2BX4.d3I(), A2BX4.d9() ]

# regex
specie_regex = compile("([A-Z][a-z]?)2([A-Z][a-z]?)([A-Z][a-z]?)4")

# Job dictionary.
jobdict = jobs.JobDict()

# loop over materials.
for material in input.materials:

  # creates dictionary to replace A2BX4 with meaningfull species.
  match = specie_regex.match(material)
  assert match != None, RuntimeError("Incorrect material " + material + ".")
  # checks species are known to vasp functional
  for i in range(1, 4):
    assert match.group(i) in input.vasp.species,\
           RuntimeError("%s not in specie dictionary." % (match.group(i)))
  # actually creates dictionary.
  species_dict = {"A": match.group(1), "B": match.group(2), "X": match.group(3)}

  # loop over lattices. 
  for lattice in lattices:

    # creates a structure.
    structure = fill_structure(lattice.cell, lattice)
    # changes atomic types.
    for atom in structure.atoms:
      atom.type  = species_dict[atom.type]
    # set volume
    structure.scale = input.volume(structure)

    # job dictionary for this lattice.
    lat_jobdict = jobdict / material / lattice.name / "non-magnetic"

    # sets up job parameters.
    lat_jobdict.vasp = input.relaxer
    lat_jobdict.args = structure
    lat_jobdict.add_param = "ispin", 1

    # goes through magnetic stuff
    species = specie_list(structure)
    if len([0 for u in specie if u.magnetic]) == 0: continue

    # some needed lists
    moments = [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0] # high-spin, d elctrons.
    numbers = [atom.type for atom in structure.atoms]
    numbers = [ numbers.count(u) for u in species ]
    valence = [ input.vasp.species[s] for s in species ]

    # first ferro
    magom = ""
    for n, s, val  in zip(numbers, species, valence): 
      ps = input.vasp.species[s]
      magmom +="%i*%f   " % (n, moments[ps.valence] if ps.magnetic else 0)
    job = lat_jobdict / "ferro"
    job.vasp = input.relaxer
    job.args = structure
    job.add_param = "ispin", 2
    job.add_param = "magmom", magmom
    
    # then, antiferro, one for each cation.
    if len([0 for u in specie if u.magnetic]) == 2: 
      magom, sign = "", 1e0
      for n, s, val  in zip(numbers, species, valence): 
        ps = input.vasp.species[s]
        magmom +="%i*%f   " % (n, sign * moments[ps.valence] if ps.magnetic else 0)
        if ps.magnetic: sign *= -1e0
      job = lat_jobdict / "anti-ferro-0"
      job.vasp = input.relaxer
      job.args = structure
      job.add_param = "ispin", 2
      job.add_param = "magmom", magmom

   # Then random anti-ferro.
   for i in range(input.nbantiferro):
     magom = "", 1e0
     for n, s, val  in zip(numbers, species, valence): 
       ps = input.vasp.species[s]
       if ps.magnetic:
         magmom += 
       for i in range(n):
       magmom +="%i*%f   " % (n, moments[ps.valence] if ps.magnetic else 0)
     job = lat_jobdict / "anti-ferro-%i" % (i+1)
     job.vasp = input.relaxer
     job.args = structure
     job.add_param = "ispin", 2
     job.add_param = "magmom", magmom
                                         

with Changedir(input.outdir) as pwd: jobs.save(jobdict, "jobdict")

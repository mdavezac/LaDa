""" Creates job-dictionary of GA jobs. """
def gajobs(path, inputpath = "input.py"):
  from copy import deepcopy
  import IPython.ipapi
  from numpy import array, all

  from pylada.ga.escan.elemental.functional import Darwin as Functional
  from pylada.ga.escan.elemental.evaluator  import EffectiveMass as Evaluator
  from pylada.ga.escan.elemental            import Converter
  from pylada.ga.functional                 import minimize
  from pylada.jobs import JobFolder
  from pylada.escan import read_input

  # reads input file.
  input = read_input(inputpath)

  jobfolder = JobFolder()
  for type in ["ee", "hh"]:
    for direction in ["epi", "perp"]:
      for scale in input.scales:
        for range in input.ranges: 
          nmin, nmax = min(range), max(range)
          for trial in input.trials:
            escan = deepcopy(input.escan)
            escan.vff.lattice.scale = scale
            if direction == "perp": escan.direction = input.growth_direction
            elif all(array(input.growth_direction) - (0,0,1) == (0,0,0)):
              escan.direction = (1, 0, 0), (0, 1, 0)
            elif all(array(input.growth_direction) - (1,1,0) == (0,0,0)):
              escan.direction = (0, 0, 1), (-1, 1, 0)
            else: raise ValueError("Unknown growth direction.")
            if type == "ee": escan.type = "e"
            elif type == "hh" or type == "lh": escan.type = "h"
    
            converter = Converter(growth=input.growth_direction, lattice=escan.vff.lattice)
            
            kwargs = { "popsize": input.population_size, "rate": input.offspring_rate,
                       "max_gen": input.max_generations, "pools": input.pools,
                       "crossover_rate": input.crossover_rate, "swap_rate": input.swap_rate, 
                       "growth_rate": input.growth_rate, "nmin": nmin, "nmax": nmax,
                       "dosym": input.dosym, "rootworkdir": "$GSCRATCH", "comparison": minimize }
            kwargs.update(getattr(input, 'kwargs', {}))
            evaluator = Evaluator( converter = converter, 
                                   functional = escan, 
                                   **getattr(input, 'kwargs', {}) )
            functional = Functional(evaluator, **kwargs)
            functional.n = slice(2,4) if type == "lh" else slice(0,2)
    
            name = "{0[0]}{0[1]}{0[2]}/{1}/{2}/"\
                   .format(input.growth_direction, direction, type)
            name += "scale_{0:.2f}/{1}_{2}/trial_{3}".format(scale, nmin, nmax, trial)
            gajob = jobfolder / name
            gajob.functional = functional 

  ip = IPython.ipapi.get()
  ip.user_ns["current_jobfolder"] = jobfolder
  ip.magic("savejobs " + path)

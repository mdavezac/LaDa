def create_jobs(path, inputpath="input.py", **kwargs):
  """ Creates GA-jobs. 
  
      :Parameters:
        path 
          Path where the job-dictionary will be saved. Calculations will be
          performed in the parent directory of this file. Calculations will be
          performed in same directory as this file.
        inputpath
          Path to an input file. Defaults to input.py. 
        kwargs
          Any keyword/value pair to take precedence over anything in the input
          file.
  """
  from IPython.ipapi import get as get_ipy

  from lada.jobs import JobDict
  from lada.ga.escan.nanowires import read_input
  from lada.ga.escan.nanowires.converter import Converter
  from lada.ga.escan.nanowires.functional import Darwin

  input = read_input(inputpath)

  lattice = input.vff.lattice
  types = set([u for a in lattice.sites for u in a.type]) - set([input.passivant])

  jobdict = JobDict()
  for core_type in input.core_types:
    for nmin, nmax in input.ranges:
      for trial in xrange(input.nb_trials):
        # create conversion object from bitstring to nanowires and back.
        converter = Converter( input.vff.lattice, growth=input.growth_direction,
                               core_radius=input.core_radius, core_type=core_type,
                               types = list(types), thickness=input.thickness,
                               sep = input.separation )
        # create objective fucntion.
        evaluator = input.Evaluator( converter, input.escan,
                                     degeneracy = input.__dict__.get('degeneracy', 1e-3) )
        
        # create GA functional itself.
        functional = Darwin(evaluator, nmin=nmin, nmax=nmax, **input.__dict__)
        functional.popsize     = input.population_size
        functional.max_gen     = input.max_generations
        functional.rate        = input.offspring_rate
        functional.pools       = input.pools
        if input.rootworkdir != None: functional.rootworkdir = input.rootworkdir
    
        gajob = jobdict / "{0}_core/{2}_{3}/trial_{1}".format(core_type, trial, nmin, nmax)
        gajob.functional = functional



  # saves jobs.
  ip = get_ipy()
  ip.user_ns["current_jobdict"] = jobdict
  ip.magic("savejobs " + path)


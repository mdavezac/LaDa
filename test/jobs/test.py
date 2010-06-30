""" Tests job dictionary and semaphores. """ 
import random
from os import getcwd
from optparse import OptionParser
from boost.mpi import world, broadcast
from lada import jobs
from functional import Functional

# Reads program options.
parser = OptionParser()
parser.add_option( "--save", dest="saveme", help="Pickle job dictionary.",\
                   metavar="FILE", type="str")
parser.add_option( "--load", dest="loadme", help="Loads a pickled job dictionary.", \
                   metavar="FILE", type="str")
parser.add_option( "--pools", dest="pools", help="Number of mpi pools.", \
                   metavar="N", type="int", default=1)
parser.add_option( "--pbs", dest="pbs", help="Directory where to write pbs stuff.", \
                   metavar="DIRECTORY", type="str")
(options, args) = parser.parse_args()

# takes care of mpi pools
local_comm = world.split(world.rank % options.pools)

# creates a functional instance.
functional = Functional()

# Creates dictionary
if options.loadme == None or options.saveme != None:
  job_dictionary = jobs.JobDict()
  for caps in ["A", "B", "C", "D"]:
    capsjob = job_dictionary / caps

    if caps == "B": 
      capsjob.vasp = functional
      capsjob.args = "beast", 666
      capsjob["indict"] = "indicted"

    for numbers in ["1", "2", "3"]: 
      numbersjob = capsjob / numbers

      if numbers == "1" and caps == "A": continue
      numbersjob.vasp = functional
      numbersjob.args = caps, numbers, caps in ["B", "D"]
      numbersjob["indict"] = "indicted"
      numbersjob["something"] = "else"
# loads dictionary
else: job_dictionary = jobs.load(options.loadme, comm=world)

# saves dictionary
if options.saveme != None:
  jobs.save(jobdict=job_dictionary, path=options.saveme, overwrite=True, comm=world)
# writes pbs stuff.
if options.pbs != None and world.rank == 0:
  jobs.pbs_script( outdir="results", jobdict=job_dictionary, pools=options.pools,\
                   queue="regular", mppwidth=world.size, python_path=getcwd() ) 

# Computes all jobs.
if options.loadme == None and options.saveme == None and options.pbs == None:
  for job, outdir in job_dictionary.walk_through("results"):
    # launch jobs and stores result
    result = job.compute(outdir=outdir)
    # Root process of pool prints result.
    if local_comm.rank == 0: print result, "\n"
# Executes jobs using jobs.bleed
elif options.loadme != None: 
  for job, outdir in jobs.bleed(options.loadme, outdir=options.pbs, comm=local_comm):
    # decides on the same waittime for all processes.
    waittime = broadcast(local_comm, random.randint(0,2), 0)
    # launch jobs and stores result
    result = job.compute(outdir=outdir, wait=waittime)
    # Root process of pool prints result.
    if local_comm.rank == 0: print result, "\n"

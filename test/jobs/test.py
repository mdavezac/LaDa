""" Tests job dictionary and semaphores. """ 
import random
from optparse import OptionParser
from boost.mpi import world, broadcast
from lada import jobs

# Reads program options.
parser = OptionParser()
parser.add_option( "--save", dest="saveme", help="Pickle job dictionary.", metavar="FILE", type="str")
parser.add_option( "--load", dest="loadme", help="Loads a pickled job dictionary.", metavar="FILE", type="str")
parser.add_option( "--mpipools", dest="mpipools", help="Number of mpipools.", metavar="N", type="int", default=1)
(options, args) = parser.parse_args()

# takes care of mpi pools
if options.mpipools > world.size: options.mpipools = world.size
local_comm = world.split(world.rank % options.mpipools)

# Functional which will be called
class Functional(object):
  def __call__(self, *args, **kwargs):
    from time import sleep
    
    if "wait" in kwargs: sleep(kwargs["wait"])
    return "Rank: %i\nargs: %s\nkwargs: %s" % (world.rank, args, kwargs)

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
if options.saveme != None: jobs.save(path=options.saveme, overwrite=True, comm=world)

# Computes all jobs.
if options.loadme == None and options.saveme == None:
  for job, outdir in job_dictionary.walk_through("results"):
    # launch jobs and stores result
    result = job.compute(outdir=outdir)
    # Root process of pool prints result.
    if local_comm.rank == 0: print result, "\n"
# Executes jobs using jobs.bleed
elif options.loadme != None: 
  for job, outdir in jobs.bleed(options.loadme, outdir="results", comm=local_comm):
    # decides on the same waittime for all processes.
    waittime = broadcast(local_comm, random.randint(0,2), 0)
    # launch jobs and stores result
    result = job.compute(outdir=outdir, wait=waittime)
    # Root process of pool prints result.
    if local_comm.rank == 0: print result, "\n"

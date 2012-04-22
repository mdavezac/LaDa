class ExtractSingle(object): 
  def __init__(self, outdir):
    from os.path import exists, isfile, isdir, join
    from re import search
    from lada.misc import RelativePath
    super(ExtractSingle, self).__init__()

    outdir = RelativePath(outdir).path
    if isdir(outdir): outdir = join(outdir, 'stdout')
    self.directory = outdir
    self.success = False
    if not exists(outdir): return
    if not isfile(outdir): return

    try:
      with open(outdir, 'r') as file:
        line = file.next()
        regex = search( "pi to order (\d+) is approximately (\S+), "
                        "Error is (\S+) "\
                        "\s*-- slept (\d+) seconds at each iteration -- "
                        "\s*mpi world size is (\d+)", line )
        if regex is None: return 
        self.order = int(regex.group(1))
        self.pi    = float(regex.group(2))
        self.error = float(regex.group(3))
        self.sleep = int(regex.group(4))
        self.comm  = {'n': int(regex.group(5))}
  
        line = file.next()
        self.system = line[line.find(':')+1:].rstrip().lstrip()
        line = file.next()
        self.nodename = line[line.find(':')+1:].rstrip().lstrip()
        line = file.next()
        self.release = line[line.find(':')+1:].rstrip().lstrip()
        line = file.next()
        self.version = line[line.find(':')+1:].rstrip().lstrip()
        line = file.next()
        self.machine = line[line.find(':')+1:].rstrip().lstrip()
    except: self.success = False
    else: self.success = True

class ExtractMany(object):
  def __init__(self, outdir, order=None):
    from glob import iglob
    from os.path import isfile, join, basename
    from lada.misc import RelativePath

    super(ExtractMany, self).__init__()

    self.directory = RelativePath(outdir).path
    checkdir = order is None

    if checkdir:
      order = []
      for file in iglob(join(self.directory, 'stdout*')):
        o = int(basename(file)[6:])
        extract = ExtractSingle(file)
        if not extract.success: continue
        if extract.order == o: order.append(o)

    self.order   = []
    self.sleep   = []
    self.pi      = []
    self.error   = []
    self.system  = []
    self.version = []
    self.machine = []
    self.release = []
    self.comm    = []
    self.success = False
    error = False
    try: 
      for o in order:
        extract = ExtractSingle(join(self.directory, 'stdout')+str(o))
        if not extract.success: error = True; continue
        if extract.order != o: error = True; continue
        self.order.append(o)
        self.pi.append(extract.pi)
        self.error.append(extract.error)
        self.system.append(extract.system)
        self.version.append(extract.version)
        self.machine.append(extract.machine)
        self.release.append(extract.release)
        self.comm.append(extract.comm)
    except: self.success = False
    else: self.success = not error



from lada.functools import stateless
class Functional(object):
  def __init__(self, program, order=4, sleep=0):
    super(Functional, self).__init__()
    self.program = program
    self.order = order
    self.sleep = sleep

  def iter(self, outdir=None, sleep=None, overwrite=False, comm=None):
    from copy import deepcopy
    from os.path import join
    from lada.process.program import ProgramProcess
    from lada.misc import RelativePath
    self = deepcopy(self)
    outdir = RelativePath(outdir).path
    if sleep is not None: self.sleep = sleep
    order = self.order if hasattr(self.order, '__iter__') else [self.order]
    for o in order:
      stdout = join(outdir, 'stdout{0}'.format(o))
      stderr = join(outdir, 'stderr{0}'.format(o))
      if overwrite == False: 
        extract = ExtractSingle(stdout)
        if extract.success:
          yield extract
          continue
      yield ProgramProcess( self.program, cmdline=['--order', str(o), '--sleep', str(self.sleep)],
                            outdir=outdir, stdout=stdout, stderr=stderr, dompi=True, comm=comm)
  
  def __call__(self, outdir=None, sleep=None, overwrite=False, comm=None):
    from lada.misc import RelativePath
    outdir = RelativePath(outdir).path
    for program in self.iter(outdir, sleep, overwrite):
      if getattr(program, 'success', False) == False:
        program.wait()
    return self.Extract(outdir)

  def Extract(self, outdir):
    order = self.order if hasattr(self.order, '__iter__') else [self.order]
    return ExtractMany(outdir, order=order) 

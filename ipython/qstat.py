__docformat__ = "restructuredtext en"
__all__ = ['qstat']
from lada import lada_with_slurm

if lada_with_slurm:
  def qstat(self, arg):
    """ squeue --user=`whoami` -o "%7i %.3C %3t  --   %50j" """
    from subprocess import Popen, PIPE
    from IPython.genutils import SList
    from getpass import getuser

    # finds user name.
    whoami = getuser()
    squeue = Popen(["squeue", "--user=" + whoami, "-o", "\"%7i %.3C %3t    %j\""], stdout=PIPE)
    result = squeue.stdout.read().rstrip().split('\n')
    result = SList([u[1:-1] for u in result[1:]])
    return result.grep(str(arg[1:-1]))
else:
  def qstat(self, arg):
    """ Prints jobs of current user. """
    from subprocess import Popen, PIPE
    from getpass import getuser
    from IPython.genutils import SList
    from re import compile
    # get user jobs ids
    ids = Popen('qstat -u{0}'.format(getuser()).split(), stdout=PIPE).communicate()[0].split('\n')
    ids = SList(ids).grep(getuser()).fields(0)
  
    result = SList()
    name_re = compile("Job_Name\s*=\s*(.+)")
    state_re = compile("job_state\s*=\s*(\S+)")
    mppwidth_re = compile("Resource_List.mppwidth\s*=\s*(\d+)")
    for id in ids:
      full = Popen('qstat -f {0}'.format(id).split(), stdout=PIPE).communicate()[0]
      name = name_re.search(full).group(1)
      state = state_re.search(full).group(1)
      mppwidth = mppwidth_re.search(full).group(1)
      result.append("{0:>10} {1:>4} {2:>3} -- {3}".format(id, mppwidth, state, name))
    return result


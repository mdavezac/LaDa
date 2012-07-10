#! /work/y07/y07/nag/python-shared-xe6/2.6.5/bin/python
#PBS -e global_comm_err
#PBS -o global_comm_out
#PBS -N global_comm
#PBS -l mppwidth=64
#PBS -l walltime=00:02:00
#PBS -A e05-qmdev-nic
#PBS -V 

def test(nbprocs, ppn):
  import lada
  from lada.process.mpi import create_global_comm

  lada.default_comm['ppn'] = ppn
  lada.default_comm['n'] = nbprocs
  print 'EXPECTED N={0}, PPN={1}'.format(nbprocs, ppn)
  create_global_comm(nbprocs)
  print 'FOUND'
  for u in lada.default_comm.iteritems(): print u[0], u[1]
  print 'MACHINES'
  for u in lada.default_comm.machines.iteritems(): print u[0], u[1]

if __name__ == '__main__': 
  test(nbprocs=64, ppn=32)

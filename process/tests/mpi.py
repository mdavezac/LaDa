def test():
  """ Test MPI Communicator. """
  from lada.process.mpi import Communicator, ProcessNumberError

  root = Communicator(n=32)
  for i in xrange(4): root.machines["node0{0}".format(i)] = 8

  newcomm = root.lend(5)
  assert newcomm['n'] == 5
  assert newcomm.parent() is root
  assert len(newcomm.machines) == 1
  assert root.machines[newcomm.machines.keys()[0]] == 3
  assert root['n'] == 27
  newcomm.cleanup()
  assert newcomm['n'] == 0
  assert len(newcomm.machines) == 0
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())

  newcomm = root.lend(8)
  assert newcomm['n'] == 8
  assert sum(newcomm.machines.itervalues()) == newcomm['n']
  assert newcomm.parent() is root
  assert len(newcomm.machines) == 1
  key = newcomm.machines.keys()[0]
  assert key not in root.machines
  assert newcomm.machines[key] == 8
  assert root['n'] == 24
  newcomm.cleanup()
  assert newcomm['n'] == 0
  assert len(newcomm.machines) == 0
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())

  newcomm = root.lend(12)
  assert newcomm['n'] == 12
  assert sum(newcomm.machines.itervalues()) == newcomm['n']
  assert newcomm.parent() is root
  assert len(newcomm.machines) == 2
  key0, key1 = newcomm.machines.keys()
  if newcomm.machines[key0] != 8: key0, key1 = key1, key0
  assert newcomm.machines[key0] == 8
  assert newcomm.machines[key1] == 4
  assert key0 not in root.machines
  assert root.machines[key1] == 4
  assert root['n'] == 20
  newcomm.cleanup()
  assert newcomm['n'] == 0
  assert len(newcomm.machines) == 0
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())

  comms = root.split(4)
  assert root['n'] == 0
  assert len(root.machines) == 0
  machines = []
  for comm in comms: 
    assert comm['n'] == 8
    assert sum(comm.machines.itervalues()) == comm['n']
    assert len(comm.machines) == 1
    assert comm.machines.keys()[0] not in machines
    machines.append(comm.machines.keys()[0])
  for comm in comms: comm.cleanup()
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())
  
  comms = root.split(5)
  assert root['n'] == 0
  assert len(root.machines) == 0
  machines = {}
  for comm in comms: 
    assert comm['n'] in [6, 7]
    assert sum(comm.machines.itervalues()) == comm['n']
    for key, value in comm.machines.iteritems():
      if key not in machines: machines[key] = value
      else: machines[key] += value
  assert sum(machines.itervalues()) == 32
  assert all(u == 8 for u in machines.itervalues())
  for comm in comms: comm.cleanup()
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())

  comms = root.split(3)
  assert root['n'] == 0
  assert len(root.machines) == 0
  machines = {}
  for comm in comms: 
    assert comm.parent() is root
    assert comm['n'] in [10, 11]
    assert sum(comm.machines.itervalues()) == comm['n']
    for key, value in comm.machines.iteritems():
      if key not in machines: machines[key] = value
      else: machines[key] += value
  assert sum(machines.itervalues()) == 32
  assert all(u == 8 for u in machines.itervalues())

  machines = comms[0].machines.copy()
  for key, value in comms[1].machines.iteritems():
    if key in machines: machines[key] += value
    else: machines[key] = value
  comm = comms.pop(0)
  comms[0].acquire(comm)
  assert comm.parent is None
  assert comm['n'] == 0
  assert len(comm.machines) == 0
  assert comms[0].parent() is root
  assert comms[0]['n'] == sum(machines.itervalues())
  assert comms[0]['n'] == sum(comms[0].machines.itervalues())
  for key in machines:
    assert machines[key] == comms[0].machines[key]
  for key in comms[0].machines:
    assert machines[key] == comms[0].machines[key]
  for comm in comms: comm.cleanup()
  assert root['n'] == 32
  assert all(u == 8 for u in root.machines.itervalues())

  try: comm.lend(33)
  except ProcessNumberError: pass
  else: raise Exception()
  try: comm.split(33)
  except ProcessNumberError: pass
  else: raise Exception()

if __name__ == "__main__":
  from sys import argv, path
  from os.path import abspath
  if len(argv) > 1: path.extend(argv[1:])
  test()

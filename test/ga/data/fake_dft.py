import os
from boost import mpi
import sys
import os


# load correct path on lester.
if os.uname()[1] == "head": 
  sys.path.append("/uhome/mdavezac/usr/lib/python2.6/site-packages")

from lada import ce, crystal


functional = ce.Cubic()
# try:
functional.load( "ce.xml" )
functional.set_mpi(mpi.world)
# except: exit



for filename in os.listdir("."):
  if filename == "ce.xml": continue
  if filename == ".fake_dft.py.swp": continue
  if filename == "fake_dft.py": continue
  if filename == "LDAs.dat": continue
  try:
   structure = crystal.read_structure( filename )
   print   "%45s %8.4f " % ( structure.name, functional.chemical(structure) )
  except: continue
      

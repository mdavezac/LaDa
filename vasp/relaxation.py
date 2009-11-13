#
#  Version: $Id$
#

def relaxation( _structure, _vasp, tolerance=1e-3 ):
  from math import fabs as abs
  import os

  indir = str(_vasp.indir)
  assert os.path.exists(indir) and os.path.isdir(indir),\
        "%s does not exist or is not a directory.\n" % (indir)
  index = 1
  _vasp.indir  = "%s/step_%i" % ( indir, 0 )
  _vasp.relaxation = "global volume ionic cellshape"
  if not os.path.exists(_vasp.indir): 
    os.makedirs(_vasp.indir)
  else: assert os.path.isdir(_vasp.indir), "%s exists and is not a directory.\n" % (_vasp.indir)
  
  oldenergy = float( _structure.energy )
  yield _structure
  structure = _vasp.final()

  tol = tolerance * float( len(_structure.atoms) )
  while( abs(oldenergy - structure.energy) > tol ):
    _vasp.restart = str( _vasp.indir )
    _vasp.indir  = "%s/step_%i" % ( indir, index )
    if not os.path.exists(_vasp.indir): 
      os.makedirs(_vasp.indir)
    yield structure
    index += 1
    structure = _vasp.final()
    oldenergy = float(structure.energy)

  _vasp.indir  = "%s/final_static" % ( indir )
  if not os.path.exists(_vasp.indir): 
    os.makedirs(_vasp.indir)
  _vasp.relaxation = "static"
  yield structure
  structure = _vasp.final()


def main():
  from lada import atat, crystal
  import os.path
  import shutil
  import subprocess
  from  vasp import Vasp, Specie

  structure = crystal.sStructure()
  structure.cell = atat.rMatrix3d( [[1,0,0],[0,1,0],[0,0,1]] )
  atom = crystal.StrAtom()
  atom.type = "Rb"
  atom.pos = atat.rVector3d([0,0,0])
  structure.atoms.append( crystal.StrAtom(atom) )
  atom.type = "K"
  atom.pos = atat.rVector3d([0.5,0.5,0.5])
  structure.atoms.append( crystal.StrAtom(atom) )
  structure.name = "KRb_nosym"
  structure.scale = 6

  K = Specie( "K", "~/AtomicPotentials/pseudos/K_s" )
  Rb = Specie( "Rb", "~/AtomicPotentials/pseudos/Rb_s" )

  vasp = Vasp()
  vasp.isym = "off"
  vasp.indir = str(structure.name)
  if os.path.exists( vasp.indir ):  shutil.rmtree( vasp.indir )
  script = vasp.prepare( structure, wpath="KRb" )
  file = open( os.path.join( vasp.indir, 'script' ), 'w' )
  print >>file, script
  file.close()
  subprocess.call( ["bash", os.path.join(vasp.indir, "script") ] )

  for structure in relaxation( structure, vasp ):
    print structure.name, structure.energy, "\n", structure
    if os.path.exists( vasp.indir ): continue
    script = vasp.prepare( structure, wpath=vasp.indir )
    file = open( os.path.join( vasp.indir, 'script' ), 'w' )
    print >>file, script
    print script
    file.close()
    subprocess.call( ["bash", os.path.join(vasp.indir, "script") ] )
  

if __name__ == "__main__":
  main()

    

    
 

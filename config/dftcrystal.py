CRYSTAL_geom_blocks = set(['CRYSTAL', 'SLAB', 'POLYMER', 'HELIX', 'MOLECULE'])
""" List of starting blocks in CRYSTAL input.

    CRYSTAL input does not differentiate between its block and keyword inputs.
    As such, to parse its input. the name of the blocks must be known
    explicitely. 
    This particular set is used to figure out where the input starts in an
    output file.
"""
CRYSTAL_input_blocks = set([ 'MARGINS', 'BIDIERD', 'CPHF', 'ELASTCON', 'EOS',
                             'SYMMWF', 'LOCALWF', 'ANISOTRO', 'ECH3', 'EDFT',
                             'EIGSHROT', 'OPTGEOM', 'FIXINDEX', 'GRID3D',
                             'MAPNET', 'POT3', 'DFT', 'PRINTOUT',
                             'REFLECTANCE', 'ROTCRY' ])            \
                       | CRYSTAL_geom_blocks
""" List of blocks in CRYSTAL input.

    CRYSTAL input does not differentiate between its block and keyword inputs.
    As such, to parse its input. the name of the blocks must be known
    explicitely.
"""

CRYSTAL_filenames = { 'crystal.out':   '{0}.out',      # output file
                      'crystal.err':   '{0}.err',      # error file
                      'crystal.d12':   '{0}.d12',      # input file
                      'fort.9':        '{0}.f9',       # binary wave-functions
                      'fort.98':       '{0}.f98',      # formatted wave-functions
                      'GAUSSIAN.DAT':  '{0}.gjf',      # Gaussian 94/98 input
                      'MOLDRAW.DAT':   '{0}.mol',      # MOLDRAW input
                      'fort.33':       '{0}.xyz',      # xyz/Xmol input
                      'fort.34':       '{0}.gui',      # DLV input
                      'FINDSYM.DAT':   '{0}.findsym',  # Findsym input
                      'OPTHESS.DAT':   '{0}.opthess',  # formatted hessian
                      'OPTINFO.DAT':   '{0}.optinfo',  # restart info
                      'PPAN.DAT':      '{0}.ppan'      # Muliken population analysis
                    }
""" Filnames of crystal programs. 

    Map from fortran output filename to desired name. That desired format
    should always be the same to make it easy to find directories with CRYSTAL
    outputs.
""" 
CRYSTAL_propnames = { 'fort.25':      '{0}.f25',      # bands, maps, doss data
                      'GRED.DAT':     '{0}.GRED',     # direct lattice
                      'KIBZ.DAT':     '{0}.KIBZ',     # reciprocal lattice, IBZ
                      'KRED.DAT':     '{0}.KRED',     # reciprocal lattice, full BZ
                      'LINEA.DAT':    '{0}.LINEA',    # EMD line
                      'PROF.DAT':     '{0}.PROF',     # EMD in a plane
                      'DIEL.DAT':     '{0}.DIEL',     # dielectric constant
                      'POTC.DAT':     '{0}.POTC',     # exact electrostatic potential
                      'fort.31':      '{0}.prop3d',   # charge/spin density/potential
                      'fort.8':       '{0}.localwf',  # wannier function
                      'freqinfo.DAT': '{0}.freqinfo'  # info for freq restart
                    }
crystal_program ='crystal'
""" Path to crystal executable. """

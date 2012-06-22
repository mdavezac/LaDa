CRYSTAL_geom_blocks = set(['CRYSTAL', 'SLAB', 'POLYMER', 'HELIX', 'MOLECULE'])
""" List of starting blocks in CRYSTAL input.

    CRYSTAL input does not differentiate between its block and keyword inputs.
    As such, to parse its input. the name of the blocks must be known
    explicitely. 
    This particular set is used to figure out where the input starts in an
    output file.
"""
CRYSTAL_input_blocks = set([ 'MARGINS', 'CORRELAT', 'EXCHANGE', 'BIDIERD',
                             'CPHF', 'ELASTCON', 'EOS', 'SYMMWF', 'LOCALWF',
                             'ANISOTRO', 'ECH3', 'EDFT', 'EIGSHROT', 'OPTGEOM',
                             'FIXINDEX', 'GRID3D', 'MAPNET', 'POT3',
                             'PRINTOUT', 'REFLECTANCE', 'ROTCRY' ])            \
                       | CRYSTAL_geom_blocks
""" List of blocks in CRYSTAL input.

    CRYSTAL input does not differentiate between its block and keyword inputs.
    As such, to parse its input. the name of the blocks must be known
    explicitely.
"""

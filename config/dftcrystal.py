CRYSTAL_input_blocks = set([ 'CRYSTAL', 'SLAB', 'POLYMER', 'HELIX', 'MOLECULE',
                             'MARGINS', 'CORRELAT', 'EXCHANGE', 'BIDIERD',
                             'CPHF', 'ELASTCON', 'EOS', 'SYMMWF', 'LOCALWF',
                             'ANISOTRO', 'ECH3', 'EDFT', 'EIGSHROT',
                             'FIXINDEX', 'GRID3D', 'MAPNET', 'POT3',
                             'PRINTOUT', 'REFLECTANCE', 'ROTCRY' ])
""" List of blocks in CRYSTAL input.

    CRYSTAL input does not differentiate between its block and keyword inputs.
    As such, to parse its input. the name of the blocks must be known
    explicitely.
"""

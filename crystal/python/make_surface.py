def make_surface(structure, miller, nlayers=5, vacuum=15, acc=5):
    """ Returns a slab from the 3D structure 
        
        Takes a structure and makes a slab defined by the miller indices 
        with nlayers number of layers and vacuum defining the size 
        of the vacuum thickness. Variable acc determines the number of 
        loops used to get the direct lattice vectors perpendicular 
        and parallel to miller. For high index surfaces use larger acc value 

        .. warning: #. cell is always set such that miller is alogn z-axes
                    #. nlayers and vacuum are always along z-axes.

        :param structure: 
            Bulk :py:class:`structure <lada.crystal.Structure>` from which to
            create a surface.
        :param miller: 3x1 float64 array
            Miller indices defining the slab    
        :param nlayers: integer
            Number of layers in the slab
        :param vacuum: real
            Vacuum thicness in angstroms
        :param acc: integer
            number of loops for finding the cell vectors of the slab structure
    """
    from lada.crystal import fill_structure, Structure
    from numpy import arange, sqrt, array, transpose, dot, cross
    from numpy.linalg import inv

    direct_cell=transpose(structure.cell)

    orthogonal = []  # lattice vestors orthogonal to miller

    for n1 in arange(-acc,acc+1):
        for n2 in arange(-acc,acc+1):
            for n3 in arange(-acc,acc+1):
                
                pom = array([n1,n2,n3])
                if dot(pom,miller) == 0 and dot(pom,pom) != 0:
                    orthogonal.append(array([n1,n2,n3]))

    # chose the shortest parallel and set it to be a3 lattice vector
    norm_orthogonal=[sqrt(dot(dot(x,direct_cell),dot(x,direct_cell))) for x in orthogonal]
    a1 = orthogonal[norm_orthogonal.index(min(norm_orthogonal))] 

    # chose the shortest orthogonal to miller and not colinear with a1 and set it as a2
    in_plane=[]

    for x in orthogonal:
        if dot(x,x)>0.:
            v=cross(dot(x,direct_cell),dot(a1,direct_cell))
            v=sqrt(dot(v,v))
            if v>0.:
                in_plane.append(x)

    norm_in_plane = [sqrt(dot(dot(x,direct_cell),dot(x,direct_cell))) for x in in_plane]
    a2 = in_plane[norm_in_plane.index(min(norm_in_plane))]

    a1 = dot(a1,direct_cell)
    a2 = dot(a2,direct_cell)

    # new cartesian axes z-along miller, x-along a1, and y-to define the right-hand orientation
    e1 = a1/sqrt(dot(a1,a1))
    e2 = a2 - dot(e1,a2)*e1
    e2 = e2/sqrt(dot(e2,e2))
    e3 = cross(e1,e2)

    # find vectors parallel to miller and set the shortest to be a3
    parallel = []

    for n1 in arange(-acc,acc+1):
        for n2 in arange(-acc,acc+1):
            for n3 in arange(-acc,acc+1):
                pom = dot(array([n1,n2,n3]),direct_cell)
                if sqrt(dot(pom,pom))-dot(e3,pom)<1e-8 and sqrt(dot(pom,pom))>0.:
                    parallel.append(pom)

    norm_parallel = [sqrt(dot(x,x)) for x in parallel]
    a3 = parallel[norm_parallel.index(min(norm_parallel))]

    # making a structure in the new unit cell - defined by the a1,a2,a3
    new_direct_cell = array([a1,a2,a3])
    structure = fill_structure(transpose(new_direct_cell),structure.to_lattice())

    # transformation matrix to new coordinates x' = dot(m,x)
    m = array([e1,e2,e3])

    # seting output structure
    out_structure = Structure()
    out_structure.name = structure.name + 'SLAB'
    out_structure.scale = structure.scale
    out_structure.set_cell = transpose(dot(new_direct_cell,transpose(m)))

    for atom in structure.atoms:
        out_structure.add_atom = [dot(m,atom.pos),atom.type]

    # repaeting to get nlayers and vacuum
    repeat_cell   = dot(out_structure.cell,array([[1.,0.,0.],[0.,1.,0.],[0.,0.,nlayers]]))
    out_structure = fill_structure(repeat_cell,out_structure.to_lattice())

    # checking whether there are atoms close to the cell faces and putting them back to zero
    for i in range(len(out_structure.atoms)):
        scaled_pos = dot(out_structure.atoms[i].pos,inv(transpose(out_structure.cell)))
        for j in range(3):
            if abs(scaled_pos[j]-1.)<1e-5:
                scaled_pos[j]=0.
        out_structure.atoms[i].pos = dot(scaled_pos,transpose(out_structure.cell))

    # adding vaccum to the cell
    out_structure.set_cell = out_structure.cell + array([[0.,0.,0.],[0.,0.,0.],[0.,0.,float(vacuum)/float(out_structure.scale)]])

    # translating atoms so that center of the slab and the center of the cell along z-axes coincide
    max_z = max([x.pos[2] for x in out_structure.atoms])
    min_z = min([x.pos[2] for x in out_structure.atoms])
    center_atoms = 0.5*(max_z+min_z)
    center_cell  = 0.5*out_structure.cell[2][2]

    for i in range(len(out_structure.atoms)):
        out_structure.atoms[i].pos = out_structure.atoms[i].pos + array([0.,0.,center_cell-center_atoms])

    # exporting the final structure
    return out_structure


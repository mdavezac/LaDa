//
//  Version: $Id$
//

#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/structure.h>
#include <lamarck/atom.h>

#include "layered.h"

namespace Layered
{


  bool Physics :: Load( const TiXmlElement &_node )
  {
    if ( not lattice.Load( _node ) )
    {
      std::cerr << " Could not load lattice type from input!! " << std::endl; 
      return false;
    }

    Ising_CE::Structure::lattice = &lattice;
    if (     lattice
         and parent->Attribute("direction") 
         and parent->Attribute("multiplicity") )
    {
      if ( not Load_Structure( *_parent ) ) 
      {
        std::cerr << "Found attributes for constructing a layered structure...\n" 
                  << "But something went wrong.\n"
                  << "Continuing with standard load.\n"; 
      }
      if ( not structure.Load( *_parent ) ) return false;
    }

    if ( not consistency_check() )  return false;

    return true;
  }

  bool Physics :: Load_Structure( const TiXmlElement& &_node )
  {
    rVector3d cdir;
    iVector3d direction;
    types :: t_unsigned multiplicity; 
    rVector3d cell = &structure.cell;
    
    // First, Load Attributes 
    std::istringstream sstr = _node.Attribute("direction");
    sstr >> direction[0]; if ( sstr.fail() ) goto failure;
    sstr >> direction[1]; if ( sstr.fail() ) goto failure;
    sstr >> direction[2]; if ( sstr.fail() ) goto failure;

    if ( atat::norm2( direction ) == 0 ) goto failure;
    
    unsigned u;
    _node.Attribute( "multiplicity", &u );
    if ( not u ) goto failure;
    multiplicity = (types::t_unsigned) u;


    // Then constructs unit cell
    cell = lattice->cell;
    cell(0,0) =   (types::t_real) ( d(0) * multiplicity ) * lattice->cell(0,0)
                + (types::t_real) ( d(1) * multiplicity ) * lattice->cell(0,1)
                + (types::t_real) ( d(2) * multiplicity ) * lattice->cell(0,2);
    cell(1,0) =   (types::t_real) ( d(0) * multiplicity ) * lattice->cell(1,0)
                + (types::t_real) ( d(1) * multiplicity ) * lattice->cell(1,1)
                + (types::t_real) ( d(2) * multiplicity ) * lattice->cell(1,2);
    cell(2,0) =   (types::t_real) ( d(0) * multiplicity ) * lattice->cell(2,0)
                + (types::t_real) ( d(1) * multiplicity ) * lattice->cell(2,1)
                + (types::t_real) ( d(2) * multiplicity ) * lattice->cell(2,2);

    // Checks that cell is not singular
    if ( std::abs( cell.det() ) < types::tolerance )
      cell.set_column(1, lattice->cell.get_column( 0 ) );
    if ( std::abs( cell.det() ) < types::tolerance )
    {
      std::cerr << "Could not construct unit-cell\n" << cell << std::endl;
      return false;
    }

    // Makes sure the triad is direct
    if ( std::abs( cell.det() ) < types::tolerance )
    {
      rVector3d d = cell.get_column(2);
      cell.set_column(2, cell.get_column(1) );
      cell.set_column(1, d);
    }

    // Now adds atoms to structure.
    structure.atoms.clear();

    // Atoms of site 0 can be added in the exact same way that kvectors are
    // found.  Indeed, lets fool that particular routine
    rMatrix3d holdcell = cell;
    cell =  !( ~(cell) );
    structure.find_k_vectors();
    // And undo the cell 
    cell = holdcell;
    
    // We now sort the real space atoms according to layer
    std::sort( structure.k_vecs.begin(), structure.k_vecs.end(),
               Depth( cell.get_column(0) ) );


    // Finally, we copy the kvector positions as atoms, and the related sites
    Ising_CE::Structure::t_kAtoms::const_iterator i_kvec = structure.k_vecs.begin();
    Ising_CE::Structure::t_kAtoms::const_iterator i_kvec_end = structure.k_vecs.end();
    structure.atoms.clear(); structure.atoms.reserve( structure.k_vecs.size() );
    bool only_one_site = lattice->sites.size() == 1;
    rVector3d origin = lattice->sites.front().pos;
    for(; i_kvec != i_kvec_end; ++i_kvec )
    {
      Ising_CE::Structure::t_Atom atom;
      atom.site = 0; atom.pos = i_kvec->pos;
      atom.type = Ising_CE::Structure::t_Atom::t_Type(0);
      atom.freeze = lattice->sites.front().freeze;
      structure.atoms.push_back(atom);

      if( only_one_site ) continue;

      Ising_CE::Lattice::t_Sites::const_iterator i_site = lattice->sites.begin();
      Ising_CE::Lattice::t_Sites::const_iterator i_site_end = lattice->sites.end();
      for(++i_site; i_site != i_site_end; ++i_site )
      {
        Ising_CE::Structure::t_Atom atom2(atom);
        atom2.pos += ( i_site->pos - origin );
        atom2.freeze = i_site->freeze;
      }
    }
   
    structure.k_vecs.clear();
    // More of the same... but for kvecs
    structure.find_k_vecs();
  }

  inline bool Depth::operator()( const rVector3d &_1, const rVector3d &_2 )
  {
    types::t_real a =   _1(0) * depth(0) + _1(1) * depth(1) + _1(2) * depth(2);
    types::t_real b =   _2(0) * depth(0) + _2(1) * depth(1) + _2(2) * depth(2);
    return t_QuantityTraits<double> :: less( a, b );
  }

} // namespace Layered


#endif

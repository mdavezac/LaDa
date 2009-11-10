//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <iomanip>
#include <algorithm>

#include <boost/lexical_cast.hpp>

#include <crystal/which_site.h>
#include <opt/types.h>
#include <opt/fuzzy.h>

#include "mlcluster.h"


namespace LaDa
{
  namespace CE 
  {
    size_t size_t_cast( types::t_int i ) 
    { 
      LADA_ASSERT( i >= 0, "Symmetry operation is not of lattice.\n" );
      return size_t(std::abs(i));
    }

    void MLCluster :: apply_symmetry( Crystal::SymmetryOperator const &_op )
    {
      LADA_ASSERT(Crystal::Structure::lattice != NULL, "No lattice is set.\n");
      Crystal::Lattice const &lattice = *Crystal::Structure::lattice;
      atat::rMatrix3d const inv_cell( !lattice.cell );

      // finds new origin.
      atat::rVector3d const transformed_pos(_op(origin.pos)); 
      origin.site = size_t_cast(which_site(transformed_pos, inv_cell, lattice.sites));
      origin.pos = lattice.sites[origin.site].pos;

      if ( size() == 0 ) return;

      t_Spins :: iterator i_spin = begin();
      t_Spins :: iterator const i_last = end();
      for(; i_spin != i_last; ++i_spin)
      {
        i_spin->pos = _op(i_spin->pos);
        i_spin->site = size_t_cast(which_site(i_spin->pos+origin.pos, inv_cell, lattice.sites));
      }
    }


    bool MLCluster :: load( const TiXmlElement &_node )
    {

      clear();
      LADA_ASSERT( _node.Attribute("site"), "Missing site attribute.\n" );
      origin.site = boost::lexical_cast<size_t>(_node.Attribute("site"));
      LADA_ASSERT( origin.site < Crystal::Structure::lattice->sites.size(), 
                   "Site index out-of-range.\n" );
      origin.pos = Crystal::Structure::lattice->sites[origin.site].pos;

      const TiXmlElement *child = _node.FirstChildElement( "spin" );
      for ( ; child; child=child->NextSiblingElement( "spin" ) )
      {
        LADA_ASSERT( child->Attribute("x"), "Missing x attribute.\n" );
        LADA_ASSERT( child->Attribute("y"), "Missing y attribute.\n" );
        LADA_ASSERT( child->Attribute("z"), "Missing z attribute.\n" );
        LADA_ASSERT( child->Attribute("site"), "Could not find site attribute.\n" );
        Spin const spin = {
                            boost::lexical_cast<size_t>( child->Attribute("site") ),
                            atat::rVector3d
                            (
                              boost::lexical_cast<types::t_real>( child->Attribute("x") ),
                              boost::lexical_cast<types::t_real>( child->Attribute("y") ),
                              boost::lexical_cast<types::t_real>( child->Attribute("z") ) 
                            )
                          };
        push_back( spin );
      }

      return true;
    }

    bool operator==(MLCluster::Spin const &_a, MLCluster::Spin const &_b) 
    {
      if( _a.site != _b.site ) return false;
      if( not Fuzzy::is_zero(_a.pos(0)-_b.pos(0)) ) return false;
      if( not Fuzzy::is_zero(_a.pos(1)-_b.pos(1)) ) return false;
      return  Fuzzy::is_zero(_a.pos(2)-_b.pos(2));
    }
    
    bool MLCluster::operator==( MLCluster const & _c ) const
    {
      if( origin != _c.origin ) return false;
      LADA_ASSERT( Fuzzy::is_zero( atat::norm2(origin.pos-_c.origin.pos) ), 
                   "Inconsistent origin positions.\n" )
      if( size() != _c.size() ) return false;

      t_Spins :: const_iterator i_spin = _c.begin();
      t_Spins :: const_iterator const i_spin_end = _c.end();
      for(; i_spin != i_spin_end; ++i_spin)
        if( end() == std::find(begin(), end(), *i_spin) )
          return false;
      return true;
    }

    types::t_real MLCluster::operator()( Crystal::Structure const &_str, 
                                         std::vector< std::vector<size_t> > const &_map,
                                         Crystal::t_SmithTransform const &_transform ) const
    {
      types::t_real result(0);
      Crystal::Structure::t_Atoms::const_iterator i_atom = _str.atoms.begin();
      Crystal::Structure::t_Atoms::const_iterator const i_atom_end = _str.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom)
        result += operator()(*i_atom, _str, _map, _transform);
      return result;
    }

    // Finds spin for a given position.
    types::t_real find_spin( atat::rVector3d const &_vec, 
                             size_t _site, Crystal::Structure const &_str, 
                             std::vector< std::vector<size_t> > const &_map,
                             Crystal::t_SmithTransform const &_transform )
    {
      LADA_ASSERT
      (
        _site == Crystal::which_site(_vec, !_str.lattice->cell, _str.lattice->sites),
        "Could not find corresponding spin.\n"
      );
      size_t const smith_index = Crystal::get_linear_smith_index(_transform, _vec);
      LADA_ASSERT
      (
        smith_index < _map[_site].size(),
        "Atomic map and smith index are not coherent.\n"
      );
      LADA_ASSERT
      (
        _map[_site][smith_index] < _str.atoms.size(),
        "Atomic map, smith index, and structure are not coherent.\n"
      );
      return _str.atoms[ _map[_site][smith_index] ].type;
    }

    types::t_real MLCluster::operator()( Crystal::Structure::t_Atom const &_atom,
                                         Crystal::Structure const &_str, 
                                         std::vector< std::vector<size_t> > const &_map,
                                         Crystal::t_SmithTransform const &_transform ) const
    {
      LADA_ASSERT(_str.lattice, "Lattice not set.\n");
      LADA_ASSERT(_map.size() == _str.lattice->sites.size(),
                  "Atomic map and lattice are not coherent.\n");
      Crystal::Lattice const &lattice = *_str.lattice;
      types::t_real result(0);
      if( _atom.site == origin.site ) // atom is  at cluster origin.
      {
        result = _atom.type;
        t_Spins::const_iterator i_spin = begin();
        t_Spins::const_iterator const i_spin_end = end();
        for(; i_spin != i_spin_end; ++i_spin)
        {
          LADA_ASSERT(i_spin->site < lattice.sites.size(), "Index out-of-range.\n");
          atat::rVector3d const vec(_atom.pos+i_spin->pos-lattice.sites[i_spin->site].pos);
          result *= find_spin(vec, i_spin->site, _str, _map, _transform);
        }
      } 

      // atom may be at cluster branch.
      t_Spins::const_iterator i_orig = begin();
      t_Spins::const_iterator const i_orig_end = end();
      for(; i_orig != i_orig_end; ++i_orig)
      {
        LADA_ASSERT(i_orig->site < lattice.sites.size(), "Index out-of-range.\n")
        if( _atom.site != i_orig->site ) continue;

        atat::rVector3d const center(_atom.pos-i_orig->pos);
        types::t_real intermediate
          = find_spin(center-lattice.sites[origin.site].pos, origin.site, _str, _map, _transform);
        t_Spins::const_iterator i_spin = begin();
        t_Spins::const_iterator const i_spin_end = end();
        for(; i_spin != i_spin_end; ++i_spin)
        {
          LADA_ASSERT(i_spin->site < lattice.sites.size(), "Index out-of-range.\n");
          atat::rVector3d const vec( center + i_spin->pos - lattice.sites[i_spin->site].pos );
          intermediate *= find_spin(vec, i_spin->site, _str, _map, _transform);
        }
        result += intermediate;
      }

      return result;
    }

    std::ostream& operator<<( std::ostream &_stream,  MLCluster::Spin const &_spin )
    {
      return _stream << "@" << _spin.site
                     << " " << _spin.pos(0) << " " << _spin.pos(1) << " " << _spin.pos(2);
    }
    
    std::ostream& operator<<( std::ostream &_stream,  MLCluster const &_cls )
    {
      _stream << std::fixed << std::setprecision(5) << std::setw(9);
      _stream << "Cluster@" << _cls.origin.site;
      
      if (_cls.size() == 0 ) return _stream << " J1\n";
      else _stream << "\n";

      MLCluster :: t_Spins :: const_iterator i_spin = _cls.begin();
      MLCluster :: t_Spins :: const_iterator const i_last = _cls.end();
      for ( ; i_spin != i_last; ++i_spin) _stream << " " << *i_spin << "\n";
      return _stream;
    }
  } // namespace CE
}

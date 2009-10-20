//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <set>
#include <iterator>
#include <algorithm>

#include "boost/lexical_cast.hpp"

#include <crystal/lattice.h>

#include "numeric_type.h"

namespace LaDa
{
  namespace enumeration
  {
    //! Create a flavor basis.
    boost::shared_ptr<FlavorBase> create_flavor_base( size_t _card, size_t _nflavor )
    {
      boost::shared_ptr<FlavorBase> result( new FlavorBase() );
      result->reserve( _card );
      t_uint k(1);
      for(size_t i(0); i < _card; ++i, k*=_nflavor) result->push_back( k );
      return result;
    }

    void check_size( size_t _card, size_t _nflavor )
    {
      t_uint first(1), second(_nflavor);
      for( size_t i(0); i < _card; ++i )
      {
        first *= _nflavor;
        second *= _nflavor;
        if( first > second ) BOOST_THROW_EXCEPTION( supercell_too_large() );
      } 
    }

    boost::shared_ptr<Database> create_database( size_t _card, size_t _nflavor )
    {
      check_size( _card, _nflavor );
      t_uint k = std::pow(_nflavor, _card);
      boost::shared_ptr<Database> result(new Database(k));
      result->flip();
      return result;
    }

    t_uint count_flavors( Crystal::Lattice const &_lattice)
    {
      std::set<Crystal::Lattice::t_Site::t_Type::value_type> set_;
      foreach( Crystal::Lattice::t_Site const &site, _lattice.sites )
        if( site.type.size() > 1 )
          std::copy( site.type.begin(), site.type.end(), std::inserter(set_, set_.begin()));
      return t_uint(set_.size());
    }

    std::string integer_to_bitstring( t_uint _x, FlavorBase const& _fl)
    {
      if( _fl.size() < 2 ) return "";
#     ifdef LADA_DEBUG
        if( _x >= _fl.back() * _fl[1] )
          BOOST_THROW_EXCEPTION( internal() << error_string("Argument _x is out of range.") );
#     endif

      std::string result;
      FlavorBase::const_reverse_iterator i_flavor = _fl.rbegin();
      FlavorBase::const_reverse_iterator i_flavor_end = _fl.rend();
      for(;i_flavor != i_flavor_end; ++i_flavor)
      {
        t_uint const flavor( _x / (*i_flavor) );
        _x %= (*i_flavor);

        result += boost::lexical_cast<std::string>(flavor);
      } // c
      return result;
    }

    void integer_to_structure( Crystal::Structure &_out,
                               t_uint _x, FlavorBase const &_fl, size_t const _card )
    {
      std::string const string( integer_to_bitstring(_x, _fl) );
      Crystal::t_SmithTransform const
        transform( Crystal::get_smith_transform(_out.lattice->cell, _out.cell) );
      Crystal::Structure::t_Atoms::iterator i_atom = _out.atoms.begin();
      Crystal::Structure::t_Atoms::iterator i_atom_end = _out.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom)
      {
        const size_t index = Crystal::get_linear_smith_index(transform, *i_atom);
        i_atom->type = (string[index] == '0') ? -1.0: 1.0;
      }
    }

    void integer_to_structure( Crystal::TStructure<std::string> &_out,
                               t_uint _x, FlavorBase const &_fl, size_t const _card )
    {
      LADA_ASSERT(_out.lattice != NULL, "lattice not set.")
      std::vector<size_t> bitstring;
      integer_to_vector(_x, _fl, bitstring);
      Crystal::t_SmithTransform const
        transform( Crystal::get_smith_transform(_out.lattice->cell, _out.cell) );
      Crystal::TStructure<std::string>::t_Atoms::iterator i_atom = _out.atoms.begin();
      Crystal::TStructure<std::string>::t_Atoms::iterator i_atom_end = _out.atoms.end();
      for(; i_atom != i_atom_end; ++i_atom)
      {
        LADA_ASSERT( i_atom->site != -1, "Site index not set.\n");
        size_t const index = Crystal::get_linear_smith_index(transform, *i_atom);
        size_t const type_index = bitstring[index];
        LADA_ASSERT( i_atom->site >= 0 and i_atom->site < _out.lattice->sites.size(),
                     "site index out of range.\n" ); 
        LADA_ASSERT(     type_index >= 0
                     and type_index < _out.lattice->sites[i_atom->site].type.size(),
                     "Type index out of range.\n" ); 
        i_atom->type = _out.lattice->sites[i_atom->site].type[type_index];
      }
    }
  }
}

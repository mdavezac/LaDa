#include "LaDaConfig.h"

#include <set>
#include <iterator>
#include <algorithm>

#include "boost/lexical_cast.hpp"

#include <crystal/lattice.h>
#include <opt/pow.h>

#include "numeric_type.h"

namespace LaDa
{
  namespace enumeration
  {
    //! Create a flavor basis.
    boost::shared_ptr<FlavorBase> create_flavor_base( size_t _card, size_t _nflavor )
    {
      if( not check_size( _card, _nflavor ) ) 
        BOOST_THROW_EXCEPTION( supercell_too_large() ); 
      boost::shared_ptr<FlavorBase> result( new FlavorBase() );
      result->reserve( _card );
      t_uint k(1);
      for(size_t i(0); i < _card; ++i, k*=_nflavor) result->push_back( k );
      return result;
    }

    boost::shared_ptr<Database> create_database( size_t _card, size_t _nflavor )
    {
      if( not check_size( _card, _nflavor ) ) 
        BOOST_THROW_EXCEPTION( supercell_too_large() ); 
      t_uint k = opt::pow(_nflavor, _card);
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


    template<class T_STRUCTURE, class T_CONVERT>
      void integer_to_structure_impl_( T_STRUCTURE &_out, t_uint _x, 
                                       FlavorBase const &_fl, T_CONVERT const &_conv )
      {
        LADA_ASSERT(_out.lattice != NULL, "Lattice not set.\n");
        size_t const N(_out.lattice->sites.size()); 
        std::vector<size_t> site_index_map;
        Crystal::Lattice::t_Sites::iterator i_site = _out.lattice->sites.begin();
        Crystal::Lattice::t_Sites::iterator i_site_end = _out.lattice->sites.end();
        for(size_t i(0); i_site != i_site_end; ++i_site)
          if( i_site->type.size() < 2 ) site_index_map.push_back(N);
          else
          {
            site_index_map.push_back(i);
            ++i;
          }
  
        std::vector<size_t> bitstring;
        integer_to_vector(_x, _fl, bitstring);
        Crystal::t_SmithTransform const
          transform( Crystal::get_smith_transform(_out.lattice->cell, _out.cell) );
        typename T_STRUCTURE::t_Atoms::iterator i_atom = _out.atoms.begin();
        typename T_STRUCTURE::t_Atoms::iterator i_atom_end = _out.atoms.end();
        for(; i_atom != i_atom_end; ++i_atom)
        {
          LADA_ASSERT( i_atom->site != -1, "Site index not set.\n");
          LADA_ASSERT( i_atom->site >= 0 and i_atom->site < N, "site index out of range.\n" ); 
          if( site_index_map[i_atom->site] == N ) continue;
  
          size_t const enum_index = site_index_map[i_atom->site];
          math::rVector3d const pos( i_atom->pos - _out.lattice->sites[i_atom->site].pos );
          size_t const index = Crystal::get_linear_smith_index(transform, enum_index, pos);
          LADA_ASSERT
          (
                bitstring[index] >= 0
            and bitstring[index] < _out.lattice->sites[i_atom->site].type.size(),
            "Index out of range.\n" 
          );
          _conv(*i_atom, bitstring[index]);
        }
      }

    struct Conv
    {
      Crystal::Lattice const * lattice;
      void operator()(Crystal::Structure::t_Atom &_out, size_t const _i) const 
        { _out.type = types::t_real( types::t_int(_i << 1) - 1); }
      void operator()(Crystal::TStructure<std::string>::t_Atom &_out, size_t const _i) const 
        { _out.type = lattice->sites[_out.site].type[_i]; }
    };
    void integer_to_structure( Crystal::Structure &_out,t_uint _x, FlavorBase const &_fl)
      { integer_to_structure_impl_( _out, _x, _fl, Conv() ); }
    void integer_to_structure( Crystal::TStructure<std::string> &_out,
                               t_uint _x, FlavorBase const &_fl)
    {
      Conv c; c.lattice = _out.lattice;
      integer_to_structure_impl_( _out, _x, _fl, c );
    }

  }
}

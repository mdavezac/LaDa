//
//  Version: $Id$
//
#ifndef LADA_ENUM_MUMERIC_TYPE_H_
#define LADA_ENUM_MUMERIC_TYPE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/dynamic_bitset.hpp>
#include <boost/shared_ptr.hpp>

#include <crystal/structure.h>
#include <crystal/smith.h>

#include "exceptions.h"

namespace LaDa
{
  //! \cond
  namespace Crystal { class Lattice; }
  //! \endcond

  namespace enumeration
  {
    //! A database of structures.
    typedef boost::dynamic_bitset<long long unsigned int> Database;

    //! \brief Type of the database integers.
    //! \warning We cannot consider spaces larger than this type can hold.
    //!          In other words, this unsigned integer must be able to
    //!          container (nflavors) ^ (nsites in supercell).
    typedef Database :: size_type t_uint;

    //! Type of the elements DxG.
    typedef std::pair<size_t, math::iVector3d> t_Element;

    //! A base of nflavor^m.
    typedef std::vector<t_uint> FlavorBase;


    //! Create a flavor basis.
    boost::shared_ptr<FlavorBase> create_flavor_base( size_t _card, size_t _nflavor );

    inline t_uint get_index(size_t const &_d, math::iVector3d const &_g, 
                            math::iVector3d const &_smith, size_t const &_card) throw()
      { return _card-1- Crystal::get_linear_smith_index(_smith, _d, _g); }

    //! Throws when supercell is too large for t_uint integer type.
    void check_size( size_t _card, size_t _nflavor );

    //! Returns an initialized database for specified size. 
    boost::shared_ptr<Database> create_database( size_t _card, size_t _nflavor );

    //! Returns the number of stomic species in the lattice.
    t_uint count_flavors( Crystal::Lattice const &_lattice);

    //! Transforms an integer to a string.
    std::string integer_to_bitstring( t_uint _x, FlavorBase const& _fl);
    //! Transforms an integer to a string.
    template< class T_VECTOR>
      void integer_to_vector(t_uint _x, FlavorBase const& _fl, T_VECTOR &_out)
      {
        if( _fl.size() < 2 ) return;
#       ifdef LADA_DEBUG
          if( _x >= _fl.back() * _fl[1] )
            BOOST_THROW_EXCEPTION( internal() << error_string("Argument _x is out of range.") );
#       endif
        
        _out.resize( _fl.size() );
        FlavorBase::const_reverse_iterator i_flavor = _fl.rbegin();
        FlavorBase::const_reverse_iterator i_flavor_end = _fl.rend();
        typename T_VECTOR::iterator i_result = _out.begin();
        for(;i_flavor != i_flavor_end; ++i_flavor, ++i_result)
        {
          *i_result = _x / (*i_flavor);
          _x %= (*i_flavor);
        }
      }

    //! Transform an integer to a structure.
    void integer_to_structure( Crystal::Structure &_out,
                               t_uint _x, FlavorBase const &_fl);
    //! Transform an integer to a structure.
    void integer_to_structure( Crystal::TStructure<std::string> &_out,
                               t_uint _x, FlavorBase const &_fl);
  }
}

#endif

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

#include <atat/vectmac.h>

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
    typedef std::pair<size_t, atat::iVector3d> t_Element;

    //! A base of nflavor^m.
    typedef std::vector<t_uint> FlavorBase;


    //! Create a flavor basis.
    boost::shared_ptr<FlavorBase> create_flavor_base( size_t _card, size_t _nflavor );

    inline t_uint get_index(size_t const &_d, atat::iVector3d const &_g, 
                            atat::iVector3d const &_smith, size_t const &_card) throw()
      { return _card-1- _g(2) - _smith(2) * ( _g(1) + _smith(1) * ( _g(0) + _smith(0) *_d ) ); }

    //! Throws when supercell is too large for t_uint integer type.
    void check_size( size_t _card, size_t _nflavor ) throw(boost::exception);

    //! Returns an initialized database for specified size. 
    boost::shared_ptr<Database> create_database( size_t _card, size_t _nflavor );

    //! Returns the number of stomic species in the lattice.
    t_uint count_flavors( Crystal::Lattice const &_lattice);

    //! Transforms an integer to a string.
    std::string integer_to_bitstring( t_uint _x, FlavorBase const& _fl);

  }
}

#endif

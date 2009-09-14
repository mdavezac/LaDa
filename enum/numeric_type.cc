//
//  Version: $Id$
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <set>
#include <iterator>
#include <algorithm>

#include <crystal/lattice.h>

#include "exceptions.h"
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

    void check_size( size_t _card, size_t _nflavor ) throw(boost::exception)
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
      return boost::shared_ptr<Database>(new Database(std::pow(_nflavor, _card), true));
    }

    t_uint count_flavors( Crystal::Lattice const &_lattice)
    {
      std::set<Crystal::Lattice::t_Site::t_Type::value_type> set_;
      foreach( Crystal::Lattice::t_Site const &site, _lattice.sites )
        std::copy( site.type.begin(), site.type.end(), std::inserter(set_, set_.begin()));
      return t_uint(set_.size());
    }

  }
}

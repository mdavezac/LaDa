
//
//  Version: $Id$
//
#ifndef _LADA_COMPARE_SITES_H_
#define _LADA_COMPARE_SITES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <set>
#include "lattice.h"


namespace LaDa
{
  namespace Crystal 
  {

    struct CompareSites
    {
      CompareSites( Lattice::t_Site const &_site )
      {
        pos = _site.pos;
        std::copy(_site.type.begin(), _site.type.end(), std::inserter(set_, set_.begin()));
      }
      CompareSites( CompareSites const &_c ): set_(_c.set_) {}
      bool operator()( Lattice::t_Site::t_Type const &_types ) const
      {
        Lattice::t_Site::t_Type::const_iterator i_first = _types.begin();
        Lattice::t_Site::t_Type::const_iterator const i_end = _types.end();
        for(; i_first != i_end; ++i_first )
          if( set_.find(*i_first) == set_.end() ) break;
        return i_first == i_end;
      }
      bool operator()( Lattice::t_Site const& _site ) const
      {
        return not(    Fuzzy::neq(_site.pos(0), pos(0)) 
                    or Fuzzy::neq(_site.pos(1), pos(1)) 
                    or Fuzzy::neq(_site.pos(2), pos(2)) );
      }

      std::set<Lattice::t_Site::t_Type::value_type> set_;
      atat::rVector3d pos;
    };


  } // namespace Crystal
} // namespace LaDa
#endif

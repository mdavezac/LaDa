
//
//  Version: $Id$
//
#ifndef _LADA_COMPARE_SITES_H_
#define _LADA_COMPARE_SITES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <set>
#include <opt/types.h>
#include "lattice.h"


namespace LaDa
{
  namespace Crystal 
  {

    struct CompareSites
    {
      CompareSites( Lattice::t_Site const &_site, types::t_real _tolerance = -1e0 )
      {
        tolerance = _tolerance < 0e0 ? types::tolerance: _tolerance;
        pos = _site.pos;
        std::copy(_site.type.begin(), _site.type.end(), std::inserter(set_, set_.begin()));
      }
      CompareSites( CompareSites const &_c ): set_(_c.set_), pos(_c.pos), tolerance(_c.tolerance) {}
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
        return     std::abs(_site.pos(0)-pos(0)) < tolerance 
               and std::abs(_site.pos(1)-pos(1)) < tolerance 
               and std::abs(_site.pos(2)-pos(2)) < tolerance;
      }

      std::set<Lattice::t_Site::t_Type::value_type> set_;
      math::rVector3d pos;
      types::t_real tolerance;
    };


  } // namespace Crystal
} // namespace LaDa
#endif

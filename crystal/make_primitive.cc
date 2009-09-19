//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cmath>

#include <atat/is_int.h>
#include <opt/fuzzy.h>

#include "lattice.h"
#include "compare_sites.h"


namespace LaDa
{
  namespace Crystal 
  {
    atat::rVector3d into_cell( atat::rVector3d const &_vec, 
                               atat::rMatrix3d const &_cell, 
                               atat::rMatrix3d const &_inv)
    {
      atat::rVector3d result( _inv * _vec );
      result(0) -= std::floor(result(0));
      result(1) -= std::floor(result(1));
      result(2) -= std::floor(result(2));
      return _cell * result;
    }

    bool Lattice::make_primitive()
    {
      if( sites.size() == 0 ) return true;
      // copies lattice.
      Lattice copy(*this);
      bool is_primitive = true;

      // moves sites into unit-cell.
      atat::rMatrix3d const inv( !cell );
      foreach(t_Site &_site, copy.sites) into_cell(_site.pos, cell, inv);

      // Then compares fractional translations from site 0 to sites of same type.
      std::vector<atat::rVector3d> translations;
      CompareSites compsites(copy.sites.front());
      t_Sites :: const_iterator i_site = copy.sites.begin();
      t_Sites :: const_iterator const i_site_end = copy.sites.end();
      for(i_site; i_site != i_site_end; ++i_site )
      {
        // Translations are created from equivalent sites only.
        if( not compsites(i_site->type) ) continue;

        // creates translation.
        atat::rVector3d translation( into_cell(i_site->pos - compsites.pos, cell, inv) );

        // checks that it leaves the lattice invariant.
        bool has_mapping(true);
        foreach( t_Site const &site, copy.sites )
        {
          CompareSites check_translation(site);
          check_translation.pos = into_cell( site.pos + translation, cell, inv );

          t_Sites::const_iterator i_found
            = std::find_if( copy.sites.begin(), copy.sites.end(), check_translation );

          if( i_found != copy.sites.end() ) continue;
          has_mapping = false;
          break;
        }

        if( not has_mapping ) continue;

        // adds translation to vector. This lattice is not primitive.
        translations.push_back( translation );
        is_primitive = false;
      }

      // This lattice is primitive.
      if( is_primitive ) return true;

      // adds original translations.
      translations.push_back( cell.get_column(0) );
      translations.push_back( cell.get_column(1) );
      translations.push_back( cell.get_column(2) );

      // Loops over possible primitive cells.
      typedef std::vector<atat::rVector3d> :: const_iterator t_cit;
      t_cit const i_vec_begin( translations.begin() );
      t_cit const i_vec_end( translations.end() );
      atat::rMatrix3d new_cell = cell;
      types::t_real volume = std::abs(atat::det(new_cell));
      for( t_cit i_first(i_vec_begin); i_first != i_vec_end; ++i_first )
        for( t_cit i_second(i_vec_begin); i_second != i_vec_end; ++i_second )
        {
          if( i_first == i_second ) continue;
          for( t_cit i_third(i_vec_begin); i_third != i_vec_end; ++i_third )
          {
            if( i_first == i_third or i_second == i_third ) continue;
            // construct new cell.
            atat::rMatrix3d trial;
            trial.set_column(0, *i_first);
            trial.set_column(1, *i_second);
            trial.set_column(2, *i_third);

            // singular matrix?
            types::t_real const det( atat::det(trial) );
            if( Fuzzy::is_zero(det) ) continue;
            // Volume smaller than current new_cell?
            if( Fuzzy::gt( std::abs(det), volume) ) continue;
            // Direct matrix?
            if( det < 0e0 )
            {
              trial.set_column(2, *i_second);
              trial.set_column(3, *i_third);
              LADA_ASSERT(atat::det(trial) > 0, "Shouldn't happen.\n");
            }

            // Checks that all lattice sites are integers.
            atat::rMatrix3d const trial_inv( !trial );
            bool all_integer = true;
            foreach( t_Site const site, copy.sites )
            {
              atat::rVector3d vec( trial_inv * site.pos );
              if( atat::is_integer(vec) ) continue;
              all_integer = false;
              break;
            }
            if( not all_integer ) continue;
            volume = std::abs(det);
            new_cell = trial;
          }
        }

      // Found the new cell with smallest volume (e.g. primivite)
      LADA_ASSERT
      (
        not Fuzzy::is_zero( volume - atat::det(new_cell) ), 
        "Could not find primitive cell.\n"
      );

      // now creates new lattice.
      sites.clear();
      cell = new_cell;
      i_site = copy.sites.begin();
      sites.push_back(*i_site);
      for(; i_site != i_site_end; ++i_site)
      {
        compsites = *i_site;
        t_Sites::const_iterator i_found = std::find_if(sites.begin(), sites.end(), compsites);
        if( i_found == sites.end() ) sites.push_back(*i_site);
      }
      
      return true;
    }

  } // namespace Crystal

} // namespace LaDa

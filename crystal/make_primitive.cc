#include "LaDaConfig.h"

#include <cmath>

#include <Eigen/LU> 

#include <math/misc.h>
#include <math/fuzzy.h>
#include <limits>

#include <math/fuzzy.h>
#include "lattice.h"
#include "compare_sites.h"


namespace LaDa
{
  namespace Crystal 
  {
    bool Lattice::make_primitive(types::t_real _tolerance)
    {
      if( _tolerance < 0e0 ) _tolerance = types::tolerance;
      if( sites.size() == 0 ) return true;
      // copies lattice.
      Lattice copy(*this);
      bool is_primitive = true;

      // moves sites into unit-cell.
      math::rMatrix3d const inv(cell.inverse());
      foreach(t_Site &_site, copy.sites) _site.pos = into_cell(_site.pos, cell, inv);

      // Then compares fractional translations from site 0 to sites of same type.
      std::vector<math::rVector3d> translations;
      CompareSites compsites(copy.sites.front(), _tolerance);
      t_Sites :: const_iterator i_site = copy.sites.begin();
      t_Sites :: const_iterator const i_site_end = copy.sites.end();
      for(; i_site != i_site_end; ++i_site )
      {
        // Translations are created from equivalent sites only.
        if( not compsites(i_site->type) ) continue;

        // creates translation.
        math::rVector3d const translation( into_cell(i_site->pos - compsites.pos, cell, inv) );
        
        // loop on null translation.
        if(     math::is_null(translation(0))
            and math::is_null(translation(1)) 
            and math::is_null(translation(2)) ) continue;

        // checks that it leaves the lattice invariant.
        t_Sites :: const_iterator i_mapping = copy.sites.begin();
        for(; i_mapping != i_site_end; ++i_mapping)
        {
          CompareSites check_translation(*i_mapping, _tolerance);
          check_translation.pos = into_cell( i_mapping->pos + translation, cell, inv );

          t_Sites::const_iterator i_found
            = std::find_if( copy.sites.begin(), copy.sites.end(), check_translation );

          if(i_found == copy.sites.end()) break;
          if(not check_translation(i_found->type)) break;
        }

        if( i_mapping != i_site_end ) continue;

        // adds translation to vector. This lattice is not primitive.
        translations.push_back( translation );
        is_primitive = false;
      }

      // This lattice is primitive.
      if( is_primitive ) return true;

      // adds original translations.
      translations.push_back( cell.col(0) );
      translations.push_back( cell.col(1) );
      translations.push_back( cell.col(2) );

      // Loops over possible primitive cells.
      typedef std::vector<math::rVector3d> :: const_iterator t_cit;
      t_cit const i_vec_begin( translations.begin() );
      t_cit const i_vec_end( translations.end() );
      math::rMatrix3d new_cell = cell;
      types::t_real volume = std::abs(new_cell.determinant());
      for( t_cit i_first(i_vec_begin); i_first != i_vec_end; ++i_first )
        for( t_cit i_second(i_vec_begin); i_second != i_vec_end; ++i_second )
        {
          if( i_first == i_second ) continue;
          for( t_cit i_third(i_vec_begin); i_third != i_vec_end; ++i_third )
          {
            if( i_first == i_third or i_second == i_third ) continue;
            // construct new cell.
            math::rMatrix3d trial;
            trial.col(0) = *i_first;
            trial.col(1) = *i_second;
            trial.col(2) = *i_third;

            // singular matrix?
            types::t_real const det(trial.determinant());
            if( math::is_null(det) ) continue;
            // Volume smaller than current new_cell?
            if( math::geq( std::abs(det), volume) ) continue;
            // Direct matrix?
            if( det < 0e0 )
            {
              trial.col(2) = *i_second;
              trial.col(1) = *i_third;
              LADA_ASSERT(trial.determinant() > 0, "Shouldn't happen.\n");
            }
            // Checks that original cell is a supercell.
            if( not math::is_integer( ((!new_cell) * cell).eval() ) ) continue;

            // Checks that all lattice sites are integers.
            volume = std::abs(det);
            new_cell = trial;
          }
        }

      // Found the new cell with smallest volume (e.g. primivite)
      LADA_ASSERT
      (
        not math::is_null(cell.determinant() - new_cell.determinant() ), 
        "Could not find primitive cell.\n"
      );

      // now creates new lattice.
      sites.clear();
      cell = new_cell;
      math::rMatrix3d const inv_cell(cell.inverse());
      foreach(t_Site &s, copy.sites) s.pos = into_cell(s.pos, cell, inv_cell);
      sites.push_back(copy.sites.front());
      for(i_site = copy.sites.begin(); i_site != i_site_end; ++i_site)
      {
        CompareSites check(*i_site);
        t_Sites::const_iterator i_found = std::find_if(sites.begin(), sites.end(), check);
        if( i_found == sites.end() ) sites.push_back(*i_site);
        else 
        { LADA_DOASSERT( check(i_found->type), "Two inequivalent sites at same position.\n" ); }
      }

      while(not make_primitive());
      
      return false;
    }

  } // namespace Crystal

} // namespace LaDa

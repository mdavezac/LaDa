#ifndef LADA_CRYSTAL_PRIMITIVE_H
#define LADA_CRYSTAL_PRIMITIVE_H
#include "LaDaConfig.h"

#include <cmath>

#include <Eigen/LU> 

#include <math/fuzzy.h>
#include <math/gruber.h>

#include "structure.h"
#include "compare_sites.h"
#include "exceptions.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {
    //! Returns the primitive unit structure. 
    template<class T_TYPE>
      TemplateStructure<T_TYPE> primitive( TemplateStructure<T_TYPE> const &_structure,
                                           types::t_real _tolerance = -1e0 )
      {
        if( _tolerance < 0e0 ) _tolerance = types::tolerance;
        if( _structure.size() == 0 ) BOOST_THROW_EXCEPTION(error::empty_structure());

        // copies lattice.
        typedef TemplateStructure<T_TYPE> t_Sites;
        TemplateStructure<T_TYPE> result(_structure.copy());
        math::rMatrix3d cell = math::gruber(result.cell());
        bool is_primitive = true;
  
        // moves sites into unit-cell.
        math::rMatrix3d const inv(cell.inverse());
        typename TemplateStructure<T_TYPE>::iterator i_atom = result.begin();
        typename TemplateStructure<T_TYPE>::iterator const i_atom_end = result.end();
        for(; i_atom != i_atom_end; ++i_atom)
          i_atom->pos() = into_cell(i_atom->pos(), cell, inv);
  
        // Then compares fractional translations from site 0 to sites of same type.
        std::vector<math::rVector3d> translations;
        CompareOccupations<T_TYPE> const compsites(result.front()->type);
        math::rVector3d const center = result.front()->pos;
        typename t_Sites :: const_iterator i_site = _structure.begin();
        typename t_Sites :: const_iterator const i_site_end = _structure.end();
        for(; i_site != i_site_end; ++i_site )
        {
          // Translations are created from equivalent sites only.
          if( not compsites(i_site->type()) ) continue;
  
          // creates translation.
          math::rVector3d const translation = into_voronoi(i_site->pos() - center, cell, inv);
          
          // loop on null translation.
          if( math::is_null(translation, _tolerance) ) continue;
  
          // checks that it leaves the lattice invariant.
          typename t_Sites :: const_iterator i_mapping = _structure.begin();
          typename t_Sites :: const_iterator const i_fend = result.end(); 
          for(; i_mapping != i_site_end; ++i_mapping)
          {
            math::rVector3d const pos = into_cell(i_mapping->pos() + translation, cell, inv);
            CompareSites<T_TYPE> const cmp(pos, i_mapping->type(), _tolerance);
            typename t_Sites::iterator i_found = result.begin(); 
            for(; i_found != i_fend; ++i_found) if(cmp(*i_found)) break;
            if(i_found == i_fend) break;
          }
  
          if( i_mapping != i_site_end ) continue;
  
          // adds translation to vector. This lattice is not primitive.
          translations.push_back(into_voronoi(translation, cell));
          is_primitive = false;
        }
  
        // This lattice is primitive.
        if( is_primitive ) return _structure;
  
        // adds original translations.
        translations.push_back( cell.col(0) );
        translations.push_back( cell.col(1) );
        translations.push_back( cell.col(2) );
  
        // Loops over possible primitive cells.
        typedef std::vector<math::rVector3d> :: const_iterator t_cit;
        t_cit const i_vec_begin( translations.begin() );
        t_cit const i_vec_end( translations.end() );
        math::rMatrix3d new_cell = result.cell();
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
              if( math::is_null(det, 3e0*_tolerance) ) continue;
              // Volume smaller than current new_cell?
              if( math::geq(std::abs(det), volume, 3e0 * _tolerance) ) continue;
              // Direct matrix?
              if( det < 0e0 )
              {
                trial.col(2) = *i_second;
                trial.col(1) = *i_third;
#               ifdef LADA_DEBUG
                  if(trial.determinant() < types::tolerance)
                    BOOST_THROW_EXCEPTION(error::internal() << error::string("Negative volume."));
#               endif
              }
              // Checks that original cell is a supercell.
              if( not math::is_integer(trial.inverse() * result.cell(), _tolerance) ) continue;
  
              // Checks that all lattice sites are integers.
              volume = std::abs(det);
              new_cell = trial;
            }
          }
  
        // Found the new cell with smallest volume (e.g. primivite)
        if(math::eq(_structure.volume(), new_cell.determinant()))
          BOOST_THROW_EXCEPTION(
              error::internal() << error::string("Found translation but no primitive cell."));
  
        // now creates new lattice.
        result.clear();
        result.cell() = math::gruber(new_cell);
        math::rMatrix3d const inv_cell(result.cell().inverse());
        for(i_site = _structure.begin(); i_site != i_site_end; ++i_site)
        {
          math::rVector3d const pos = into_cell(i_site->pos(), result.cell(), inv_cell); 
          CompareSites<T_TYPE> const check(pos, i_site->type(), _tolerance);
          typename t_Sites::const_iterator i_found = result.begin(); 
          typename t_Sites::const_iterator const i_fend = result.end(); 
          for(; i_found != i_fend; ++i_found)
          {
            if(check(*i_found)) break;
          }
          if( i_found == i_fend )
          {
            result.push_back(i_site->copy());
            result.back()->pos = pos;
          }
        }
        if(_structure.size() % result.size() != 0)
          BOOST_THROW_EXCEPTION
          ( 
             error::internal() 
               << error::string("Nb of atoms in output not multiple of input.")
          );
        if(math::neq(types::t_real(_structure.size()/result.size()), _structure.volume()/result.volume()))
          BOOST_THROW_EXCEPTION(error::internal() << error::string("Size and volumes do not match."));
  
        return result;
      }

  } // namespace crystal
} // namespace LaDa
#endif

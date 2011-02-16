#include "LaDaConfig.h"

#include <cstdlib>

#include <algorithm>
#include <functional>
#include <boost/filesystem/operations.hpp>

#include <physics/physics.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>
#include <opt/ndim_iterator.h>

#include "atomic_center.h"
  
namespace LaDa
{
  namespace vff
  { 
    types::t_int AtomicCenter :: add_bond( t_BondRefd _bond, 
                                            const types::t_real _cutoff ) 
    {
      bool found_bond = false;
      opt::NDimIterator< types::t_int, std::less_equal<types::t_int> > period;
      period.add(-1,1);
      period.add(-1,1);
      period.add(-1,1);

      do // goes over periodic image
      {
        // constructs perdiodic image of atom *i_bond
        math::rVector3d frac_image, image;
        frac_image[0] =  (types::t_real) period.access(0);
        frac_image[1] =  (types::t_real) period.access(1);
        frac_image[2] =  (types::t_real) period.access(2);
        image = _bond->i_atom_->pos + structure->cell * frac_image;

        // checks if within 
        if( (image - i_atom_->pos).squaredNorm() < _cutoff  )
        {
          bonds.push_back( _bond );
          translations.push_back( frac_image );
          do_translates.push_back( not math::is_zero(frac_image.squaredNorm()) );
          found_bond = true;
        }
      } while ( ++period ); 

      return found_bond ? (types::t_int) bonds.size() : -1; // not a bond
    }
#   ifdef LADA_DEBUG
      void AtomicCenter :: const_iterator :: check() const
      {
        LADA_NASSERT(not parent,
                 "Pointer to parent atom is invalid.\n")
        LADA_NASSERT( parent->bonds.size() == 0,
                  "The number of bond is zero.\n")
        LADA_NASSERT( parent->translations.size() == 0,
                  "The number of translations is zero.\n")
        LADA_NASSERT( parent->do_translates.size() == 0,
                  "The number of translation switches is zero.\n")
        LADA_NASSERT( i_bond - parent->bonds.end() > 0,
                  "The bond iterator is beyond the last bond.\n")
        if( i_bond ==  parent->bonds.end() ) return;
        LADA_NASSERT( i_bond - parent->bonds.begin() < 0,
                  "The bond iterator is before the first bond.\n")
        types::t_int pos = i_bond -  parent->bonds.begin();
        LADA_NASSERT( i_translation - parent->translations.end() > 0,
                  "The translation iterator is beyond the last bond.\n")
        LADA_NASSERT( i_translation - parent->translations.begin() < 0,
                  "The translation iterator is before the first bond.\n")
        LADA_NASSERT( i_translation != parent->translations.begin() + pos,
                     "The bond iterator and the translation "
                  << "iterator are out of sync.\n")
        LADA_NASSERT( i_do_translate - parent->do_translates.end() > 0,
                  "The do_translate iterator is beyond the last bond.\n")
        LADA_NASSERT( i_do_translate - parent->do_translates.begin() < 0,
                  "The do_translate iterator is before the first bond.\n")
        LADA_NASSERT( i_do_translate != parent->do_translates.begin() + pos,
                     "The bond iterator and the do_translate "
                  << "iterator are out of sync.\n")
      }
      void AtomicCenter :: const_iterator :: check_valid() const
      {
        check();
        LADA_NASSERT(i_bond == parent->bonds.end(), "Invalid iterator.\n";);
        LADA_NASSERT(not parent->i_atom_, "Origin of the parent atom is invalid.\n");
        LADA_NASSERT(not (*i_bond)->i_atom_, "Origin of the bond atom is invalid.\n");
      }
#   endif
  } // namespace vff
} // namespace LaDa

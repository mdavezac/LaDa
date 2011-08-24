#ifndef LADA_CRYSTAL_UTILITIES_H
#define LADA_CRYSTAL_UTILITIES_H

#include "LaDaConfig.h"


#include <opt/types.h>
#include <math/misc.h>
#include "atom.h"

namespace LaDa
{
  namespace crystal 
  {

    //! Refolds a periodic vector into the unit cell.
    math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell, 
                               math::rMatrix3d const &_inv);
    //! Refolds a periodic vector into the unit cell.
    inline math::rVector3d into_cell( math::rVector3d const &_vec, 
                               math::rMatrix3d const &_cell )
      { return into_cell(_vec, _cell, _cell.inverse()); }
    
   
    //! Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell, 
                                  math::rMatrix3d const &_inv);
    //! \brief Refolds a periodic vector into the voronoi cell (eg first BZ or WignerSeitz).
    //! \details May fail if the matrix is a weird parameterization of the
    //!          lattice. It is best to use a grubber(...) cell . 
    inline math::rVector3d into_voronoi( math::rVector3d const &_vec, 
                                  math::rMatrix3d const &_cell )
      { return into_voronoi(_vec, _cell, _cell.inverse()); }

    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell, 
                                   math::rMatrix3d const &_inv);
    //! \brief Refolds a periodic vector into a cell centered around zero (in
    //!        fractional coordinates).
    //! \details Since the vector is refolded in fractional coordinates, it may
    //!          or may not be the vector with smallest norm. Use math::rVector3d
    //!          into_voronoi() to get the equivalent vector with smallest norm.
    inline math::rVector3d zero_centered( math::rVector3d const &_vec, 
                                   math::rMatrix3d const &_cell )
      { return zero_centered(_vec, _cell, _cell.inverse()); }

    //! \brief Compares species and positions in a structure.
    //! \details Expects the atomic types to be sequences.
    template<class T_TYPE>
      struct CompareSites
      {
        //! \typedef the type of the atoms held by a structure.
        typedef Atom<T_TYPE> t_Atom;
        //! \brief Constructor
        //! \param[in] _site Site containing the species against which to check.
        //! \param[in] _tolerance when comparing atomic positions. If negative,
        //!                       uses standard tolerance defined in opt/types.h.
        CompareSites(t_Atom const &_site, types::t_real _tolerance = -1e0 )
        {
          tolerance_ = _tolerance < 0e0 ? types::tolerance: _tolerance;
          pos_ = _site.pos;
          std::copy(_site.type.begin(), _site.type.end(), std::inserter(set_, set_.begin()));
        }
        //! Copy constructor
        CompareSites(CompareSites const &_c): set_(_c.set_), pos_(_c.pos_), tolerance_(_c.tolerance_) {}
        //! Returns true if all species in the container are in the input set.
        template<class T_CONTAINER>
          bool operator()(T_CONTAINER const &_types ) const
          {
            typename T_CONTAINER::const_iterator i_first = _types.begin();
            typename T_CONTAINER::const_iterator const i_end = _types.end();
            for(; i_first != i_end; ++i_first )
              if( set_.find(*i_first) == set_.end() ) return false;
            return true;
          }
        //! Compares the position of site to input site.
        bool operator()( t_Atom const& _site ) const
          { return math::eq(_site.pos, pos_, tolerance_); }
      
        private:
          //! \var set containing input species.
          std::set<typename T_TYPE::value_type> set_;
          //! \var input position.
          math::rVector3d pos_;
          //! \var tolerance when comparing positions.
          types::t_real tolerance_;
      };

  } // namespace Crystal
} // namespace LaDa
  
#endif

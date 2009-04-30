//
//  Version: $Id: ideal_lattice.h 1075 2009-04-29 17:34:00Z davezac $
//

#ifndef _LADA_CRYSTAL_DIVIDE_AND_CONQUER_H_
#define _LADA_CRYSTAL_DIVIDE_AND_CONQUER_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <boost/tuple/tuple.hpp>

#include <opt/types.h>

#include "lattice.h"
#include "structure.h"


namespace LaDa 
{
  namespace Crystal 
  {
    template< class T_TYPE >
      class ConquerBox
      {
        public:
          //! Type of the divide and conquer box, with cell and origin.
          typedef boost::tuples::tuple< atat::rMatrix3d, atat::rVector3d > t_Box;
          //! Returns the cell of the divide and conquer box.
          const atat::rMatrix3d& cell() const { return bt::get<0>(box_); }
          //! Returns the cell origin of the divide and conquer box.
          const atat::rVector3d& origin() const { return bt::get<1>(box_); }
          //! Returns the index of and atom inside this box, w.r.t. the original structure.
          size_t index( size_t _i ) const { return indices[_i]; }
          //! \brief Returns true if \a _vec points inside the box.
          //! \details \a _vec is in cartesian coordinates.
          bool is_inside_box_cartesian( const atat::rVector3d &_vec ) const;
          //! \brief Returns true if \a _vec points inside the box.
          //! \details \a _vec is in fractional coordinates.
          bool is_inside_box_fractional( const atat::rVector3d &_vec ) const;

        protected:
          t_Box overlap_box_;
          t_Box box_;
          Crystal::Structure &structure_;
          std::vector< const Crystal::Structure<T_TYPE>::t_Atom& > atoms_;
          std::vector< size_t > indices_;
      };

    //! Type of the return of the divide and conquer method.
    template< class T_TYPE >
      struct t_ConquerBoxes 
      {
        //! Type of the container returned of the divide and conquer method.
        typedef std::vector< ConquerBox<T_TYPE> > type;
        //! Type of the return of the divide and conquer method. A shared pointer.
        typedef std::vector< ConquerBox<T_TYPE> > shared_ptr;
      };
    
    //! \brief Returns divide and conquer boxes.
    //! \details The return should not exceed the lifetime of \a _structure. 
    //!          In fact, \a _structure should not be changed at all, until the
    //!          return is destroyed. 
    //! \params[in] \a _structure is the structure to divide into manageable subsets.
    //! \params[in] \a _n is a vector of integers denoting how to divide each
    //!                   axis of \a _structure.
    //! \params[in] \a _overlap_distance is the size of the larger box for
    //!                                  which to include further atoms.
    //!                                 (Cartesian coordinates).
    typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::Structure &_structure, 
                                const atat::iVector3d &_n,
                                const types::t_real _overlap_distance );
                           

  } // namespace Crystal

} // namespace LaDa

#include "divide_and_conquer.impl.h"

#endif

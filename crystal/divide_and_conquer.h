#ifndef _LADA_CRYSTAL_DIVIDE_AND_CONQUER_H_
#define _LADA_CRYSTAL_DIVIDE_AND_CONQUER_H_

#include "LaDaConfig.h"

#include <vector>

#include <boost/tuple/tuple.hpp>
#include <boost/shared_ptr.hpp>

#include <opt/types.h>

#include "lattice.h"
#include "structure.h"


namespace LaDa 
{
  namespace Crystal 
  {
    //! \cond 
    template< class T_TYPE > class ConquerBox;
    //! \endcond

    //! Type of the return of the divide and conquer method.
    template< class T_TYPE >
      struct t_ConquerBoxes 
      {
        //! Type of the container returned of the divide and conquer method.
        typedef std::vector< ConquerBox<T_TYPE> > type;
        //! Type of the return of the divide and conquer method. A shared pointer.
        typedef boost::shared_ptr< type > shared_ptr;
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
    template< class T_TYPE > typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                                const math::iVector3d &_n,
                                const types::t_real _overlap_distance );
                           
    //! \brief Tries to guess adequate parameters for divide and conquer.
    //! \details the conditions are that the number of boxes is a multiple of
    //!          the number of processor, the number of atoms per box, the
    //!          squatness of the box.
    template< class T_TYPE > 
      math::iVector3d guess_dnc_params( const Crystal::TStructure<T_TYPE> &_structure, 
                                        size_t _nperbox );
    
    //! A small box for divide and conquer algorithms.
    template< class T_TYPE >
      class ConquerBox
      {
          friend  typename t_ConquerBoxes<T_TYPE>::shared_ptr  
            divide_and_conquer_boxes<T_TYPE>( const Crystal::TStructure<T_TYPE> &_structure, 
                                              const math::iVector3d &_n,
                                              const types::t_real _overlap_distance );
          //! Type of the atoms.
          typedef typename Crystal::TStructure<T_TYPE>::t_Atoms t_Atoms;
          //! Type of an atom.
          typedef typename Crystal::TStructure<T_TYPE>::t_Atom t_Atom;
        public:
          //! Type of the divide and conquer box, with cell, origin, and overlap.
          typedef boost::tuples::tuple
          < 
            math::rMatrix3d, 
            math::rVector3d,
            math::rVector3d
          > t_Box;
          //! \brief Type of the atomic state.
          //! \detail first item contains index of atom in structure.
          //!         second item is true if atom is in small box.
          typedef boost::tuples::tuple< size_t, bool > t_State;
          //! A vector of states.
          typedef std::vector< t_State > t_States;

          //! Returns the cell of the divide and conquer box.
          const math::rMatrix3d& cell() const { return boost::tuples::get<0>(box_); }
          //! Returns the cell origin of the divide and conquer box.
          const math::rVector3d& origin() const { return boost::tuples::get<1>(box_); }
          //! Returns the index of and atom inside this box, w.r.t. the original structure.
          const t_State& state( size_t _i ) const
          { 
            LADA_NASSERT( _i < states_.size(), "Index out-of-range.\n" )
            return states_[_i]; 
          }
          //! Iterator to the first state in the (large) box.
          typename t_States::const_iterator begin() const { return states_.begin(); }
          //! Iterator to the end of the states in the (large) box.
          typename t_States::const_iterator end() const { return states_.end(); }
          //! Number of states in the (large) box.
          size_t size() const { return states_.size(); }
          //! checks whether a state is already counted
          bool is_counted(size_t i) const
          {
            foreach(t_State const &state, states_)
              if( boost::tuples::get<0>(state) == i ) return true;
            return false;
          }

        protected:
          //! Holds boxes sizes and position.
          t_Box box_;
          //! Holds states located in the box.
          t_States states_;
      };


  } // namespace Crystal

} // namespace LaDa

#include "divide_and_conquer.impl.h"

#endif

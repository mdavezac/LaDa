#ifndef LADA_CRYSTAL_PERIODIC_DNC_H
#define LADA_CRYSTAL_PERIODIC_DNC_H

#include "LaDaConfig.h"

#include <vector>

#include <boost/shared_ptr.hpp>

#include <opt/types.h>

#include "lattice.h"
#include "structure.h"


namespace LaDa 
{
  namespace Crystal 
  {

    //! Container of divide and conquer boxes.
    class DnCBoxes
    {
      public:
        //! References a point in a periodic divide and conquer box.
        struct Point
        {
          //! index in given structure.
          size_t index; 
          //! Translation in cartesian coordinates.
          math::rVector3d translation;
          //! If true, then position is in small box.
          bool in_small_box;
        };
        //! Defines a divide and conquer box.
        typedef std::vector<Point> t_Box;
        //! Defines a divide and conquer box.
        typedef std::vector<t_Box> t_Boxes;

        //! Does not initialize to avoid throwing.
        DnCBoxes() : n_(0,0,0) {}

        //! Creates box.
        template<class T_TYPE>
          void init( Crystal::TStructure<T_TYPE> const &_structure,
                     math::iVector3d const &_n,
                     types::t_real _overlap);

        //! Returns size of mesh.
        math::iVector3d mesh() const { return n_; }
        //! Returns size of overlap.
        types::t_real overlap() const { return overlap_; }
        
        //! Iterator to the first divide and conquer box.
        t_Boxes::const_iterator begin() const { return container_.begin(); }
        //! Iterator to the last divide and conquer box.
        t_Boxes::const_iterator end() const { return container_.end(); }
        //! Number of boxes.
        size_t size() const { return container_.size(); }

        //! \brief Guesses parameters of divide and conquer mesh base on the
        //!        desired number of atoms per box.
        template<class T_TYPE>
          math::iVector3d guess_mesh( const Crystal::TStructure<T_TYPE> &_structure, 
                                      size_t _nperbox ) const; 

      protected:
        //! Definition of the mesh.
        math::iVector3d n_;
        //! Overlap between small boxes.
        types::t_real overlap_;
        //! Small cell of a divide and conquer box.
        math::rMatrix3d small_cell_;
        //! Mesh of small-cells.
        t_Boxes container_;
    };

    inline bool operator==(DnCBoxes::Point const &_a, DnCBoxes::Point const &_b)
    {
      return _a.index == _b.index and _a.in_small_box == _b.in_small_box 
             and math::is_zero( (_a.translation - _b.translation).squaredNorm() );
    }


  } // namespace Crystal

} // namespace LaDa

#include "periodic_dnc.impl.h"

#endif

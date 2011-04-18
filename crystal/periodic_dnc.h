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
        typedef std::vector<Point> DnCBox;

        //! Does not initialize to avoid throwing.
        DnCBox() : n_(0,0,0) {}

        //! Creates box.
        template<class T_TYPE>
          void init( Crystal::TStructure<T_TYPE> const &_structure,
                     math::iVector3d const &_n,
                     types::t_real _overlap);

        //! Returns size of mesh.
        math::iVector3d n() const { return n_; }
        //! Returns size of overlap.
        math::iVector3d overlap() const { return overlap_; }



      protected:
        //! Definition of the mesh.
        math::iVector3d n_;
        //! Overlap between small boxes.
        types::t_real overlap_;
        //! Small cell of a divide and conquer box.
        math::rMatrix3d small_cell_;
        //! Type of the container of small-boxes.
        typedef std::vector<DnCBox> container_;
        //! Mesh of small-cells.
        std::vector<DnCBox> container;
    } 
    //! Defines a mesh of DnCBoxes.
    typedef std::vector<DnCBox> DnCBoxes;

    //! \cond 
    class DnCBox;
    //! \endcond

    //! Type of the return of the divide and conquer method.
    struct t_DnCBoxes 
    {
      //! Type of the container returned of the divide and conquer method.
      typedef std::vector< DnCBox > type;
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
    template< class T_TYPE > typename t_DnCBoxes::shared_ptr  
      create_dnc_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                        const math::iVector3d &_n,
                        const types::t_real _overlap_distance );
                           
    class DnCBox
    {
       template< class T_TYPE >
         friend  typename t_DnCBoxes::shared_ptr  
           create_dnc_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                             const math::iVector3d &_n,
                             const types::t_real _overlap_distance );
      public:
        //! Instantiates a divide and conquer box.
        DnCBox() : cell_(), points_() {};
        //! Instantiates a divide and conquer box.
        DnCBox(math::rMatrix3d const &_cell) : cell_(_cell), points_() {};
        //! Copies a DnCBox
        DnCBox( DnCBox const &_c) : cell_(_c.cell_), points_(_c.points_) {};
        //! Description of a point in the box.
        struct Point
        {
          //! index in given structure.
          size_t index; 
          //! Translation in cartesian coordinates.
          math::rVector3d translation;
          //! If true, then position is in small box.
          bool in_small_box;
        };
        //! A vector of points.
        typedef std::vector< Point > t_Points;

        //! Returns the cell of the divide and conquer box.
        const math::rMatrix3d& cell() const { return cell_; }
        //! Returns the index of and atom inside this box, w.r.t. the original structure.
        Point const& operator[](size_t _i) const { return points_[_i]; }
        //! Iterator to the first state in the (large) box.
        t_Points::const_iterator begin() const { return points_.begin(); }
        //! Iterator to the end of the states in the (large) box.
        t_Points::const_iterator end() const { return points_.end(); }
        //! Number of states in the (large) box.
        size_t size() const { return points_.size(); }
        //! checks whether a state is already counted
        bool is_counted(size_t _i) const
        {
          foreach(Point const &point, points_)
            if( point.index == _i and point.in_small_box) return true;
          return false;
        }

        //! Returns iterator to point if found.
        t_Points::const_iterator  find(Point const &_point) const
        {
          t_Points::const_iterator i_first = points_.begin();
          t_Points::const_iterator const i_end = points_.end();
          for(; i_first != i_end; ++i_first)
            if(i_first->index != _point.index)                    continue;
            else if(i_first->in_small_box != _point.in_small_box) continue;
            else if(math::is_zero(i_first->translation-_point->translation))      break;
          return i_first;
        }

      protected:
        //! Cell of this box.
        math::rMatrix3d cell_;
        //! Holds points located in the box.
        t_Points points_;
    };


  } // namespace Crystal

} // namespace LaDa

#include "periodic_dnc.impl.h"

#endif

#ifndef LADA_CRYSTAL_PERIODIC_DNC_H
#define LADA_CRYSTAL_PERIODIC_DNC_H

#include "LaDaConfig.h"

#include <vector>

#include <math/gruber.h>
#include <math/misc.h>

#include "structure.h"


namespace LaDa 
{
  namespace crystal 
  {

    //! Container of divide and conquer boxes.
    class DnCBoxes
    {
      public:
        //! References a point in a periodic divide and conquer box.
        struct Point
        {
          //! Translation in cartesian coordinates.
          math::rVector3d translation;
          //! index in given structure.
          size_t index; 
          //! If true, then position is in small box.
          bool in_small_box;
        };
        //! Defines a divide and conquer box.
        typedef std::vector<Point> t_Box;
        //! Defines a divide and conquer box.
        typedef std::vector<t_Box> t_Boxes;
        //! Type of the iterators over boxes.
        typedef t_Boxes::const_iterator const_iterator;
        //! Type of the constant reference to the boxes.
        typedef t_Boxes::const_reference const_reference;
        //! Value type associated with this container.
        typedef t_Boxes::value_type value_type;

        //! Does not initialize to avoid throwing.
        DnCBoxes() : n_(0,0,0) {}

        //! Creates box.
        template<class T_TYPE>
          void init( TemplateStructure<T_TYPE> const &_structure,
                     math::iVector3d const &_n, types::t_real _overlap);

        //! Returns size of mesh.
        math::iVector3d mesh() const { return n_; }
        //! Returns size of overlap.
        types::t_real overlap() const { return overlap_; }
        
        //! Iterator to the first divide and conquer box.
        const_iterator begin() const { return container_.begin(); }
        //! Iterator to the last divide and conquer box.
        const_iterator end() const { return container_.end(); }
        //! Number of boxes.
        size_t size() const { return container_.size(); }

        //! Returns nth box.
        const_reference operator[](size_t _i) const { return container_[_i]; }

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
      return     _a.index == _b.index 
             and math::is_null( (_a.translation - _b.translation).squaredNorm() );
    }


    //! \brief Guesses parameters of divide and conquer mesh base on the
    //!        desired number of atoms per box.
    math::iVector3d guess_mesh(math::rMatrix3d const &_cell, size_t _N, size_t _nperbox);
    //! \brief Guesses parameters of divide and conquer mesh base on the
    //!        desired number of atoms per box.
    template<class T_TYPE>
      math::iVector3d guess_mesh(TemplateStructure<T_TYPE> const &_structure, size_t _nperbox) 
      { return guess_mesh(_structure.cell(), _structure.size(), _nperbox); }


#   ifdef LADA_INDEX
#     error LADA_INDEX already defined.
#   endif
#   define LADA_INDEX(a, b) (a(0) * b(1) + a(1)) * b(2) + a(2);

    template<class T_TYPE>
      void DnCBoxes::init( const TemplateStructure<T_TYPE> &_structure, 
                           const math::iVector3d &_n,
                           const types::t_real _overlap )
      {
        namespace bt = boost::tuples;
        typedef math::iVector3d iVector3d;
        typedef math::rVector3d rVector3d;
        typedef math::rMatrix3d rMatrix3d;
        typedef typename crystal::TemplateStructure<T_TYPE> :: const_iterator const_iterator;

        // constructs cell of small small box
        math::rMatrix3d const strcell( math::gruber(_structure.cell()) );
        math::rMatrix3d cell(strcell);
        for( size_t i(0); i < 3; ++i ) cell.col(i) /= types::t_real( _n(i) );
        
        // Inverse matrices.
        math::rMatrix3d const invstr(strcell.inverse()); // inverse of structure 
        math::rMatrix3d const invbox(cell.inverse());   // inverse of small box.

        // Edges
        rVector3d const sb_edges
          (
            _overlap / std::sqrt(cell.col(0).squaredNorm()),
            _overlap / std::sqrt(cell.col(1).squaredNorm()),
            _overlap / std::sqrt(cell.col(2).squaredNorm())
          );

        // Constructs mesh of small boxes.
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        container_.clear(); container_.resize(Nboxes);

        // Now adds points for each atom in each box.
        const_iterator i_atom = _structure.begin();
        const_iterator const i_atom_end = _structure.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          // Gets coordinate in mesh of small-boxes. Only works because cell
          // and boxes are commensurate.
          rVector3d const rfrac(invbox*i_atom->pos);
          iVector3d const _ifrac(math::floor_int(rfrac));
          iVector3d const __ifrac(_ifrac(0) % _n(0), _ifrac(1) % _n(1), _ifrac(2) % _n(2));
          iVector3d const ifrac
            (
              __ifrac(0) < 0 ? __ifrac(0) + _n(0): __ifrac(0), 
              __ifrac(1) < 0 ? __ifrac(1) + _n(1): __ifrac(1), 
              __ifrac(2) < 0 ? __ifrac(2) + _n(2): __ifrac(2)
            );
          // Computes index within cell of structure.
          types::t_int const u = LADA_INDEX(ifrac, _n);
#         ifdef LADA_DEBUG
            for(size_t i(0); i < 3; ++i)
              if(ifrac(i) < 0 or ifrac(i) >= _n(i))
                BOOST_THROW_EXCEPTION(error::out_of_range());
            if(u < 0 or u >= Nboxes) BOOST_THROW_EXCEPTION(error::out_of_range());
#         endif

          // creates apropriate point in small-box. 
          DnCBoxes::Point const orig = {cell * (ifrac - _ifrac).cast<math::rMatrix3d::Scalar>(), index, true};
          container_[u].push_back(orig);

          // Finds out which other boxes it is contained in, including periodic images.
          for( types::t_int i(-1 ); i <= 1; ++i )
            for( types::t_int j(-1 ); j <= 1; ++j )
              for( types::t_int k(-1 ); k <= 1; ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;

                // First checks if on edge of small box.
#               ifndef LADA_WITH_EIGEN3
                  rVector3d displaced = rfrac + sb_edges.cwise()*rVector3d(i,j,k);
#               else
                  rVector3d displaced =   rfrac
                                        + (sb_edges.array()*rVector3d(i,j,k).array()).matrix();
#               endif
                iVector3d const boxfrac
                  ( 
                    math::floor_int(displaced(0)) == _ifrac(0) ? ifrac(0): ifrac(0) + i,
                    math::floor_int(displaced(1)) == _ifrac(1) ? ifrac(1): ifrac(1) + j, 
                    math::floor_int(displaced(2)) == _ifrac(2) ? ifrac(2): ifrac(2) + k
                  );
                // Now checks if small box is at edge of periodic structure. 
                iVector3d const strfrac
                  (
                    boxfrac(0) < 0 ? 1: (boxfrac(0) >= _n(0) ? -1: 0),
                    boxfrac(1) < 0 ? 1: (boxfrac(1) >= _n(1) ? -1: 0),
                    boxfrac(2) < 0 ? 1: (boxfrac(2) >= _n(2) ? -1: 0)
                  );
                bool const is_edge(strfrac(0) != 0 or strfrac(1) != 0 or strfrac(2) != 0);

                // Computes index of box where displaced atom is located.
                iVector3d const modboxfrac
                  (
                    boxfrac(0) + strfrac(0) * _n(0),
                    boxfrac(1) + strfrac(1) * _n(1),
                    boxfrac(2) + strfrac(2) * _n(2)
                  );
                types::t_int const uu = LADA_INDEX(modboxfrac, _n);
#               ifdef LADA_DEBUG
                  if(uu < 0 or uu >= container_.size()) BOOST_THROW_EXCEPTION(error::out_of_range());
#               endif

                // Don't need to go any further: not an edge state of either
                // small box or structure.
                if(u == uu  and not is_edge) continue;

                DnCBoxes::Point const overlap 
                  = {
                      orig.translation + strcell*strfrac.cast<math::rMatrix3d::Scalar>(),
                      index,
                      false
                    };
                DnCBoxes::t_Box &box = container_[uu];
                if( box.end() == std::find(box.begin(), box.end(), overlap) ) 
                  box.push_back(overlap);
              }
        }
        n_ = _n;
        overlap_ = _overlap;
      }
#   undef LADA_INDEX

  } // namespace crystal

} // namespace LaDa


#endif

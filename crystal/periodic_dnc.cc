#include "LaDaConfig.h"

#include <limits>
#include <set>

#include <boost/numeric/conversion/converter.hpp>

#include <opt/debug.h>
#include <math/misc.h>

#include "periodic_dnc.h"

namespace LaDa 
{
  namespace Crystal 
  {
      math::iVector3d DnCBoxes::guess_mesh( const Crystal::TStructure< std::string > &_structure, 
                                            size_t _nperbox ) const
      {
        types::t_real const alpha = 0.5;
        types::t_real const c0 = std::sqrt( _structure.cell.col(0).squaredNorm() );
        types::t_real const c1 = std::sqrt( _structure.cell.col(1).squaredNorm() );
        types::t_real const c2 = std::sqrt( _structure.cell.col(2).squaredNorm() );
        types::t_real const Natoms( _structure.atoms.size() );

        int const Nboxes = int(std::floor(types::t_real(Natoms)/types::t_real(_nperbox) + 0.5));
        if(Nboxes == 0) return math::iVector3d::Ones();

        math::iVector3d result(1,1,1);
        types::t_real mini = Natoms / _nperbox * (c0+c1+c2);
        for(size_t n0(1); n0 <= Nboxes; ++n0)
          for(size_t n1(1); n1 <= Nboxes; ++n1)
          {
            size_t n2 = Nboxes / n0 / n1;
            if(n2 == 0) continue;

            types::t_real const a =
                 std::abs(c0/types::t_real(n0) - c1/types::t_real(n1))
               + std::abs(c0/types::t_real(n0) - c2/types::t_real(n2))
               + std::abs(c1/types::t_real(n1) - c2/types::t_real(n2));
            if(a < mini) 
            {
              result(0) = n0;
              result(1) = n1;
              result(2) = n2;
              mini = a;
            }
          }
        return result;
      }

#   ifdef LADA_INDEX
#     error LADA_INDEX already defined.
#   endif
#   define LADA_INDEX(a, b) (a(0) * b(1) + a(1)) * b(2) + a(2);

      void DnCBoxes::init( const Crystal::TStructure< std::string > &_structure, 
                           const math::iVector3d &_n,
                           const types::t_real _overlap )
      {
        namespace bt = boost::tuples;
        typedef math::iVector3d iVector3d;
        typedef math::rVector3d rVector3d;
        typedef math::rMatrix3d rMatrix3d;
        typedef Crystal::TStructure< std::string > :: t_Atoms t_Atoms;


        // constructs cell of small small box
        math::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i ) cell.col(i) *= 1e0 / types::t_real( _n(i) );
        
        // Inverse matrices.
        math::rMatrix3d const invstr( _structure.cell.inverse() ); // inverse of structure 
        math::rMatrix3d const invbox( cell.inverse() );   // inverse of small box.

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
        t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          // Position inside structure cell.
          rVector3d const incell( into_cell(i_atom->pos, _structure.cell, invstr) );
          rVector3d const fracincell(invbox * incell);
          LADA_DOASSERT(math::is_integer(invstr * (incell - i_atom->pos)), "WTF\n");
          // Gets coordinate in mesh of small-boxes.
          iVector3d const __ifrac(math::floor_int(fracincell));
          iVector3d const ifrac
            (
              __ifrac(0) + (__ifrac(0) < 0 ? _n(0): (__ifrac(0) >= _n(0)? -_n(0): 0)), 
              __ifrac(1) + (__ifrac(1) < 0 ? _n(1): (__ifrac(1) >= _n(1)? -_n(1): 0)), 
              __ifrac(2) + (__ifrac(2) < 0 ? _n(2): (__ifrac(2) >= _n(2)? -_n(2): 0))
            );
          // Computes index within cell of structure.
          types::t_int const u = LADA_INDEX(ifrac, _n);
          LADA_ASSERT( u >= 0 and u < Nboxes,
                       "Index out-of-range.\n" << u << " >= " << Nboxes << "\n" );

          // creates apropriate point in small-box. 
          DnCBoxes::Point const orig = {incell - i_atom->pos, index, true};
          container_[u].push_back(orig);

          // Finds out which other boxes it is contained in, including periodic images.
          for( types::t_int i(-1 ); i <= 1; ++i )
            for( types::t_int j(-1 ); j <= 1; ++j )
              for( types::t_int k(-1 ); k <= 1; ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;

                // First checks if on edge of small box.
                rVector3d displaced = fracincell + sb_edges.cwise()*rVector3d(i,j,k);
                iVector3d const __boxfrac( math::floor_int(displaced) );
                iVector3d const boxfrac
                  ( 
                     ifrac(0) + (__boxfrac(0) == ifrac(0) ? 0: i),
                     ifrac(1) + (__boxfrac(1) == ifrac(1) ? 0: j), 
                     ifrac(2) + (__boxfrac(2) == ifrac(2) ? 0: k)
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
                LADA_ASSERT(uu >= 0 and uu < container_.size(), "ERROR");

                // Don't need to go any further: not an edge state of either
                // small box or structure.
                if(u == uu  and not is_edge) continue;

                DnCBoxes::Point const overlap 
                  = {
                      incell + _structure.cell*strfrac.cast<types::t_real>() - i_atom->pos,
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

  } // namespace Crystal

} // namespace LaDa


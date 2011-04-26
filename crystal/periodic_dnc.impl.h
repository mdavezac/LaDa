#include <limits>
#include <set>

#include <opt/debug.h>
#include <math/misc.h>

namespace LaDa 
{
  namespace Crystal 
  {
    template< class T_TYPE >
      math::iVector3d DnCBoxes::guess_mesh( const Crystal::TStructure<T_TYPE> &_structure, 
                                            size_t _nperbox ) const
      {
        const types::t_real c1 = std::sqrt( _structure.cell.col(0).squaredNorm() );
        const types::t_real c2 = std::sqrt( _structure.cell.col(1).squaredNorm() );
        const types::t_real c3 = std::sqrt( _structure.cell.col(2).squaredNorm() );
        const types::t_real Natoms( _structure.atoms.size() );
        const types::t_real Nperbox( _nperbox );

        types::t_real n1, n2, n3;
        n1 =  std::pow( Natoms / Nperbox * c1*c1/c2/c3, 1e0/3e0 );
        n2 =  n1 * c2 / c1;
        n3 =  n1 * c3 / c1;
        if( n1 <= 0.5 ) n1 = 1;
        if( n2 <= 0.5 ) n2 = 1;
        if( n3 <= 0.5 ) n3 = 1;
        return math::iVector3d( rint( n1 ), rint( n2 ), rint( n3 ) );
      }

#   ifdef LADA_INDEX
#     error LADA_INDEX already defined.
#   endif
#   define LADA_INDEX(a, b) (a(0) * b(1) + a(1)) * b(2) + a(2);

    template< class T_TYPE > void
      DnCBoxes::init( const Crystal::TStructure<T_TYPE> &_structure, 
                      const math::iVector3d &_n,
                      const types::t_real _overlap )
      {
        std::cout << "Begin" << std::endl;
        namespace bt = boost::tuples;
        typedef math::iVector3d iVector3d;
        typedef math::rVector3d rVector3d;
        typedef math::rMatrix3d rMatrix3d;
        typedef typename Crystal::TStructure<T_TYPE> :: t_Atoms t_Atoms;
        rVector3d const roundoff( 1e1 * std::numeric_limits<types::t_real>::epsilon()
                                      * rVector3d::Ones() );


        // constructs cell of small small box
        math::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i ) cell.col(i) *= 1e0 / types::t_real( _n(i) );
        std::cout << _n(0) << " " << _n(1) << " " << _n(2) << "\n";
        std::cout << "0\n" << cell << std::endl;
        
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
        rVector3d const lb_edges(invstr*_overlap*rVector3d::Ones() + roundoff);

        // Constructs mesh of small boxes.
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        container_.clear(); container_.resize(Nboxes);

        // Now adds points for each atom in each box.
        typename t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        std::cout << "1" << std::endl;
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          // Position inside structure cell.
          rVector3d const incell( into_cell(i_atom->pos, _structure.cell, invstr) );
          rVector3d const fracincell(invbox * incell);
          LADA_DOASSERT(math::is_integer(invstr * (incell - i_atom->pos)), "WTF\n");
          // Gets coordinate in mesh of small-boxes.
          iVector3d const ifrac
            ( math::floor((fracincell + roundoff).eval()).cast<types::t_int>() ); 
          // Computes index within cell of structure.
          types::t_int const u = LADA_INDEX(ifrac, _n);
          LADA_ASSERT( u >= 0 and u < Nboxes,
                       "Index out-of-range.\n" << u << " >= " << Nboxes << "\n" );

          // creates apropriate point in small-box. 
          DnCBoxes::Point const orig = {incell - i_atom->pos, index, true};
        std::cout << "2 " << u << std::endl;
        std::cout << incell(0) << " " << incell(1) << " " << incell(2) << std::endl;
        std::cout << i_atom->pos(0) << " " << i_atom->pos(1) << " " << i_atom->pos(2) << std::endl;
        std::cout << fracincell(0) << " " << fracincell(1) << " " << fracincell(2) << std::endl;
        std::cout << ifrac(0) << " " << ifrac(1) << " " << ifrac(2) << std::endl;
          container_[u].push_back(orig);
        std::cout << "3" << std::endl;

          // Finds out which other boxes it is contained in, including periodic images.
          for( types::t_int i(-1 ); i <= 1; ++i )
            for( types::t_int j(-1 ); j <= 1; ++j )
              for( types::t_int k(-1 ); k <= 1; ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;

                // First checks if on edge of small box.
                rVector3d displaced = fracincell + sb_edges.cwise()*rVector3d(i,j,k) + roundoff;
                iVector3d const __boxfrac( math::floor(displaced).cast<types::t_int>() );
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
                LADA_DOASSERT(math::is_integer(invstr * overlap.translation), "WTF2\n");
                DnCBoxes::t_Box &box = container_[uu];
        std::cout << "4" << std::endl;
                if( box.end() == std::find(box.begin(), box.end(), overlap) ) 
                  box.push_back(overlap);
        std::cout << "5" << std::endl;
              }
        }
        std::cout << "6" << std::endl;
        n_ = _n;
        overlap_ = _overlap;

      }
#   undef LADA_INDEX

  } // namespace Crystal

} // namespace LaDa


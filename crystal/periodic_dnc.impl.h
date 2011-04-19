#include <limits>
#include <set>

#include <opt/debug.h>

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
        namespace bt = boost::tuples;
        typedef math::iVector3d iVector3d;
        typedef math::rVector3d rVector3d;
        typedef math::rMatrix3d rMatrix3d;
        typedef typename Crystal::TStructure<T_TYPE> :: t_Atoms t_Atoms;
        const types::t_real roundoff = 1e1 * std::numeric_limits<types::t_real>::epsilon();
        const rVector3d rOnes( rVector3d::Ones() );


        // constructs cell of small small box
        math::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i ) cell.col(i) *= 1e0 / types::t_real( _n(i) );

        // Constructs mesh of small boxes.
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        container_.clear(); container_.resize(Nboxes);

        // Now adds points for each atom in each box.
        math::rMatrix3d const inv_small_cell( cell.inverse() );
        math::rMatrix3d const inv_str_cell( _structure.cell.inverse() );
        typename t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          // Position inside structure cell.
          rVector3d const incell( into_cell(i_atom->pos, _structure.cell, inv_str_cell) );
          // Gets coordinate in mesh of small-boxes.
          iVector3d const ifrac( inv_small_cell*incell + roundoff*rOnes ); 
          // Computes index within cell of structure.
          types::t_int const u = LADA_INDEX(ifrac, _n);
          LADA_ASSERT( u < Nboxes, "Index out-of-range.\n" << u << " >= " << Nboxes << "\n" );

          // creates apropriate point in small-box. 
          DnCBoxes::Point const orig = {index, incell - i_atom->pos, true};
          container_[u].push_back(orig);

          // Finds out which other boxes it is contained in, including periodic images.
          for( types::t_int i(-1 ); i <= 1; ++i )
            for( types::t_int j(-1 ); j <= 1; ++j )
              for( types::t_int k(-1 ); k <= 1; ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                rVector3d const displaced(incell + rVector3d(i,j,k) * _overlap);
                // Gets coordinate in mesh of small-boxes. 
                iVector3d const oifrac(inv_small_cell*incell + roundoff*rOnes); 
                // Computes index within cell of structure.
                types::t_int uu = LADA_INDEX(oifrac, _n);

                rVector3d overlaptrans(incell-i_atom->pos);
                // If following test is true, then overlap point still in same small box.
                if( u == uu ) continue;
                // If following test is true, then we are looking at a periodic image.
                else if(uu < 0 or uu >= Nboxes) 
                {
                  rVector3d const period
                  (
                    oifrac(0) < 0 ? 1: (oifrac(0) >= _n(0) ? -1: 1),
                    oifrac(1) < 0 ? 1: (oifrac(1) >= _n(1) ? -1: 1),
                    oifrac(2) < 0 ? 1: (oifrac(2) >= _n(2) ? -1: 1)
                  );
                  overlaptrans += _structure.cell*period;
                  iVector3d const f(oifrac(0) % _n(0), oifrac(1) % _n(1), oifrac(2) % _n(2)); 
                  iVector3d const ff
                  (
                     f(0) < 0 ? f(0) + _n(0): f(0),
                     f(1) < 0 ? f(1) + _n(1): f(1),
                     f(2) < 0 ? f(2) + _n(2): f(2)
                  );
                  uu = LADA_INDEX(ff, _n)
                }
                DnCBoxes::Point const overlap = {index, overlaptrans, false};
                DnCBoxes::t_Box &box = container_[uu];
                if( box.end() != std::find(box.begin(), box.end(), overlap) ) 
                  box.push_back(overlap);
              }
        }
        n_ = _n;
        overlap_ = _overlap;
      }
#   undef LADA_INDEX

  } // namespace Crystal

} // namespace LaDa


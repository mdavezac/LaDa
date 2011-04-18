#include <limits>
#include <set>

#include <opt/debug.h>

namespace LaDa 
{
  namespace Crystal 
  {
#   ifdef LADA_INDEX
#     error LADA_INDEX already defined.
#   endif
#   define LADA_INDEX(a, b) (a(0) * b(1) + a(1)) * b(2) + a(2);

    template< class T_TYPE > typename t_DnCBoxes::shared_ptr  
      periodic_dnc( const Crystal::TStructure<T_TYPE> &_structure, 
                    const math::iVector3d &_n,
                    const types::t_real _overlap_distance )
      {
        namespace bt = boost::tuples;
        typedef math::iVector3d iVector3d;
        typedef math::rVector3d rVector3d;
        typedef math::rMatrix3d rMatrix3d;
        typedef typename Crystal::TStructure<T_TYPE> :: t_Atoms t_Atoms;
        const types::t_real roundoff = 1e1 * std::numeric_limits<types::t_real>::epsilon();


        // constructs cell of small small box
        math::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i ) cell.col(i) *= 1e0 / types::t_real( _n(i) );

        // constructs template small box.
        DnCBox template_box;
        template_box.cell_ = cell;
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        
        // Constructs mesh of small boxes.
        t_DnCBoxes::shared_ptr result( new t_DnCBoxes::type(Nboxes, template_box) );

        // Now adds points for each atom in each box.
        math::rMatrix3d const inv_small_cell( cell.inverse() );
        math::rMatrix3d const inv_str_cell( structure.cell.inverse() );
        typename t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          // Position inside structure cell.
          rVector3d const incell( into_cell(i_atom->pos, _structure.cell, inv_str_cell) );
          // Gets coordinate in mesh of small-boxes. //(0)+roundoff, rfrac(1)+roundoff, rfrac(2)+roudoff);
          iVector3d const ifrac( inv_small_cell * incell + roundoff); 
          // Computes index within cell of structure.
          types::t_int const u( LADA_INDEX(ifrac, _n) );
          LADA_ASSERT( u < Nboxes, "Index out-of-range.\n" << u << " >= " << Nboxes << "\n" );

          // creates apropriate point in small-box. 
          DnCBox::Point const orig = {index, incell - i_atom->pos, true};
          (*result)[u].points_.push_back(orig);

          // Finds out which other boxes it is contained in, including periodic images.
          for( types::t_int i(-extent(0) ); i <= extent(0); ++i )
            for( types::t_int j(-extent(1) ); j <= extent(1); ++j )
              for( types::t_int k(-extent(2) ); k <= extent(2); ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                rVector3d const displaced(incell + rVector3d(i,j,k) * _overlap_distance);
                // Gets coordinate in mesh of small-boxes. //(0)+roundoff, rfrac(1)+roundoff, rfrac(2)+roudoff);
                iVector3d const oifrac(inv_small_cell * incell + roundoff); 
                // Computes index within cell of structure.
                types::t_int uu( LADA_INDEX(oifrac, _n) );

                rVector3d overlaptrans(incell-i_atom->pos);
                // If following test is true, then overlap point still in same small box.
                if( u == uu ) continue
                // If following test is true, then we are looking at a periodic image.
                else if(uu < 0 or uu >= Nboxes) 
                {
                  rVector3d const period
                  (
                    oifrac(0) < 0 ? 1: (oifrac(0) >= _n(0) ? -1: 1),
                    oifrac(1) < 0 ? 1: (oifrac(1) >= _n(1) ? -1: 1),
                    oifrac(2) < 0 ? 1: (oifrac(2) >= _n(2) ? -1: 1)
                  );
                  overlaptrans += structure.cell*period;
                  iVector3d const f(oifrac(0) % _n(0), oifrac(1) % _n(1), oifrac(2) % _n(2)); 
                  iVector3d const ff
                  (
                     f(0) < 0 ? f(0) + _n(0): f(0),
                     f(1) < 0 ? f(1) + _n(1): f(1),
                     f(2) < 0 ? f(2) + _n(2): f(2)
                  );
                  uu = LADA_INDEX(ff, _n)
                }
                DnCBox::Point const overlap = {index, overlaptrans, false};
                DnCBox const &box = (*result)[uu];
                if( box.end() != std::find(box.begin(), box.end(), overlap) ) 
                  box.points_.push_back(overlap);
              }
        }
        return result;
      }
#   undef LADA_INDEX

  } // namespace Crystal

} // namespace LaDa


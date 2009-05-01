//
//  Version: $Id$
//

#include <boost/tuple/tuple_io.hpp>

namespace LaDa 
{
  namespace Crystal 
  {

    template< class T_TYPE > typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                                const atat::iVector3d &_n,
                                const types::t_real _overlap_distance )
      {
        namespace bt = boost::tuples;
        typedef typename Crystal::TStructure<T_TYPE> :: t_Atoms t_Atoms;
        typedef typename t_ConquerBoxes<T_TYPE>::type t_result;
        typedef typename t_ConquerBoxes<T_TYPE>::shared_ptr t_result_ptr;

        // constructs cell of small small box
        atat::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i )
          cell.set_column(i, cell.get_column(i) * ( 1e0 / types::t_real( _n(i) ) ) );
        // constructs cell of large box.
        atat::rMatrix3d large_cell( cell );
        atat::rVector3d odist;
        for( size_t i(0); i < 3; ++i )
        {
          const atat::rVector3d column( cell.get_column(i) );
          const types::t_real a( std::sqrt( atat::norm2( column ) ) );
          odist(i) = _overlap_distance/a;
          large_cell.set_column(i, (1e0 + 2e0*odist(i))  * column );
        }


        // constructs template box.
        ConquerBox<T_TYPE> d_n_c_box;
        bt::get<0>( d_n_c_box.box_ ) = cell;
        bt::get<2>( d_n_c_box.box_ ) = large_cell;
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        t_result_ptr result( new t_result( Nboxes, d_n_c_box ) );
        

        // Adds translation.
        typename t_result::iterator i_box( result->begin() );
        for( size_t i(0); i < _n(0); ++i )
          for( size_t j(0); j < _n(1); ++j )
            for( size_t k(0); k < _n(2); ++k, ++i_box )
            {
              const atat::iVector3d ivec( i, j, k );
              const atat::rVector3d rvec
              ( 
                types::t_real(i) + 0.5,
                types::t_real(j) + 0.5,
                types::t_real(k) + 0.5
              );
              bt::get<1>( i_box->box_ ) = cell * rvec;
            }

        // adds atoms to each box.
        const atat::rMatrix3d inv_str( !_structure.cell );
        const atat::rMatrix3d inv_lcell( !large_cell );
        const atat::rMatrix3d inv_cell( !cell );
        typename t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        std::cout << cell  << "\n\n" << large_cell << "\n\n" << odist << "\n\n";
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          const atat::rVector3d rfrac( inv_str * i_atom->pos );
          const atat::rVector3d ifrac
          (
            rfrac(0) - std::floor( rfrac(0) ),
            rfrac(1) - std::floor( rfrac(1) ),
            rfrac(2) - std::floor( rfrac(2) )
          );
          const atat::rVector3d in_para( _structure.cell * ifrac );
          const atat::rVector3d lfrac( inv_lcell * ( in_para + odist ) );
          const atat::rVector3d sfrac( inv_cell * in_para );
          const types::t_int i( lfrac(0)  );
          const types::t_int j( lfrac(1)  );
          const types::t_int k( lfrac(2)  );
          const types::t_int u( ( i * _n(1) + j ) * _n(2) +  k );
          __ASSERT( u >= Nboxes, "Index out-of-range.\n" )
          const bool is_in_small_box
          (
                i == types::t_int( sfrac(0) )
            and j == types::t_int( sfrac(1) )
            and k == types::t_int( sfrac(2) )
          );
          std::cout << i_atom->pos << "  " 
                    << bt::make_tuple( i,j,k ) << " =? " 
                    << bt::make_tuple( types::t_int( sfrac(0) ), 
                                       types::t_int( sfrac(1) ), 
                                       types::t_int( sfrac(2) ) )
                    << " -> " << ( is_in_small_box ? "true": "false" ) <<  "\n";
          (*result)[u].states_.push_back( bt::make_tuple( index, is_in_small_box ) );
        }
        return result;
      }
   
    template< class T_TYPE >
      atat::iVector3d guess_dnc_params( const Crystal::TStructure<T_TYPE> &_structure, 
                                        size_t _nperbox ) 
      {
        const types::t_real c1 = std::sqrt( atat::norm2( _structure.cell.get_column(0) ) );
        const types::t_real c2 = std::sqrt( atat::norm2( _structure.cell.get_column(1) ) );
        const types::t_real c3 = std::sqrt( atat::norm2( _structure.cell.get_column(2) ) );
        const types::t_real Natoms( _structure.atoms.size() );
        const types::t_real Nperbox( _nperbox );

        types::t_real n1, n2, n3;
        n1 =  std::pow( Natoms / Nperbox * c1*c1/c2/c3, 1e0/3e0 );
        n2 =  n1 * c2 / c1;
        n3 =  n1 * c3 / c1;
        if( n1 <= 0.5 ) n1 == 1;
        if( n2 <= 0.5 ) n2 == 1;
        if( n3 <= 0.5 ) n3 == 1;
        return atat::iVector3d( rint( n1 ), rint( n2 ), rint( n3 ) );
      }

  } // namespace Crystal

} // namespace LaDa


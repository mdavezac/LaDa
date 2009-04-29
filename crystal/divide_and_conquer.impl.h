
namespace LaDa 
{
  namespace Crystal 
  {

    typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                                const atat::iVector3d &_n,
                                const types::t_real _overlap_distance )
      {
        typedef typename t_ConquerBoxes<T_TYPE>::shared_ptr t_result;
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


        { // constructs template box.
          ConquerBox d_n_c_box;
          bt::get<0>( d_n_c_box.box_ ) = cell;
          bt::get<2>( d_n_c_box.box_ ) = large_cell;
          const Nboxes( _n(0) * _n(1) * _n(2) );
          t_result_ptr result( new t_result( Nboxes, d_n_c_box ) );
        }

        // Adds translation.
        t_result::iterator i_box( result->begin() );
        for( size_t i(0); i < _n(0); ++i )
          for( size_t j(0); j < _n(1); ++j )
            for( size_t j(0); k < _n(2); ++k, ++i_box )
            {
              const atat::iVector3d ivec( i, j, k );
              const atat::rVector3d rvec
              ( 
                types::t_real(i) + 0.5,
                types::t_real(j) + 0.5,
                types::t_real(k) + 0.5
              );
              bt::get<0>( i_box->box_ ) = rvec * cell;
            }

        // adds atoms to each box.
        const atat::rMatrix3d inv_str( !_structure.cell );
        const atat::rMatrix3d to_cell_index( (!large_cell) * _structure.cell );
        t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          const atat::rVector3d rfrac( inv_str * ( i_atom->pos + odist ) );
          const atat::rVector3d ifrac
          (
            rfrac(0) - std::floor( rfrac(0) ),
            rfrac(1) - std::floor( rfrac(1) ),
            rfrac(2) - std::floor( rfrac(2) )
          );
          const atat::rVector3d lfrac( inv_lcell * ( i_atom->pos + odist ) );
          const types::t_real i( lfrac(0)  );
          const types::t_real j( lfrac(1)  );
          const types::t_real k( lfrac(2)  );
          const types::t_real u( ( i * _n(1) + j ) * _n(2) +  k );
          __ASSERT( u < Nboxes, "Index out-of-range.\n" )
          const atat::rVector3d sfrac( inv_cell * i_atom->pos );
          const bool is_in_small_box
          (
                i == types::t_int( sfrac(0) )
            and j == types::t_int( sfrac(1) )
            and k == types::t_int( sfrac(2) )
          );
          (*result)[u].states_.push_back
            ( 
               bt::make_tuple( boost::cref( *i_atom ), index, is_in_small_box ) 
             );
        }
      }
   
    typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      atat::iVector3d guess_dnc_params( const Crystal::TStructure<T_TYPE> &_structure, 
                                        size_t _nperbox )
      {
      }

  } // namespace Crystal

} // namespace LaDa



namespace LaDa 
{
  namespace Crystal 
  {

    typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::Structure &_structure, 
                                const atat::iVector3d &_n,
                                const types::t_real _overlap_distance )
      {
        typedef typename t_ConquerBoxes<T_TYPE>::shared_ptr t_result;
        typedef typename t_ConquerBoxes<T_TYPE>::shared_ptr t_result_ptr;
        typename t_ConquerBoxes<T_TYPE>::shared_ptr result( new t_result );

        // constructs cell of small small box
        atat::rMatrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i )
          cell.set_column(i, cell.get_column(i) * ( 1e0 / types::t_real( _n(i) ) ) );
        // constructs cell of large box.
        atat::rMatrix3d large_cell( cell );
        for( size_t i(0); i < 3; ++i )
        {
          const atat::rVector3d column( cell.get_column(i) );
          large_cell.set_column(i,  * ( 1e0 / types::t_real( _n(i) ) ) );

      }
   

  } // namespace Crystal

} // namespace LaDa


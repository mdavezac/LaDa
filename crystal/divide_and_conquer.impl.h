//
//  Version: $Id$
//

#include <boost/tuple/tuple_io.hpp>
#include <set>

namespace LaDa 
{
  namespace Crystal 
  {
    namespace details
    {
      inline size_t box_index( const Eigen::Matrix3d &_inv_str,
                               const Eigen::Matrix3d &_str,
                               const Eigen::Matrix3d &_inv_cell,
                               const Eigen::Vector3d &_pos,
                               const Eigen::Vector3i &_n )
      {
        // Puts atom back into parallelogram.
        const Eigen::Vector3d rfrac( _inv_str * _pos );
        const Eigen::Vector3d ifrac
        (
          rfrac(0) - std::floor( rfrac(0) ),
          rfrac(1) - std::floor( rfrac(1) ),
          rfrac(2) - std::floor( rfrac(2) )
        );
        const Eigen::Vector3d in_para( _str * ifrac );
        // Then finds out which box it belongs to.
        const Eigen::Vector3d frac( _inv_cell * in_para );
        const types::t_int i( frac(0)  );
        const types::t_int j( frac(1)  );
        const types::t_int k( frac(2)  );
        return types::t_int( ( i * _n(1) + j ) * _n(2) +  k );
      }
    }

    template< class T_TYPE > typename t_ConquerBoxes<T_TYPE>::shared_ptr  
      divide_and_conquer_boxes( const Crystal::TStructure<T_TYPE> &_structure, 
                                const Eigen::Vector3i &_n,
                                const types::t_real _overlap_distance )
      {
        namespace bt = boost::tuples;
        typedef typename Crystal::TStructure<T_TYPE> :: t_Atoms t_Atoms;
        typedef typename t_ConquerBoxes<T_TYPE>::type t_result;
        typedef typename t_result::value_type::t_State t_State;
        typedef typename t_ConquerBoxes<T_TYPE>::shared_ptr t_result_ptr;

        // constructs cell of small small box
        Eigen::Matrix3d cell( _structure.cell );
        for( size_t i(0); i < 3; ++i )
          cell.set_column(i, cell.col(i) * ( 1e0 / types::t_real( _n(i) ) ) );
        // constructs cell of large box.
        Eigen::Vector3d odist;
        for( size_t i(0); i < 3; ++i )
        {
          const Eigen::Vector3d column( cell.col(i) );
          const types::t_real a( std::sqrt( column.squaredNorm() ) );
          odist(i) = _overlap_distance/a;
        }

        // constructs template box.
        ConquerBox<T_TYPE> d_n_c_box;
        bt::get<0>( d_n_c_box.box_ ) = cell;
        bt::get<2>( d_n_c_box.box_ ) = odist;
        const size_t Nboxes( _n(0) * _n(1) * _n(2) );
        t_result_ptr result( new t_result( Nboxes, d_n_c_box ) );
        

        // Adds translation.
        typename t_result::iterator i_box( result->begin() );
        for( size_t i(0); i < _n(0); ++i )
          for( size_t j(0); j < _n(1); ++j )
            for( size_t k(0); k < _n(2); ++k, ++i_box )
            {
              const Eigen::Vector3i ivec( i, j, k );
              const Eigen::Vector3d rvec
              ( 
                types::t_real(i) + 0.5,
                types::t_real(j) + 0.5,
                types::t_real(k) + 0.5
              );
              bt::get<1>( i_box->box_ ) = cell * rvec;
            }

        // Periodic images are added to large box only if directions where
        // box is not larger than the cell. Otherwise, when looping over all
        // states, we would go over the same state twice. The following defines
        // the direction for which we can look for periodic images.
        const Eigen::Vector3i extent( _n(0)>1 ? 1:0, _n(1)>1 ? 1:0, _n(2)>1 ? 1:0 );
        // adds atoms to each box.
        const Eigen::Matrix3d inv_str( !_structure.cell );
        const Eigen::Matrix3d inv_cell( !cell );
        typename t_Atoms :: const_iterator i_atom = _structure.atoms.begin();
        typename t_Atoms :: const_iterator i_atom_end = _structure.atoms.end();
        for( size_t index(0); i_atom != i_atom_end; ++i_atom, ++index )
        {
          const types::t_int u
          (
            details::box_index
            (
              inv_str, _structure.cell,
              inv_cell, i_atom->pos, _n 
            ) 
          );
          __ASSERT( u >= Nboxes, "Index out-of-range.\n" )
          (*result)[u].states_.push_back( bt::make_tuple( index, true ) );

          // Finds out which large box it is contained in.
          std::set< size_t > lb_set;
          for( types::t_int i(-extent(0) ); i <= extent(0); ++i )
            for( types::t_int j(-extent(0) ); j <= extent(0); ++j )
              for( types::t_int k(-extent(0) ); k <= extent(0); ++k )
              {
                if( i == 0 and j == 0 and k == 0 ) continue;
                const Eigen::Vector3d displaced
                (
                  i_atom->pos(0) + types::t_real(i) * odist(0),
                  i_atom->pos(1) + types::t_real(j) * odist(1),
                  i_atom->pos(2) + types::t_real(k) * odist(2)
                );
                const types::t_int uu
                (
                  details::box_index
                  (
                    inv_str, _structure.cell,
                    inv_cell, displaced, _n 
                  ) 
                );
                // This condition stops a atom from being added simultateously
                // a box and the large box in which the box is contained.
                if( u != uu ) lb_set.insert( uu );
              }
          // inserts into large boxes.
          std::set< size_t > :: const_iterator i_lb = lb_set.begin();
          const std::set< size_t > :: const_iterator i_lb_end = lb_set.end();
          for(; i_lb != i_lb_end; ++i_lb )
          {
            __DOASSERT( *i_lb >= result->size(), "Index out of range.\n" )
            (*result)[*i_lb].states_.push_back( t_State( index, false ) );
          }
        }
        return result;
      }
   
    template< class T_TYPE >
      Eigen::Vector3i guess_dnc_params( const Crystal::TStructure<T_TYPE> &_structure, 
                                        size_t _nperbox ) 
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
        if( n1 <= 0.5 ) n1 == 1;
        if( n2 <= 0.5 ) n2 == 1;
        if( n3 <= 0.5 ) n3 == 1;
        return Eigen::Vector3i( rint( n1 ), rint( n2 ), rint( n3 ) );
      }

  } // namespace Crystal

} // namespace LaDa


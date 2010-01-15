//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/lambda/bind.hpp>

namespace LaDa
{
  namespace CE
  {
#   if defined( COLHEAD ) || defined(INCOLLAPSE) || defined(INCOLLAPSE2)
#     error "Macros with same names."
#   endif
#   define COLHEAD \
      MixedApproach<T_TRAITS> 
#    define INCOLLAPSE( var ) \
       template< class T_TRAITS >  var COLHEAD
#    define INCOLLAPSE2( code1, code2 ) \
       template< class T_TRAITS >  code1, code2 COLHEAD
 
      INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
        :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
        {
          __DEBUGTRYBEGIN
          namespace bblas = boost::numeric::ublas;
          // Loop over inequivalent configurations.
          t_Vector X( dof() );
          std::fill( _A.data().begin(), _A.data().end(), 
                     typename t_Matrix::value_type(0) );
          std::fill( _b.data().begin(), _b.data().end(), 
                     typename t_Vector::value_type(0) );
          const bblas::range colrange( 0, collapse().dof() );
          typename t_CEFit :: t_Pis :: const_iterator i_pis;
          const bool doce = cefit().dof() > 0;
          const bool dosep = collapse().dof() > 0;
          if( doce )
          {
            i_pis = cefit().pis.begin();
            __ASSERT( cefit().pis.size() != collapse().mapping().size(), 
                      "Inconsistent sizes.\n" )
          }
          for( size_t i(0); i < collapse().mapping().size(); ++i )
          {
            // allows leave-one-out, or leave-many-out.
            if( collapse().mapping().do_skip(i) ) 
             { if( doce ) ++i_pis; continue; }
        
            // create the X vector.
            bblas::vector_range< t_Vector > colX( X, colrange );
            std::fill( colX.begin(), colX.end(), typename t_Vector::value_type(0) );
            if( dosep ) collapse().create_X( i, dim, colX );
            if( doce ) std::copy( i_pis->begin(), i_pis->end(),
                                  X.begin() + collapse().dof() );
            
        
            _A += collapse().mapping().weight(i) * bblas::outer_prod( X, X ); 
            _b += collapse().mapping().weight(i) * collapse().mapping().target(i) * X;
            if( doce ) ++i_pis;
          }
          __DEBUGTRYEND(, "Error in MixedApproach::create_A_n_b()\n" )
        }
 
      INCOLLAPSE2( template< class T_MATRIX, class T_VECTOR > void )
        :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
        {
          __DEBUGTRYBEGIN
          namespace bblas = boost::numeric::ublas;
          dim = _dim;
          create_A_n_b( _A, _b );
 
          const bool dosep( collapse().dof() > 0 );
          if( collapse().dof() )
          {
            const bblas::range seprange( 0, collapse().dof() );
            bblas::matrix_range<T_MATRIX> sepmat( _A, seprange, seprange );
            bblas::vector_range<T_VECTOR> sepvec( _b, seprange );
            collapse().regularization()( collapse().separables(), _A, _b, _dim); 
          }
          if( cefit().dof() )
          {
            const bblas::range cerange( collapse().dof(), dof() );
            bblas::matrix_range<T_MATRIX> cemat( _A, cerange, cerange );
            bblas::vector_range<T_VECTOR> cevec( _b, cerange );
            cefit(). other_A_n_b( cemat, cevec );
 
            if( collapse().dof() )
            {
              bblas::matrix_column< t_Matrix > columnd( coefficients(), _dim );
              const bblas::matrix_column< const t_Matrix > column0( coefficients(), 0 );
              std::copy( column0.begin() + collapse().dof(), column0.end(), 
                         columnd.begin() + collapse().dof() );
            }
          }
          __DEBUGTRYEND(, "Error in MixedApproach::operator()()\n" )
        }
      
      INCOLLAPSE( void ) :: randomize( typename t_Vector :: value_type _howrandom )
      {
        namespace bblas = boost::numeric::ublas;
        collapse().randomize( _howrandom );
        if( cefit().dof() )
        {
          typedef typename t_Matrix :: value_type t_Type;
          bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
          typename bblas::matrix_column< t_Matrix > :: iterator i_c = column0.begin();
          typename bblas::matrix_column< t_Matrix > :: iterator i_c_end = column0.end();
          for( i_c += collapse().dof(); i_c != i_c_end; ++i_c )
            *i_c = t_Type( opt::math::rng() - 5e-1 ) * _howrandom;
          for( size_t i(1); i < coefficients().size2(); ++i )
          {
            bblas::matrix_column< t_Matrix > columnd( coefficients(), i );
            std::copy( column0.begin() + collapse().dof(), column0.end(), 
                       columnd.begin() + collapse().dof() );
          }
        }
      }
 
 
      INCOLLAPSE( void ) :: update_all()
      {
        __DEBUGTRYBEGIN
        if( cefit().dof() )
        {
          namespace bblas = boost::numeric::ublas;
          const bblas::matrix_column< const t_Matrix > column0( coefficients(), 0 );
          for( size_t i(1); i < coefficients().size2(); ++i )
          {
            bblas::matrix_column< t_Matrix > columnd( coefficients(), i );
            std::copy( column0.begin() + collapse().dof(), column0.end(), 
                       columnd.begin() + collapse().dof() );
          }
        }
        if( collapse().dof() ) collapse().update_all();
        __DEBUGTRYEND(, "Error in MixedApproach::update_all()\n" )
      }
      INCOLLAPSE( void ) :: update( types::t_unsigned _d )
      {
        __DEBUGTRYBEGIN
        if( cefit().dof() )
        {
          namespace bblas = boost::numeric::ublas;
          const bblas::matrix_column< t_Matrix > columnd( coefficients(), dim );
          bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
          std::copy( columnd.begin() + collapse().dof(), columnd.end(), 
                     column0.begin() + collapse().dof() );
        }
        if( collapse().dof() ) collapse().update( _d );
        __DEBUGTRYEND(, "Error in MixedApproach::update()\n" )
      }
 
      INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const 
      {
        __DEBUGTRYBEGIN
        __ASSERT(    coefficients().size1() != dof() 
                  or coefficients().size2() != dimensions(),
                  "Inconsistent sizes.\n" )
        opt::ErrorTuple error;
        for(size_t n(0); n < collapse().mapping().size(); ++n )
          if( not collapse().mapping().do_skip(n) )
            error += opt::ErrorTuple(   collapse().mapping().target( n )
                                      - evaluate(n),
                                      collapse().mapping().weight( n ) );
        return error;
        __DEBUGTRYEND(, "Error in MixedApproach::evaluate()\n" )
      }
 
      INCOLLAPSE( typename COLHEAD::t_Matrix::value_type ) :: evaluate(size_t _n) const 
      {
        namespace bblas = boost::numeric::ublas;
        __DEBUGTRYBEGIN
        __ASSERT(    coefficients().size1() != dof() 
                  or coefficients().size2() != dimensions(),
                  "Inconsistent sizes.\n" )
        typename t_Matrix :: value_type intermed(0);
        // Adds Separable part.
        if( collapse().dof() ) intermed += collapse().evaluate(_n);
        // Adds CE part.
        if( cefit().dof() )
        {
          typedef const bblas::matrix_column<const t_Matrix>  t_Column;
          typedef const bblas::vector_range< t_Column > t_Ecis;
          const bblas::matrix_column< const t_Matrix> col( coefficients(), 0 );
          t_Ecis ecis( col, bblas::range( collapse().dof(), dof() ) );
                                                                              
          intermed += ( bblas::inner_prod( ecis, cefit().pis[ _n ] ) );
        }
 
        return intermed;
        __DEBUGTRYEND(, "Error in MixedApproach::evaluate()\n" )
      }
 
      INCOLLAPSE( void ) ::  init( size_t _ranks, size_t _dims )
      {
        __DEBUGTRYBEGIN
        cefit().init( clusters() );
        if ( not ( _ranks and _dims ) )
        { 
          _ranks = 0;
          _dims = 1;
        }
        namespace bblas = boost::numeric::ublas;
        separables().norms.resize( _ranks );
        coefficients().resize(   clusters().size() 
                               + _ranks * t_Separables::t_Mapping::D, _dims );
        const bblas::range rows( 0, _ranks * t_Separables::t_Mapping::D );
        const bblas::range columns( 0, _dims );
        separables().coefficients_interface().set( coefficients(), rows, columns );
        __DEBUGTRYEND(, "Error in MixedApproach::init()\n" )
      }
 
      INCOLLAPSE( void ) ::  reassign()
      {
        if( not cefit().dof() ) return;
        namespace bl = boost::lambda;
        __DEBUGTRYBEGIN
        __ASSERT( coefficients().size1() - collapse().dof() != clusters().size(),
                  "Inconsistent sizes.\n")
 
        typename t_Matrix :: const_iterator1 i_eci =   coefficients().begin1()
                                                     + collapse().dof();
        foreach( t_Clusters::value_type & _clusters, clusters() )
        {
          foreach( ::LaDa::CE::Cluster & _cluster, _clusters ) _cluster.eci = *i_eci;
          ++i_eci;
        }
        __DEBUGTRYEND(,"Error in MixedApproach::reassign().\n" )
      }
 
      INCOLLAPSE( template< class TT_TRAITS > void )
        :: operator=( const MixedApproach< TT_TRAITS >& _c )
        {
          clusters() = _c.clusters();
          coefficients() = _c.coefficients();
          collapse() = _c.collapse();
          cefit() = _c.cefit();
          init( _c.separables().ranks(), _c.separables().dimensions() );
        }
 
#    undef COLHEAD
#    undef INCOLLAPSE
#    undef INCOLLAPSE2

    template< class T_TRAITS >
    std::ostream& operator<<( std::ostream& _stream,
                              const MixedApproach<T_TRAITS> &_col )
    {
      _stream << "Mixed-Approach function: " << _col.separables() << "\n";
      typedef typename MixedApproach<T_TRAITS> :: t_Clusters :: const_iterator t_cit;
   
      t_cit i_clusters = _col.clusters().begin();
      t_cit i_clusters_end = _col.clusters().end();
      for(; i_clusters != i_clusters_end; ++i_clusters )
        _stream << i_clusters->front();
      return _stream;
    }

    template< class T_POSITIONS, class T_CLUSTERS >
      void remove_contained_clusters( const T_POSITIONS &_positions, T_CLUSTERS &_clusters )
      {
        typename T_CLUSTERS :: iterator i_clusters = _clusters.begin();
        typename T_CLUSTERS :: iterator i_clusters_end = _clusters.end();
        for(; i_clusters != i_clusters_end; ++i_clusters )
        {
          if( i_clusters->front().size() < 2 ) continue;
          typename T_CLUSTERS::value_type::iterator i_cluster =  i_clusters->begin();
          typename T_CLUSTERS::value_type::iterator i_cluster_end =  i_clusters->end();
          for(; i_cluster != i_cluster_end; ++i_cluster )
          {
            size_t n(0);
            for(; n < size_t(i_cluster->size()); ++n ) 
            {
              const Eigen::Vector3d &pos( (*i_cluster)[n] );
              if( math::is_zero( atat::norm2( pos ) ) ) continue;
              typename T_POSITIONS :: const_iterator i_pos = _positions.begin();
              typename T_POSITIONS :: const_iterator i_pos_end = _positions.end();
              for(; i_pos != i_pos_end; ++i_pos )
                if( math::is_zero( atat::norm2( *i_pos - pos ) ) ) break;
              if( i_pos == i_pos_end ) break;
            }
            if( n == size_t( i_cluster->size() ) ) break;
          }
          if( i_cluster == i_cluster_end ) continue;
          typename T_CLUSTERS :: iterator dummy = i_clusters;
          --i_clusters;
          _clusters.erase( dummy );
          i_clusters_end = _clusters.end();
        }
      }
  } // end of CE namespace.
} // namespace LaDa

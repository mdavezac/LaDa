//
//  Version: $Id$
//
#include<boost/numeric/ublas/vector_proxy.hpp>
#include<boost/numeric/ublas/matrix_proxy.hpp>
#include<boost/numeric/ublas/operation.hpp>
#include<boost/lambda/bind.hpp>

namespace CE
{
# if defined( COLHEAD ) || defined(INCOLLAPSE) || defined(INCOLLAPSE2)
#   error "Macros with same names."
# endif
# define COLHEAD \
    Mixed<T_COLBASE, T_CEBASE> 
#  define INCOLLAPSE( var ) \
     template< class T_COLBASE, class T_CEBASE >  var COLHEAD
#  define INCOLLAPSE2( code1, code2 ) \
     template< class COLBASE, class T_CEBASE >  code1, code2 COLHEAD

    INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
      :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
      {
        namespace bblas = boost::numeric::ublas;
        __ASSERT( not separables(), "Function pointer not set.\n" )
        __ASSERT( dim >= configurations_.size1(), "Inconsistent sizes.\n" )
        // Loop over inequivalent configurations.
        t_Vector X( dof() );
        std::fill( _A.data().begin(), _A.data().end(), 
                   typename t_Matrix::value_type(0) );
        std::fill( _b.data().begin(), _b.data().end(), 
                   typename t_Vector::value_type(0) );
        const bblas::range colrange( 0, t_ColBase::dof() );
        t_Pis :: const_iterator i_pis = pis.begin();
        for( size_t i(0); i < mapping.size(); ++i, ++i_pis )
        {
          // allows leave-one-out, or leave-many-out.
          if( mapping.do_skip(i) ) continue;
      
          // create the X vector.
          bblas::vector_range< t_Vector > colX( X, colrange );
          std::fill( colX.begin(), colX.end(), typename t_Vector::value_type(0) );
          t_ColBase::create_X( i, colX );
          std::copy( i_pis.begin(), i_pis.end(), X.begin() + t_ColBase::dof() );
          
      
          _A += mapping.weight(i) * bblas::outer_prod( X, X ); 
          _b += mapping.weight(i) * mapping.target(i) * X;
        }
      }

     INCOLLAPSE2( template< class T_MATRIX, class T_VECTOR > void )
       :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
       {
         dim = _dim;
         create_A_n_b( _A, _b );

         const bblas::range separange( 0, t_ColBase::dof() );
         bblas::matrix_range<T_MATRIX> sepmat( _A, seprange, seprange );
         bblas::matrix_range<T_VECTOR> sepvec( _b, seprange );
         regularization( _A, _b, _dim); 
         
         const bblas::range cerange( t_ColBase::dof(), t_FitBase::dof() );
         bblas::matrix_range<T_MATRIX> cemat( _A, cerange, cerange );
         bblas::matrix_range<T_VECTOR> cevec( _b, cerange );
         t_CEBase :: other_A_n_b( cemat, cevec );

         bblas::matrix_column< t_Matrix > columnd( coefficients(), _dim );
         const bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
         std::copy( column0.begin() + t_ColBase::dof(), column0.end(), 
                    columnd.begin() + t_ColBase::dof() );
       }
     
     INCOLLAPSE( void ) :: update_all()
     {
       bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
       for( size_t i(1); i < coefficients().size1(); ++i )
       {
         const bblas::matrix_column< t_Matrix > columnd( coefficients(), _dim );
         std::copy( column0.begin() + t_ColBase::dof(), column0.end(), 
                    columnd.begin() + t_ColBase::dof() );
       }
       t_ColBase::update_all();
     }
     INCOLLAPSE( void ) :: update( size_t _d )
     {
       const bblas::matrix_column< t_Matrix > columnd( coefficients(), _dim );
       bblas::matrix_column< t_Matrix > column0( coefficients(), 0 );
       std::copy( columnd.begin() + t_ColBase::dof(), columnd.end(), 
                  column0.begin() + t_ColBase::dof() );
       t_ColBase::update( _d );
     }

     INCOLLAPSE( opt::ErrorTuple ) :: evaluate()
     {
       namespace bblas = boost::numeric::ublas;
       opt::ErrorTuple error;
       t_Ecis eci( bblas::matrix_column< t_Matrix>( coefficients(), 0 ), 
                   bblas::range( t_ColBase::dof(), dof() ) );
       for(size_t n(0); n < t_ColBase::mapping.size(); ++n )
       {
         if( t_ColBase::mapping.do_skip(n) ) continue;
         bblas::range range( t_ColBase::mapping.range(n) );
         types::t_real intermed(0);
         // Adds Separable part.
         for( bblas::range::const_iterator j( range.begin() ); j != range.end(); ++j )
         {
           const bblas::matrix_column<const t_iMatrix> config( configurations(), *j );
           intermed +=   separables()( config )
                       * t_ColBase::mapping.eweight(n,*j - range.start() );
         }
         // Adds CE part.
         typedef const bblas::vector_range< bblas::matrix_column<t_Matrix> > t_Ecis;
         predic += ( bblas::inner_prod( ecis, t_CEBase::pis[_n] ) );

         error += opt::ErrorTuple( t_ColBase::mapping.target(n) - intermed,
                                   t_ColBase::mapping.weight(n) );
       }
       return error;
     }

     INCOLLAPSE( void ) ::  init()
     {
       coefficients_.resize( dof(), dimensions() );
       const bblas::range rows( 0, t_ColBase.dof() );
       const bblas::range column( 0, t_ColBase.dimensions() );
       t_ColBase::coefficients_interface().set( coefficients(), rows, columns );
       t_ColBase :: init( separables_ );
     }

#  undef COLHEAD
#  undef INCOLLAPSE
#  undef INCOLLAPSE2
} // end of CE namespace.


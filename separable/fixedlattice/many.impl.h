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
    Many<T_TRAITS> 
#  define INCOLLAPSE( var ) \
     template< class T_TRAITS >  var COLHEAD
#  define INCOLLAPSE2( code1, code2 ) \
     template< class T_TRAITS >  code1, code2 COLHEAD

    INCOLLAPSE2( size_t ) :: dimensions() const 
    {
      return std::accumulate( collapses_.begin(), collapses_.end(),
                              boost::bind( &t_Collapse::dimensions, _1 ) );
    } 
    INCOLLAPSE2( size_t ) :: dof() const 
    {
      return std::accumulate( collapses_.begin(), collapses_.end(),
                              boost::bind( &t_Collapse::dof, _1 ) );
    }
    INCOLLAPSE2( size_t ) :: current_dof() const 
    {
      size_t result(0);
      foreach( const t_Collapse &collapse, separables )
        if( collapse.dimensions() < dim ) result += collapse.dof();
    }

    INCOLLAPSE2(template< class T_MATRIX, class T_VECTOR > void) 
      :: create_A_n_b( T_MATRIX &_A, T_VECTOR &_b )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        __ASSERT( not t_ColBase::collapses_, "Function pointer not set.\n" )
        // Loop over inequivalent configurations.
        t_Vector X( current_dof() );
        std::fill( _A.data().begin(), _A.data().end(), 
                   typename t_Matrix::value_type(0) );
        std::fill( _b.data().begin(), _b.data().end(), 
                   typename t_Vector::value_type(0) );
        const t_Collapse &first_collapse = collapses_.front();
        for( size_t i(0); i < t_ColBase::mapping().size(); ++i )
        {
          // allows leave-one-out, or leave-many-out.
          if( t_ColBase::mapping().do_skip(i) )  continue;

          size_t rangestart(0);
          foreach( t_Collapse &collapse, collapses_ )
          {
            if( collapse.dimensions() < dim ) continue;
            const bblas::range colrange( rangestart, rangestart + collapse.dof() );
      
            // create the X vector.
            bblas::vector_range< t_Vector > colX( X, colrange );
            std::fill( colX.begin(), colX.end(), typename t_Vector::value_type(0) );
            collapse.create_X( i, colX );
            rangestart += collapse.dof();
          }
          
      
          _A += first_collapse.weight(i) * bblas::outer_prod( X, X ); 
          _b += first_collapse.weight(i) * first_collapse.mapping().target(i) * X;
        }
        __DEBUGTRYEND(, "Error in Many::create_A_n_b()\n" )
      }

    INCOLLAPSE2( template< class T_MATRIX, class T_VECTOR > void )
      :: operator()( T_MATRIX &_A, T_VECTOR &_b, types::t_unsigned _dim )
      {
        __DEBUGTRYBEGIN
        namespace bblas = boost::numeric::ublas;
        dim = _dim;
        foreach( t_Collapse &collapse, collapses_ ) collapse.dim = _dim;
        create_A_n_b( _A, _b );

        size_t rangestart(0);
        foreach( t_Collapse &collapse, collapses_ )
        {
          const bblas::range seprange( rangestart, rangestart + collapse.dof() );
          bblas::matrix_range<T_MATRIX> sepmat( _A, seprange, seprange );
          bblas::vector_range<T_VECTOR> sepvec( _b, seprange );
          collapse.regularization()( _A, _b, _dim); 
          rangestart += collapse.dof();
        }
        __DEBUGTRYEND(, "Error in Many::operator()()\n" )
      }
    
    INCOLLAPSE( void ) :: randomize( typename t_Vector :: value_type _howrandom )
    {
      foreach( t_Collapse &collapse, collapses_ )
        collapse.randomize( _howrandom );
    }


    INCOLLAPSE( void ) :: update_all()
    {
      __DEBUGTRYBEGIN
      foreach( t_Collapse &collapse, collapses_ ) collapse.update_all();
      __DEBUGTRYEND(, "Error in Many::update_all()\n" )
    }
    INCOLLAPSE( void ) :: update( types::t_unsigned _d )
    {
      __DEBUGTRYBEGIN
      foreach( t_Collapse &collapse, collapses_ ) collapse.update(_d);
      __DEBUGTRYEND(, "Error in Many::update()\n" )
    }

    INCOLLAPSE( opt::ErrorTuple ) :: evaluate() const 
    {
      __DEBUGTRYBEGIN
      __ASSERT(    coefficients().size1() != dof() 
                or coefficients().size2() != dimensions(),
                "Inconsistent sizes.\n" )
      opt::ErrorTuple error;
      const t_Collapse &first_collapse = collapses_.front();
      for(size_t n(0); n < t_ColBase::mapping().size(); ++n )
        if( not first_collapse.mapping().do_skip(n) )
          error += opt::ErrorTuple(   t_ColBase::mapping().target( n )
                                    - evaluate(n),
                                    t_ColBase::mapping().weight( n ) );
      return error;
      __DEBUGTRYEND(, "Error in Many::evaluate()\n" )
    }

    INCOLLAPSE( typename COLHEAD::t_Matrix::value_type ) :: evaluate(size_t _n) const 
    {
      __DEBUGTRYBEGIN
      typename t_Matrix :: value_type intermed(0);
      foreach( t_Collapse &collapse, collapses_ )
        intermed += collapse.evaluate(_n);

      return intermed;
      __DEBUGTRYEND(, "Error in Many::evaluate( size_t )\n" )
    }

    INCOLLAPSE( void ) ::  init( size_t _ranks, size_t _dims )
    {
      __DEBUGTRYBEGIN
      foreach( typename t_Traits :: t_Separables &separable, separables_ )
        separable.init( _ranks, _dims );
      __DEBUGTRYEND(, "Error in Many::init()\n" )
    }


#  undef COLHEAD
#  undef INCOLLAPSE
#  undef INCOLLAPSE2

  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream,
                            const MixedApproach<T_TRAITS> &_col )
  {
    foreach( typename t_Traits :: t_Separables &separable, separables_ )
      _stream << separable << "\n";
    return _stream;
  }

} // end of CE namespace.


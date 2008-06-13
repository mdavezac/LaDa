//
//  Version: $Id$
//

#include "gsl_lsq.h"

#include  <boost/type_traits/is_same.hpp>

namespace Fitting
{
  bool Gsl :: Load( const TiXmlElement &_node )
  {
    TiXmlElement *parent = &_node;
    for(; parent; parent = parent->NextSiblingElement("Fit") )
    {
      if( not parent->Attribute("type") ) continue;
      std::string name = parent->Attribute("type");
      if( name.compare( "Gsl" ) ) break;
    }
    if( not parent )
    {
      parent = _node.FirstChildElement("Fit" );
      if( not parent ) return false;
      return Load( *parent );
    }
    if( parent->Attribute( "tolerance" ) )
      parent->Attribute( "tolerance", &tolerance );
    if( parent->Attribute( "options" ) )
    {
      std::string options = parent->Attribute( "options" );
      doweights = options.find( "weight" ) != npos;
      dosvd = options.find( "svd" ) != npos;
    }
    return true;
  }

  types::t_real Gsl::operator()( t_Vector &_solution )
  {
    size_t dim2 = _solution.size();
    size_t dim1 = b.size();
    __ASSERT( dim1 * dim2 != A.size(), "Incoherent matrix/vector size.\n" )
    gsl_multifit_linear_workspace *work;
    work = gsl_multifit_linear_alloc ( dim1, dim2 );

    GslVector gslb( b );
    GslVector gslx( _solution );
    GslMatrix gslA( dim1, A );
    gsl_matrix *cov = gsl_matrix_alloc( dim2, dim2 );

    size_t rank;
    double chisq;
    if( doweights )
    {
      GslVector gslw( w );
      if( dosvd )
        gsl_multifit_wlinear_svd( gslA(), gslw(), gslb(), 
                                  (double) tolerance, &rank,
                                  gslx(), cov, &chisq, work );
      else
        gsl_multifit_wlinear( gslA(), gslw(), gslb(), 
                              gslx(), cov, &chisq, work );
    }
    else
    {
      if( dosvd )
        gsl_multifit_wlinear_svd( gslA(), gslb(), (double) tolerance, &rank,
                                  gslx(), cov, &chisq, work );
      else
        gsl_multifit_wlinear( gslA(), gslb(), gslx(), cov, &chisq, work );
    }

    gsl_multifit_linear_free( work );
    gsl_matrix_free( cov );
    return chisq;
  }

  namespace details
  {
    GslVector::GslVector( Gsl::t_Vector &_vec )
    {
      __ASSERT( boost::is_same< Gsl::t_Vector::value_type, double>::value,
                "Types are not equivalent, double Gsl::t_Vector::value_type.\n" )
      vector.size = _vec.size();
      vector.stride = 1;
      vector.data = &_vec[0];
      vector.block = &block;
      vector.owner = 0;
      bloc.size = vector.size;
      bloc.data = vector.data;
    }
    GslMatrix::GslMatrix( types::t_int _nrow, Gsl::t_Matrix &_mat )
    {
      __ASSERT( boost::is_same< Gsl::t_Matrix::value_type, double>::value,
                "Types are not equivalent, double Gsl::t_Vector::value_type.\n" )
      bloc.size = _mat.size();
      bloc.data = &_mat[0];
      matrix.size1 = _leaddim;
      matrix.size2 = _mat.size() / leaddim;
      __DOASSERT( matrix.size1 * matrix.size2 == _mat,
                  "Matrix dimensions are not coherent.\n" )
      matrix.tda = 0;
      matrix.data = &_vec[0];
      matrix.block = &block;
      matrix.owner = 0;
    }
  }
}

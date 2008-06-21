//
//  Version: $Id$
//

#include "gsl_lsq.h"

#include<boost/type_traits/is_same.hpp>

namespace Fitting
{
  bool Gsl :: Load( const TiXmlElement &_node )
  {
    const TiXmlElement *parent = &_node;
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
      doweights = options.find( "weight" ) != std::string::npos;
      dosvd = options.find( "svd" ) != std::string::npos;
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

    details::GslVector gslb( b );
    details::GslVector gslx( _solution );
    details::GslMatrix gslA( dim1, A );
    gsl_matrix *cov = gsl_matrix_alloc( dim2, dim2 );

    size_t rank;
    double chisq;
    if( doweights )
    {
      details::GslVector gslw( w );
      if( dosvd )
        gsl_multifit_wlinear_svd( (gsl_matrix*)gslA, (gsl_vector*)gslw,
                                  (gsl_vector*)gslb, (double) tolerance, &rank,
                                  (gsl_vector*)gslx, cov, &chisq, work );
      else
        gsl_multifit_wlinear( (gsl_matrix*)gslA, (gsl_vector*)gslw,
                              (gsl_vector*)gslb, (gsl_vector*)gslx, 
                              cov, &chisq, work );
    }
    else
    {
      if( dosvd )
        gsl_multifit_linear_svd( (gsl_matrix*)gslA, (gsl_vector*)gslb,
                                 (double) tolerance, &rank,
                                 (gsl_vector*)gslx, cov, &chisq, work );
      else
        gsl_multifit_linear( (gsl_matrix*)gslA, (gsl_vector*)gslb, (gsl_vector*)gslx, 
                             cov, &chisq, work );
    }

    gsl_multifit_linear_free( work );
    gsl_matrix_free( cov );
    return chisq;
  }

  namespace details
  {
    GslVector::GslVector( Gsl::t_Vector &_vec )
    {
      if( not boost::is_same< Gsl::t_Vector::value_type, double>::value )
       __THROW_ERROR("Types are not equivalent, double Gsl::t_Vector::value_type.\n" )
      vector.size = _vec.size();
      vector.stride = 1;
      vector.data = &_vec[0];
      vector.block = &bloc;
      vector.owner = 0;
      bloc.size = vector.size;
      bloc.data = vector.data;
    }
    GslMatrix::GslMatrix( types::t_int _nrow, Gsl::t_Matrix &_mat )
    {
      if( not boost::is_same< Gsl::t_Matrix::value_type, double>::value )
       __THROW_ERROR("Types are not equivalent, double Gsl::t_Vector::value_type.\n" )
      bloc.size = _mat.size();
      bloc.data = &_mat[0];
      matrix.size1 = _nrow;
      matrix.size2 = _mat.size() / _nrow;
      __DOASSERT( matrix.size1 * matrix.size2 == _mat.size(),
                  "Matrix dimensions are not coherent.\n" )
      matrix.tda = 0;
      matrix.data = &_mat[0];
      matrix.block = &bloc;
      matrix.owner = 0;
    }
  }
}

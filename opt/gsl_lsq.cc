//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "gsl_lsq.h"
#include "gsl.h"

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
    size_t dim = _solution.size();
    __ASSERT( dim * dim != A.size(),
                "Incoherent matrix/vector size: " 
             << dim << "x" << dim << "!=" << A.size() <<".\n" )
    gsl_multifit_linear_workspace *work;
    work = gsl_multifit_linear_alloc ( dim, dim );

    ::Gsl::Vector gslb( b );
    ::Gsl::Vector gslx( _solution );
    ::Gsl::Matrix gslA( dim, A );
    gsl_matrix *cov = gsl_matrix_alloc( dim, dim );

    size_t rank;
    double chisq = 0;
    if( doweights )
    {
      ::Gsl::Vector gslw( w );
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

}

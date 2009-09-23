//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/python.hpp>

#include <pyublas/numpy.hpp>
#include "../cgs.h"


namespace LaDa
{
  namespace Python
  {
    namespace details
    {
      typedef pyublas::numpy_matrix<types::t_real> t_Matrix; 
      typedef pyublas::numpy_vector<types::t_real> t_Vector; 
//     template< class T_TYPE >
//       void extract_matrix( const boost::python::list &_matrix,
//                            boost::numeric::ublas::matrix< T_TYPE > &_result )
//       {
//         namespace bp = boost::python;
//         __TRYBEGIN
//           __DOASSERT( not bp::len(_matrix), "Matrix object is of size 0.\n" )
//           __DOASSERT( not bp::len(_matrix[0]), "Matrix object is of size ?x0.\n" )
//           const size_t nbrows( bp::len(_matrix) );
//           const size_t nbcols( bp::len(_matrix[0]) );
//           _result.resize( nbrows, nbcols );
//           for( size_t i(0); i < nbrows; ++i )
//           {
//             const bp::object &row = _matrix[i];
//             __DOASSERT( bp::len(_matrix[0]) != nbcols,
//                         "Inconsistent number of columns in matrix." )
//             for( size_t j(0); j < nbcols; ++j )
//               _result(i,j) = bp::extract< types::t_real >( row[j] );
//           }
//         __TRYEND(,"Error while converting list to matrix.\n" )
//       }
//     template< class T_TYPE >
//       void extract_vector( const boost::python::list &_vector,
//                            boost::numeric::ublas::vector< T_TYPE > &_result )
//       {
//         namespace bp = boost::python;
//         __TRYBEGIN
//           const size_t s( bp::len(_vector) );
//           __DOASSERT( not s, "Vector of length 0.\n" )
//           _result.resize( s );
//           for( size_t i(0); i < s; ++i )
//             _result(i) = bp::extract< types::t_real >( _vector[i] );
//         __TRYEND(,"Error while converting list to vector.\n" )
//       }
//
      boost::python::tuple pcgs( t_Matrix const & _matrix,
                                 t_Vector &_x,
                                 t_Vector const &_b,
                                 const types::t_real _tolerance,
                                 const size_t _itermax,
                                 const bool _verbose )
      {
        namespace bp = boost::python;

        // construct cgs object.
        LaDa::Fitting::Cgs cgs;
        cgs.verbose = _verbose;
        cgs.tolerance = _tolerance;
        cgs.itermax = _itermax;

        // Performs cgs.
        LaDa::Fitting::Cgs::t_Return result = cgs( _matrix, _x, _b );

        // constructs output and returns.
        return bp::make_tuple( result.first, result.second );
      }

      boost::python::tuple plsq( t_Matrix const &_matrix, 
                                 t_Vector &_x,
                                 t_Vector const &_b,
                                 const types::t_real _tolerance,
                                 const size_t _itermax,
                                 const bool _verbose )
      {
        namespace bp = boost::python;
        namespace bnu = boost::numeric::ublas;

        // construct cgs object.
        LaDa::Fitting::Cgs cgs;
        cgs.verbose = _verbose;
        cgs.tolerance = _tolerance;
        cgs.itermax = _itermax;

        // Performs cgs.
        LaDa::Fitting::Cgs::t_Return
          result = cgs
                   ( 
                     bnu::prec_prod( bnu::trans( _matrix ), _matrix ),
                     _x,
                     bnu::prec_prod( bnu::trans( _matrix ), _b )
                   );

        // constructs output and returns.
        return bp::make_tuple( result.first, result.second );
      }

    } // namespace details

    void expose_cgs()
    {
      namespace bp = boost::python;
      bp::def
      ( 
        "cgs", 
        &details::pcgs, 
        (
          bp::arg("A"), 
          bp::arg("x"), 
          bp::arg("b"), 
          bp::arg("tolerance") = LaDa::types::tolerance,
          bp::arg("itermax") = 50,
          bp::arg("verbosity") = false
        ),
        "Solves the equation A*x=b, with A an input matrix, b an input vector,\n"
        "and x the starting guess vector.\n The number of rows of A and the size\n"
        "of the vectors x and b must match. A is encoded like a C matrix. Returns\n"
        "a tuple (x, residual, iteration), where x is the solution vector and\n"
        "iteration the number of iterations performed."
      );
    }
    void expose_llsq()
    {
      namespace bp = boost::python;
      bp::def
      ( 
        "linear_lsq", 
        &details::plsq, 
        (
          bp::arg("A"), 
          bp::arg("x"), 
          bp::arg("b"), 
          bp::arg("tolerance") = LaDa::types::tolerance,
          bp::arg("itermax") = 50,
          bp::arg("verbosity") = false
        ),
        "Solves the linear least-square problem, with A a matrix of measurement\n"
        "parameters, x a vector of coefficents, and b a vector of measurements.\n"
        "Returns a tuple (x, residual, iteration), where x is the solution vector\n"
        "and iteration the number of iterations performed. This method uses\n"
        "LaDa.cgs to obtain the solution."
      );
    }
  } // namespace Python
} // namespace LaDa

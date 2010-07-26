#include "LaDaConfig.h"

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
//         LADA_TRY_BEGIN
//           LADA_DO_NASSERT( not bp::len(_matrix), "Matrix object is of size 0.\n" )
//           LADA_DO_NASSERT( not bp::len(_matrix[0]), "Matrix object is of size ?x0.\n" )
//           const size_t nbrows( bp::len(_matrix) );
//           const size_t nbcols( bp::len(_matrix[0]) );
//           _result.resize( nbrows, nbcols );
//           for( size_t i(0); i < nbrows; ++i )
//           {
//             const bp::object &row = _matrix[i];
//             LADA_DO_NASSERT( bp::len(_matrix[0]) != nbcols,
//                         "Inconsistent number of columns in matrix." )
//             for( size_t j(0); j < nbcols; ++j )
//               _result(i,j) = bp::extract< types::t_real >( row[j] );
//           }
//         LADA_TRY_END(,"Error while converting list to matrix.\n" )
//       }
//     template< class T_TYPE >
//       void extract_vector( const boost::python::list &_vector,
//                            boost::numeric::ublas::vector< T_TYPE > &_result )
//       {
//         namespace bp = boost::python;
//         LADA_TRY_BEGIN
//           const size_t s( bp::len(_vector) );
//           LADA_DO_NASSERT( not s, "Vector of length 0.\n" )
//           _result.resize( s );
//           for( size_t i(0); i < s; ++i )
//             _result(i) = bp::extract< types::t_real >( _vector[i] );
//         LADA_TRY_END(,"Error while converting list to vector.\n" )
//       }
//
      boost::python::tuple pcgs( t_Matrix const & _matrix,
                                 t_Vector _x,
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
        return bp::make_tuple( _x, result.first, result.second );
      }

      boost::python::tuple plsq( t_Matrix const &_matrix, 
                                 t_Vector _x,
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
        return bp::make_tuple( _x, result.first, result.second );
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
        "Solves the equation A*x=b, with A an input matrix, b an input vector, "
        "and x the starting guess vector.\n\nThe number of rows of A and the size "
        "of the vectors x and b must match. A is encoded like a C matrix.\nReturns "
        "a tuple (residual, iteration)."
      );
    }
    void expose_llsq()
    {
      namespace bp = boost::python;
      std::ostringstream sst; sst << types::tolerance;
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
        ( "Solves the linear least-square problem\n\n"
          "Returns a tuple (residual, iteration).\n"
          "@param A: matrix of measurement\n"
          "@type A: numpy 2d-array.\n"
          "@param x: output parameters.\n"
          "@type x: numpy 1d-array.\n"
          "@param b: vector of measurements.\n"
          "@type b: numpy 1d-array.\n"
          "@param tolerance: tolerance of the function. Defaults to " + sst.str() + "\n."
          "@type tolerance: float\n" 
          "@param itermax: Maximum number of cgs iterations. Defaults to 50\n."
          "@type itermax: int\n" 
          "@param verbosity: chattiness. Defaults to False.\n"
          "@type itermax: bool\n" ).c_str()
      );
    }
  } // namespace Python
} // namespace LaDa

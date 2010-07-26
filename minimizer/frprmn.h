#ifndef _OPT_FRPRMN_H_
#define _OPT_FRPRMN_H_

#include "LaDaConfig.h"
#include "FCMangle.h"

#ifdef _DEBUG
#include <stdexcept>
#endif

#include <functional>

#include <boost/tuple/tuple.hpp>

#include <tinyxml/tinyxml.h>

#include <opt/debug.h>
#include <opt/types.h>

#include "variant.h"

//! \cond
extern "C" double FC_GLOBAL(frprmn, FRPRMN)
                           ( double *, const int *, const double *,
                             const double *, const double *, int *,
                             double *, void *, const int *,
                             const double * );
//! \endcond

namespace LaDa
{
  namespace Minimizer
  {
    //! \cond
    namespace details
    {
      extern void* frpr_pointer_;                           
      template< class T_DATAPAIR > double call_frpr(double* _x, double* _y);
    }
    //! \endcond

    // \brief Interface to the fortran minimizer 
    class Frpr 
    {
        
      public:
        double tolerance; //!< Convergence for overall algorithm
        double line_tolerance; //!< for linmin
        double zeps; //!< for linmin
        int itermax; //!< maximum number of iterations
        double fret; //!< minimum value of the function
        double rtol; //!< achieved tolerance

        //! Constructor 
        Frpr () : tolerance(types::tolerance),
                  line_tolerance(types::tolerance),
                  zeps(types::tolerance),
                  itermax(500), lock(false) {}
        //! Copy Constructor
        Frpr   (const Frpr& _c)
             : tolerance(_c.tolerance),
               line_tolerance(_c.line_tolerance),
               zeps(_c.zeps),
               itermax(_c.itermax), lock(false) {}

        //! Destructor
        virtual ~Frpr(){};

        //! Minimization functor
        template< class T_FUNCTION >
          typename T_FUNCTION :: t_Container :: value_type
            operator()( const T_FUNCTION &_func,
                        typename T_FUNCTION :: t_Container &_arg ) const
            { return operator_< T_FUNCTION,
                                typename T_FUNCTION::t_Container,
                                typename T_FUNCTION::t_Container :: value_type
                              >( _func, _arg ); }
        //! Minimization functor
        template< class T_FUNCTION >
          typename T_FUNCTION :: t_Return
            operator()( const T_FUNCTION &_func,
                        typename T_FUNCTION :: t_Arg &_arg ) const
            { return operator_< T_FUNCTION,
                                typename T_FUNCTION :: t_Arg,
                                typename T_FUNCTION :: t_Return
                              >( _func, _arg ); }

        //! Load minimizer parameters from XML
        bool Load( const TiXmlElement &_element );

        //! \brief Finds the node - if it is there - which describes this minimizer
        //! \details Looks for a \<Minimizer\> tag first as \a _node, then as a
        //!          child of \a _node. Different minimizer, defined by the
        //!          attribute types are allowed:
        const TiXmlElement* find_node( const TiXmlElement &_node );

        //! \brief Loads Minimizer directly from \a _node.
        //! \details If \a _node is not the correct node, the results are undefined.
        bool Load_( const TiXmlElement &_node );

        //! sets all parameters
        void set_parameters( types::t_unsigned _itermax,
                             types::t_real _tolerance,
                             types::t_real _linetolerance,
                             types::t_real _zeps )
        {
          itermax = _itermax;
          tolerance = _tolerance;
          line_tolerance = _linetolerance;
          zeps = _zeps;
        }

      protected:
        template< class T_FUNCTION, class T_CONTAINER, class T_RETURN >
          T_RETURN operator_( const T_FUNCTION &_func, T_CONTAINER &_arg ) const;

        mutable bool lock;
    };
    LADA_REGISTER_MINIMIZER_VARIANT_HEADER( Frpr, "original VFF" )
   
    namespace details
    {
      template< class T_DATAPAIR > double call_frpr(double* _x, double* _y)
      {
        namespace bt = boost::tuples;
        T_DATAPAIR &_this = *(static_cast<T_DATAPAIR*>(frpr_pointer_));
        std::copy( _x, _x + bt::get<1>(_this).size(), bt::get<1>(_this).begin() );
        bt::get<0>(_this).gradient( bt::get<1>(_this), _y );
        return bt::get<0>(_this)( bt::get<1>(_this) );
      }
    }


    template<typename T_FUNCTION, class T_CONTAINER, class T_RETURN> 
      T_RETURN Frpr :: operator_( const T_FUNCTION& _function, T_CONTAINER& _arg ) const
      {
        LADA_DO_NASSERT( lock == true, "Race condition.\n" )
        lock = true;

        int x_length = _arg.size();
        double x[x_length];
        std::copy( _arg.begin(), _arg.end(), x );
        int iter = 0;
  
#       define _MXPARM_ 10000
        LADA_DO_NASSERT( x_length > _MXPARM_,
                      "Current optimizer cannot go beyond " << _MXPARM_ << " variables\n"
                      "Change file cited above, as well as variable mxparm in "
                      "minimizer/df1dim.f90 and opt/linmin.f90 and recompile "
                      "if you need to optimize larger structures.\n";)
#       undef _MXPARM_
        double result;
  
        typedef boost::tuples::tuple< const T_FUNCTION&, T_CONTAINER& > t_DataPair;
        double( *ptr_func )( double*, double* )
          = &details::template call_frpr<t_DataPair>;
        t_DataPair data_pair(_function, _arg);
        details::frpr_pointer_ = (void*) &data_pair;
        FC_GLOBAL(frprmn, frprmn) ( x, &x_length, &tolerance,
                                    &line_tolerance, &zeps,
                                    &iter, &result, (void*)ptr_func,
                                    &itermax, &rtol );
        std::copy( x, x + _arg.size(), _arg.begin() );

        { // recomputes gradient just to make sure.
          typedef typename T_CONTAINER::value_type t_Type;
          t_Type *grad = new t_Type[ _arg.size() ];
          _function.gradient( _arg, grad );
          delete[] grad;
        }
   
  
        lock = false;
        return _function( _arg );
      }

  }
} // namespace LaDa
#endif

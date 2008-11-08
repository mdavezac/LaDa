//
//  Version: $Id$
//

#ifndef _CE_REGULARIZATION_H_
#define _CE_REGULARIZATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/tuple/tuple.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/cgs.h>
#include <opt/gsl_simplex.h>
#include <opt/errors.h>
#include <crystal/structure.h>

#include "fit.h"

namespace LaDa
{
  namespace CE
  {
    // forward declaration
    //! \cond
    class Regulated;
    //! \endcond

    //! \brief Computes CV scores and reduces number of clusters to zero.
    //! \details Regulated::clusters are unchanged at the end of the run.
    //! \brief Regulated Cluster-Expansion.
    //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
    template< class T_MINIMIZER >
    void drautz_diaz_ortiz( Regulated &_reg,
                            const T_MINIMIZER &_minimizer,
                            types::t_int _verbosity = 0,
                            types::t_real _initweights = 0e0 );


    //! \brief Regulated Cluster-Expansion.
    //! \see <A HREF="http://dx.doi.org/10.1103/PhysRevB.73.224207"> Ralf Drautz
    //!      and Alejandro Diaz-Ortiz, PRB \bf 73, 224207 (2007)</A>.
    class Regulated : public RegulatedFit< CE::FittingPolicy::Excluded<> > 
    {
        //! Type of the base class.
        typedef RegulatedFit< CE::FittingPolicy::Excluded<> > t_Base;
      public:
        //! Type of the fitting matrices.
        typedef BaseFit::t_Matrix t_Matrix;
        //! Type of the fitting target vectors.
        typedef BaseFit::t_Vector t_Vector;
        //! Type of the return.
        typedef types::t_real t_Return;
        //! Type of the argument for the minimizer.
        typedef BaseFit::t_Vector t_Arg;

      public:
        //! The clusters to fit.
        t_Clusters clusters;
        //! The fitting procedure.
        Fitting::Cgs cgs;

        //! Constructor.
        Regulated() {};
        //! Destructor.
        ~Regulated() {};

        //! Evaluates the cv score for the weights on input.
        t_Return operator()( const types::t_real * _arg ) const; 
        //! Single fit.
        opt::ErrorTuple fit( BaseFit::t_Vector &_x,
                             const types::t_real *_weights ) const
          { return t_Base::operator()( _x, _weights, cgs ); }
        //! Evaluates the gradient.
        void gradient( const types::t_real * _arg,
                       types::t_real *_gradient ) const;
        //! Reduce regulated function and cluster by 1.
        Cluster reduce();
        //! Reassigns ecis from argument values.
        void reassign( const BaseFit::t_Vector &_arg )
          { t_Base :: reassign( _arg, clusters ); }
        //! Init on owned clusters.
        void init() { t_Base :: init( clusters ); }
    };



  } // end of namespace CE
} // namespace LaDa

#include "regularization.impl.h"

#endif 

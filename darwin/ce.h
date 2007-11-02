//
//  Version: $Id$
//
#ifndef _CE_H_
#define _CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <algorithm>
#include <functional>
#ifndef __PGI
  #include<ext/algorithm>
  using __gnu_cxx::copy_n;
#endif

#include <tinyxml/tinyxml.h>
#include <eo/eoOp.h>

#include "lamarck/structure.h"
#include "lamarck/functional_builder.h"
#include "opt/opt_function_base.h"
#include "opt/types.h"

#include "evaluator.h"
#include "individual.h"
#include "single_site.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

//! all things Cluster Expansion (should be...)
namespace CE
{

  //! \ingroup Genetic
  //! \brief Interface between the genetic algorithm and the cluster expansion
  //!        functional
  //! \xmlinput The functional is expected to be found directly under the
  //!           overall \<Job\/> tag. The exact format should be described
  //!           elsewhere.
  class Darwin : public VA_CE :: Functional_Builder
  {
    protected:
      //! Type of the base class
      typedef VA_CE::Functional_Builder t_Base;
      //! single-cell-shape specialized cluster expansion functional type
      typedef t_Base::t_VA_Functional t_Functional;
      
    protected:
      //! Single-cell-shape cluster expansion functional instance
      t_Functional functional; 
      //! Structure (cell-shape) for which the CE functional is specialized
      Ising_CE::Structure &structure; 

    public:
      //! Constructor
      Darwin( Ising_CE::Structure& _str ) : functional(), structure(_str) {}; 
      //! Copy Constructor
      Darwin   ( const Darwin &_c )
             : t_Base(_c), functional(_c.functional),
               structure( _c.structure )  {}
      //! Destructor
      ~Darwin();

      //! Load the CE functional from input and creates cell-shape
      //! specialization
      bool Load( const TiXmlElement &_node );

      //! \brief Initializes the functional prior to evaluations.
      //! \details Before calling \e any evaluation function, this member must
      //!          be called. It makes sure that the functional acts upon the
      //!          correct set of variables.
      //! \param _indiv Individual to be evaluated in the near futur.
      //! \note T_INDIVIDUAL the type of the %GA individual. It is expected
      //!                     that this instance of this type can return an
      //!                     "Object" which itself can
      //!                     return the container of variables.
      template< class T_INDIVIDUAL >
      void init( T_INDIVIDUAL &_indiv )
        { functional.set_variables( &_indiv.Object().Container() ); }

      //! Returns the evaluation of the CE functional
      types::t_real evaluate() { return functional.evaluate(); }
      template< class T_QUANTITYGRADIENTS >
      //! \brief Computes the gradient of the CE functional
      //! \param _grad is a container holding the gradient. It should accept a
      //!              begin() member routine for which std::copy() can be called.
      void evaluate_gradient( T_QUANTITYGRADIENTS& _grad );
      //! \brief Computes the gradient and returns the evaluation of the CE functional
      //! \param _grad is a container holding the gradient. It should accept a
      //!              begin() member routine for which std::copy() can be called.
      template< class T_QUANTITYGRADIENTS >
      types::t_real evaluate_with_gradient( T_QUANTITYGRADIENTS& _grad );
      //! Returns the gradient of the CE functional in direction \a _pos
      types::t_real evaluate_one_gradient( types::t_unsigned _pos ) 
        {  return functional.evaluate_one_gradient( _pos ); }
  };



  template< class T_QUANTITYGRADIENTS >
  inline void Darwin::evaluate_gradient( T_QUANTITYGRADIENTS& _grad )
  {
    // note that t_Function sets gradient to 0
    types::t_unsigned N = functional.size();
    types::t_real *gradient = new types::t_real[N];
    types::t_real *keep = gradient;
    if ( not gradient ) throw std::runtime_error( "Could not allocate memory" );
    functional.evaluate_gradient( gradient );
    copy_n( gradient, N, _grad.begin() );
    delete[] keep;
  }
  template< class T_QUANTITYGRADIENTS >
  inline types::t_real Darwin::evaluate_with_gradient( T_QUANTITYGRADIENTS& _grad )
  {
    types::t_unsigned N = functional.size();
    types::t_real *gradient = new types::t_real[N];
    types::t_real *keep = gradient;
    if ( not gradient ) throw std::runtime_error( "Could not allocate memory" );
    types::t_real result = functional.evaluate_with_gradient( gradient );
    copy_n( gradient, N, _grad.begin() );
    delete[] keep;
    return result;
  }

} // namespace CE


#endif // _CE_OBJECT_H_

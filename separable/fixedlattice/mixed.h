//
//  Version: $Id$
//
#ifndef _CE_MIXED_H_
#define _CE_MIXED_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>

namespace Traits 
{
  namespace CE
  {
    template< class T_COLTRAITS, class T_CEBASE >
     class MixedApproach
     {
       public:
                        
       //! Original separable traits.
       typedef typename T_COLTRAITS :: t_Separables :: t_Traits t_OrigSepTraits;
       //! Separable function traits.
       typedef Separables
               < 
                 typename t_OrigSepTraits :: t_Mapping,
                 typename t_OrigSepTraits :: t_Policy,
                 ::CE::Policy::MatrixRangeCoefficients,
                 typename t_OrigSepTraits :: t_Vector
               > 
               t_SepTraits;
       //! Type of the separable function.
       typedef ::CE::Separables< t_SepTraits > t_Separables;
       protected:
       typedef typename T_COLTRAITS :: t_RegPolicy
                                    ::template rebind< t_Separables > :: other t_Reg;
       typedef typename T_COLTRAITS :: t_UpdatePolicy
                           ::template rebind< t_Separables,
                                              typename T_COLTRAITS :: t_Mapping, 
                                              typename T_COLTRAITS :: t_iMatrix >
                           :: other t_Up;
       public:
       //! collapse traits.
       typedef Collapse
               < 
                 t_Separables,
                 typename T_COLTRAITS :: t_Mapping,
                 t_Reg,
                 typename T_COLTRAITS :: t_iMatrix,
                 t_Up
               > 
               t_ColTraits;
       //! Type of the collapse functor base.
       typedef ::CE::Collapse< t_ColTraits > t_Collapse;
       //! CE fit base.
       typedef T_CEBASE t_CEBase;
     };

  } // end of CE namespace 

} // end of Traits namespace
namespace CE
{

  template< class T_TRAITS >
    class MixedApproach : protected T_TRAITS :: t_Collapse,
                          protected T_TRAITS :: t_CEBase
    {
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the separable collapse functor.
        typedef typename t_Traits :: t_Collapse t_ColBase;
        //! Type of the CE fitting base class.
        typedef typename t_Traits :: t_CEBase  t_CEBase;
        //! Type of the separable function.
        typedef typename t_Traits :: t_Separables  t_Separables;
        //! Type of the configuration matricex.
        typedef typename t_Traits :: t_ColTraits :: t_iMatrix t_iMatrix;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_SepTraits :: t_Matrix t_Matrix;
        //! Type of the matrices.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Matrix t_OMatrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_OrigSepTraits :: t_Vector t_OVector;

        //! Constructor.
        MixedApproach() : t_ColBase(), t_CEBase() {}
        //! Destructor.
        ~MixedApproach() {}

        //! Returns a reference to the CE fitting base.
        t_CEBase& CEFit() { return *( (t_CEBase*) this ); }
        //! Returns a constant reference to the CE fitting base.
        const t_CEBase& CEFit() const { return *( (t_CEBase*) this ); }
        //! Returns a reference to the CE fitting base.
        t_ColBase& Collapse() { return *( (t_ColBase*) this ); }
        //! Returns a constant reference to the CE fitting base.
        const t_ColBase& Collapse() const { return *( (t_ColBase*) this ); }


        //! Creates the fitting matrix and target vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &_A, T_VECTOR &_b,
                           types::t_unsigned _dim );
        //! Evaluates square errors.
        opt::ErrorTuple evaluate();

        //! Updates the separable and copies the eci from column 0 to all other columns.
        void update_all();
        //! Updates the separable and copies the eci from column d to column 0.
        void update( types::t_unsigned _d );
        //! Resets collapse functor.
        void reset() { t_ColBase::reset(); }

        //! Reference to configuration matrix.
        const t_iMatrix& configurations() const { return t_ColBase::configurations(); }
        
        //! Returns the number of dimensions.
        size_t dimensions() const { return t_ColBase::dimensions(); }
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const { return t_ColBase::dof() + t_CEBase::dof(); }
       
        //! Initializes the mixed approach.
        void init();

        //! Returns a reference to the coefficients.
        t_OMatrix& coefficients() { return coefficients_; }
        //! Returns a constant reference to the coefficients.
        const t_OMatrix& coefficients() const { return coefficients_; }
        //! Returns a reference to the mapping.
        typename t_ColBase::t_Mapping& mapping() { return t_ColBase::mapping(); }
        //! Returns a reference to the mapping.
        const typename t_ColBase::t_Mapping& mapping() const
          { return t_ColBase::mapping(); }
        //! Returns a reference to the regularization.
        typename t_ColBase::t_RegPolicy& regularization() 
          { return t_ColBase::regularization(); }
        //! Returns a constant reference to the regularization.
        const typename t_ColBase::t_RegPolicy& regularization() const
          { return t_ColBase::regularization(); }
        //! Returns a reference to the separable function;
        t_Separables& separables() { return separables_; }
        //! Returns a constant reference to the separable function;
        const t_Separables& separables() const { return separables_; }

      protected:
        //! Creates the _A and _b matrices for fitting.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );

        //! The separable function.
        t_Separables separables_;

        //! The coefficients.
        t_OMatrix coefficients_;
    };
} // end of CE namspace

#include "mixed.impl.h"

#endif

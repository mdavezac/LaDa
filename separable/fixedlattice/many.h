//
//  Version: $Id$
//
#ifndef _CE_MANY_H_
#define _CE_MANY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/ptr_container/ptr_list.hpp>
#include<boost/shared_ptr.hpp>

#include<iostream>
#include<vector>
#include<list>

#include <opt/types.h>
#include <opt/debug.h>

//! \cond 
namespace CE { template< class T_TRAITS > class Many; }
//! \endcond

namespace Traits
{
  namespace CE
  {
    //! Traits of a "Many" collapse functor.
    template< class T_SEPARABLES,
              class T_COLLAPSE,
              class T_MAPPING = ::CE::Mapping::Basic >
    struct Many 
    {
      //! Type of the Mapping.
      typedef T_SEPARABLES t_Separables;
      //! Type of the Mapping.
      typedef T_MAPPING t_Mapping;
      //! Type of the Regulations Policy
      typedef T_COLLAPSE t_Collapse;
    };

    //! \brief Traits for a combination of Cluster Expansion with Sum of
    //!       Separable Functions.
    template< class T_MANYTRAITS,  class T_CEBASE >
     struct MixedManyApproach
     {
       //! CE::Many traits.
       typedef T_MANYTRAITS t_ManyTraits;
       //! Original separable traits.
       typedef typename t_ManyTraits :: t_Separables :: t_Traits t_OrigSepTraits;
       protected:
         //! Original separable traits.
         typedef typename t_ManyTraits :: t_Collapse :: t_Traits t_CollapseTraits;
         //! Separable function traits.
         typedef Separables
                 < 
                   typename t_OrigSepTraits :: t_Mapping,
                   typename t_OrigSepTraits :: t_Policy,
                   ::CE::Policy::MatrixRangeCoefficients,
                   typename t_OrigSepTraits :: t_Vector
                 > 
                 t_SepTraits;
       public:
         //! Type of the separable function.
         typedef ::CE::Separables< t_SepTraits > t_Separables;
         //! Type of the configuration matrix.
         typedef typename t_CollapseTraits ::  t_Configurations t_Configurations;
         //! Type of the Mapping.
         typedef typename t_ManyTraits :: t_Mapping t_Mapping;
         //! Type of the Regulations Policy
         typedef typename t_CollapseTraits :: t_RegPolicy
                                      ::template rebind< t_Separables > :: other 
            t_RegPolicy;
         //! Type of the Policy.
         typedef typename t_CollapseTraits :: t_UpdatePolicy
                             ::template rebind< t_Separables,
                                             typename t_CollapseTraits :: t_Mapping, 
                                             typename t_CollapseTraits :: t_Configurations >
                             :: other t_UpdatePolicy;
       protected:
         //! collapse traits.
         typedef Collapse
                 < 
                   t_Separables,
                   typename t_CollapseTraits :: t_Mapping,
                   t_RegPolicy,
                   t_Configurations,
                   t_UpdatePolicy
                 > 
                 t_NewCollapseTraits;
         //! New Many traits.
         typedef Many< t_Separables,
                       ::CE::Collapse<t_NewCollapseTraits>,
                       t_Mapping > t_NewManyTraits;
       public:
         //! Type of the collapse functor base.
         typedef ::CE::Many< t_NewManyTraits > t_Collapse;
         //! CE fit base.
         typedef T_CEBASE t_CEBase;
     };
  }
} // end of traits namespace.

namespace CE
{
  template< class T_TRAITS >
    class Many 
    {
      public:
        //! Type of the traits.
        typedef T_TRAITS t_Traits;
        //! Type of the separable collapse functor.
        typedef typename t_Traits :: t_Collapse t_Collapse;
        //! \brief Type of the matrix range.
        //! \details Necessary interface for minimizer.
        typedef typename t_Traits :: t_Collapse :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Traits :: t_Collapse :: t_Vector t_Vector;
        //! Type of the container of separables.
        typedef boost::ptr_list< t_Collapse > t_Collapses;
        //! Type of the of separables function.
        typedef typename t_Traits :: t_Separables  t_Separables;
        //! Type of the container of separables.
        typedef boost::ptr_list< t_Separables > t_CtnrSeparables;
        //! Type of the general mapping.
        typedef typename t_Traits :: t_Mapping t_Mapping;


        //! Constructor.
        Many() : separables_( new t_CtnrSeparables ),
                 collapses_( new t_Collapses ), dim(0) {}
        //! Copy Constructor.
        Many( const Many& _c ) : separables_( _c.separables_ ),
                                 collapses_( _c.collapses ),
                                 dim( _c.dim ) {}
        //! Destructor.
        ~Many() {}

        //! Creates the fitting matrix and target vector.
        template< class T_MATRIX, class T_VECTOR >
          void operator()( T_MATRIX &_A, T_VECTOR &_b,
                           types::t_unsigned _dim );
        //! Evaluates square errors.
        opt::ErrorTuple evaluate() const;
        //! Predicts target value of a structure.
        typename t_Matrix :: value_type evaluate( size_t _n ) const;

        //! Updates the separable and copies the eci from column 0 to all other columns.
        void update_all();
        //! Updates the separable and copies the eci from column d to column 0.
        void update( types::t_unsigned _d );
        //! Resets collapse functor.
        void reset();

        //! Returns the number of dimensions.
        size_t dimensions() const;
        //! Returns the number of degrees of liberty (per dimension).
        size_t dof() const;
        //! Returns the number of configurations.
        size_t nbconfs() const;
       
        //! Randomizes both cluster energies and ecis.
        void randomize( typename t_Vector :: value_type _howrandom );

        //! Add new collapse and separables.
        size_t addone();
        
        //! Returns reference to nth separable function.
        t_Separables& separables( size_t _n ) { return separables_[_n]; }
        //! Returns constant reference to nth separable function.
        const t_Separables& separables( size_t _n ) const { return separables_[_n]; }
        //! Returns reference to nth collapse functor.
        t_Collapse& collapse( size_t _n ) { return collapses_[_n]; }
        //! Returns constant reference to nth collapse functor.
        const t_Collapse& collapse( size_t _n ) const { return collapses_[_n]; }
        //! Returns the number of collapse and separables functions.
        size_t size() const { return collapses_.size(); }
        //! Returns a reference to the mapping.
        t_Mapping mapping() { return mapping_; }
        //! Returns a constant reference to the mapping.
        const t_Mapping mapping() const { return mapping_; }

      protected:
        //! Returns the number of degrees of liberty for current dimension.
        size_t current_dof() const;
        //! Creates the _A and _b matrices for fitting.
        template< class T_MATRIX, class T_VECTOR >
          void create_A_n_b( T_MATRIX &_A, T_VECTOR &_b );
        //! The container of separable functions.
        boost::shared_ptr<t_CtnrSeparables> separables_;
        //! The collapse functor associated with the separable functions.
        boost::shared_ptr<t_Collapses> collapses_;
        //! Current dimension being updated.
        size_t dim;
        //! The mapping to the structures ( e.g. leave-one-out, leave-many-out )
        t_Mapping mapping_;
    };

  //! Prints mixed-approach description to a stream.
  template< class T_TRAITS >
  std::ostream& operator<<( std::ostream& _stream, const Many<T_TRAITS> &_col );

  //! Initializes a Many separable function depending on string input and structures.
  template< class T_STRUCTURES, class T_TRAITS >
   void init_many_collapses( const std::string &_desc, size_t _rank, types::t_real _lambda,
                             const T_STRUCTURES &_structures, Many<T_TRAITS> &_many );

} // end of CE namespace

#include "many.impl.h"

#endif

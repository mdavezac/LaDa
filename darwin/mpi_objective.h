//
//  Version: $Id$
//
#ifndef _MPIOBJECTIVE_H_
#define _MPIOBJECTIVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi/mpi_object.h>

#include "objective.h"
#include "farm.h" 


namespace Objective 
{
  namespace mpi 
  {
    template<class T_GA_TRAITS >
    class Bull : private GA::mpi::Bull,
                 public Types<T_GA_TRAITS>::Vector
    {
        //! %MPI Base of this class
        typedef GA::mpi::Bull                               t_MPIBase; 
        //! %Base of this class
        typedef Container<T_GA_TRAITS>                      t_Base; 
      public:
        typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
      protected:
        //! Type of individual in this %GA
        typedef typename t_GATraits :: t_Individual         t_Individual;
        //! Type of the fitness, as declared in the base class
        typedef typename t_Base :: t_Fitness                t_Fitness;
        //! Type of the quantity, as declared in the base class
        typedef typename t_Base :: t_Quantity               t_Quantity;
        //! Type of the scalar quantity, as declared in the base class
        typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
        //! Type of the lamarckian traits, as declared in the base class
        typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
        //! Type of the lamarckian gradients, as declared in the base class
        typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
        //! Type of the lamarckian variables, as declared in the base class
        typedef typename t_VA_Traits :: t_Type              t_VA_Type;
        //! %Traits of the quantity
        typedef typename t_GATraits :: t_QuantityTraits     t_QuantityTraits;
    
      protected:
        using t_Base :: fitness;
        
    
      public:
        //! Constructor
        Cow   ( MPI::Comm _farmercomm &_comm, MPI::Comm _cowcomm )
            : t_MPIBase( _farmercomm, _cowcomm ), t_Base() {}
        //! Copy Constructor
        Cow ( const Cow &_c ) : t_MPIBase(_c), t_Base(_c) {}
        //! Destructor
        virtual ~Cow() {}
    
        //! Returns fitness as obtained from Bull
        virtual const t_Fitness& operator()(const t_Quantity& _q);
        //! Returns fitness and gradients as obtained from Bull
        virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity & _q,
                                                         t_QuantityGradients&,
                                                         t_VA_Type *_i_grad);
        //! Requests gradients from Bull
        virtual void evaluate_gradient( const t_Quantity &_q,
                                        t_QuantityGradients &_grad,
                                        t_VA_Type *_i_grad);
        //! Requests gradient from Bull
        virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                                 t_QuantityGradients& _grad,
                                                 types::t_unsigned _n);
        //! Returns "mpi::Cow"
        virtual std::string what_is() const { return "mpi::Cow" };
       
    };

    template<class T_GA_TRAITS >
    class Cow : private GA::mpi::Cow,
                public Types<T_GA_TRAITS>::Vector
    {
        //! %MPI Base of this class
        typedef GA::mpi::Cow                                t_MPIBase; 
        //! %Base of this class
        typedef Container<T_GA_TRAITS>                      t_Base; 
      public:
        typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
      protected:
        //! Type of individual in this %GA
        typedef typename t_GATraits :: t_Individual         t_Individual;
        //! Type of the fitness, as declared in the base class
        typedef typename t_Base :: t_Fitness                t_Fitness;
        //! Type of the quantity, as declared in the base class
        typedef typename t_Base :: t_Quantity               t_Quantity;
        //! Type of the scalar quantity, as declared in the base class
        typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
        //! Type of the lamarckian traits, as declared in the base class
        typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
        //! Type of the lamarckian gradients, as declared in the base class
        typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
        //! Type of the lamarckian variables, as declared in the base class
        typedef typename t_VA_Traits :: t_Type              t_VA_Type;
        //! %Traits of the quantity
        typedef typename t_GATraits :: t_QuantityTraits     t_QuantityTraits;
    
      protected:
        using t_Base :: fitness;
        
    
      public:
        //! Constructor
        Cow( MPI::Comm _comm &_comm ) : t_MPIBase( _comm ), t_Base() {}
        //! Copy Constructor
        Cow ( const Cow &_c ) : t_MPIBase(_c), t_Base(_c) {}
        //! Destructor
        virtual ~Cow() {}
    
        //! Returns fitness as obtained from Bull
        virtual const t_Fitness& operator()(const t_Quantity& _q);
        //! Returns fitness and gradients as obtained from Bull
        virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity & _q,
                                                         t_QuantityGradients&,
                                                         t_VA_Type *_i_grad);
        //! Requests gradients from Bull
        virtual void evaluate_gradient( const t_Quantity &_q,
                                        t_QuantityGradients &_grad,
                                        t_VA_Type *_i_grad);
        //! Requests gradient from Bull
        virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                                 t_QuantityGradients& _grad,
                                                 types::t_unsigned _n);
        //! Returns "mpi::Cow"
        virtual std::string what_is() const { return "mpi::Cow" };
       
    };


  } // namespace mpi
} // namespace Objective

#include "mpi_objective.impl.h"
#endif

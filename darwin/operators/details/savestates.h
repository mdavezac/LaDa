//
//  Version: $Id: minimizergenop.h 844 2008-11-08 01:22:54Z davezac $
//
#ifndef _LADA_GA_OPERATOR_DETAILS_SAVESTATES_H_
#define _LADA_GA_OPERATOR_DETAILS_SAVESTATES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace LaDa
{
  namespace GA 
  {
    namespace Operator
    {
      namespace details
      {
        //! \brief Functor for saving current minimization state
        //! \details For any number of reasons the minimizer may want to go back
        //!          to a previous state. This functor is used by minimization
        //!          functors to save the current state and possibly reload it later
        //!          on. In practice, and in the case Virtual-Atom like minimization,
        //!          saving the "state" means saving both the quantity and the
        //!          fitness. It is assumed that the variables of the functionals for
        //!          the saved state is recovered by the minizer itself. For an
        //!          example, see minimizer::VA.
        template< class T_GATRAITS >
        class SaveStateIndividual
        {
          public:
            //! All %types relevant to %GA.
            typedef T_GATRAITS t_GATraits;
          protected:
            //! The type of the individual.
            typedef typename t_GATraits :: t_Individual                  t_Individual;
            //! All %types pertaining to the quantity
            typedef typename t_GATraits :: t_QuantityTraits              t_QuantityTraits;
            //! The type of the quantity
            typedef typename t_QuantityTraits :: t_Quantity              t_Quantity;
            //! The type of the scalar quantity
            typedef typename t_QuantityTraits :: t_ScalarQuantity        t_ScalarQuantity;
            //! The type of the scalar fitness.
            typedef typename t_GATraits :: t_Fitness :: t_ScalarFitness  t_Fitness;
 
          protected:
            //! Pointer to the individual begin  minimized
            t_Individual *current_indiv;
            //! Instance of t_Quantity for saving the state
            t_Quantity quantity;
            //! Instance of t_Fitness for saving the state
            t_Fitness fitness;
 
          public:
            //! Constructor
            SaveStateIndividual() : current_indiv(NULL) {};
            //! Constructor and Initializer
            SaveStateIndividual( t_Individual &_c) : current_indiv(_c) {};
            //! Copy Constructor
            SaveStateIndividual   ( const SaveStateIndividual &_c) 
                                : current_indiv(_c.current_indiv), quantity( _c.quantity ),
                                  fitness( _c.fitness ) {};
            //! Destructor
            ~SaveStateIndividual() {};
            
            //! Saves the state.
            void save();
            //! resets the state to the saved state.
            void reset() const;
            //! sets SaveStateIndividual::current_indiv to \a _indiv.
            void init( t_Individual &_indiv ) { current_indiv = &_indiv; }
        };

    template< class T_GATRAITS >
      inline void SaveStateIndividual<T_GATRAITS> :: save()
      {
        if ( not current_indiv ) return;
        quantity = current_indiv->const_quantities();
        fitness =  current_indiv->fitness();
      }
    template< class T_GATRAITS >
      inline void SaveStateIndividual<T_GATRAITS> :: reset() const 
      {
        if ( not current_indiv ) return;
        current_indiv->quantities() = quantity;
        current_indiv->fitness() = fitness;
      }


      } // namespace details
    } // namespace Operator
  } // namespace GA
} // namespace LaDa

#endif // _DARWIN_MINMIZER_GENOP_H_

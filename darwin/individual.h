//
//  Version: $Id$
//
#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <exception>
#include <iostream>
#include <stdexcept>       // std::runtime_error
#include <algorithm>
#include <functional>

#include <tinyxml/tinyxml.h>

#include <eo/eoObject.h>      // eoObject
#include <eo/eoPersistent.h>  // eoPersistent
#include <eo/eoOpContainer.h>  


#include "opt/opt_function_base.h"
#include "opt/types.h"

#include "objective.h"
#include "fitness.h"
#include "gatraits.h"


#ifdef _MPI
#include "mpi/mpi_object.h"
#endif 

/** \ingroup Genetic
 * @{
 * \brief Declares template classes for individuals
 * \details Individuals are those objects in Genetic Algorithms which have the
 * ability to mate, live, die, and evolve. The goal of any Genetic Algorithm is
 * to evolve a population of individuals towards the solution of the problem.
 * As such, inidividuals form one of the basis of the Algorithms. The base
 * class below take care of defining the fitness of an individual, as well as
 * ordering between individuals according to these fitnesses. The actual
 * meaning of each individual, in the sense of what is they code for is left to
 * the templated type \a Base::t_IndivTraits::t_Object. Furthermore, Individual
 * allows saving and loading to and from XML input, using templatized save and
 * load functors. 
 * \param T_INDIVTRAITS declares all the types relevant to an Individual, and then some.
 *        The two truly relevant types are 
 *        - Base::t_IndivTraits::t_Object declares the type of object the individual represents
 *        - Base::t_IndivTraits::t_Fitness declares the type of fitness to use, eg a
 *        scalar Fitness::Base, or a vectorial Fitness::Base type.
 *        - Base::t_IndivTraits::t_Quantity declares the type of the individuals quantities.
 *        .
 *
 *  \par More on Base::t_IndivTraits::t_Fitness versus Base::t_IndiviTraits::t_Quantity
 *  This last type, Base::t_IndivTraits::t_Quantity, appears when separating an
 *  individuals fitness and its  "quantities", such as formation energy,
 *  band-gap, etc... The two are linked <em>via</em>  the functors in the
 *  Objective namespace. More specifically, these "quantities" represent the
 *  raw fitness straight out of the oven which are the Evaluator classes.
 *  Base::t_IndivTraits::t_Fitness, obtained from  Base::t_IndivTraits::t_Quantity through
 *  the application of an Objective functor, represent the polished fitness
 *  against which individuals are measured and ranked.
 */
namespace Individual
{
  //! \brief Declares a standard individual
  //! \details Should be used with scalar fitnesses only.
  //! This class declares an individual, eg the union of an object, of quantities, and a
  //! fitness. The type of these depend upon the template argument \a
  //! Base::t_IndivTraits::t_Object, \a Base::t_IndivTraits::t_Quantity, and \a
  //! Base::t_IndivTraits::t_Fitness. The first type completely encodes the individuals
  //! <em>via</em>, for intance, a bitstring.  The second type represents the raw
  //! fitness, eg the formation enthalpy of an individual. The third term
  //! represent the fitness against which individuals are judged, eg the depth to
  //! the convex-hull. 
  //! \par 
  //! This class implements the behavior required by EO. It must be used for
  //! scalar fitnesses only, though it is tweaked such that we can easily derive
  //! a new class for vectorial fitnesses. In any case, it is best to acces
  //! this class <em>via</em> the templated type Inidividual::Types::Scalar.
  template<class T_INDIVTRAITS>
  class Base : public eoObject, public eoPersistent
  {
    public:
      //! Individual Traits type \sa Traits::Indiv, Individual::Types::Scalar
      typedef T_INDIVTRAITS t_IndivTraits; 
      typedef typename t_IndivTraits::t_Fitness t_Fitness; //!< Fitness type
      typedef typename t_IndivTraits::t_ScalarFitness Fitness; //!< Scalar Fitness type, for EO.
    protected:
      typedef typename t_IndivTraits :: t_Object t_Object; //!< Object type
      //! Quantity traits type, \sa Traits::Quantity
      typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
      //! Scalar Quantity type 
      typedef typename t_QuantityTraits :: t_ScalarQuantity     t_ScalarQuantity;
      //! Quantity type 
      typedef typename t_QuantityTraits :: t_Quantity           t_Quantity;  

    private:
      typedef Base<t_IndivTraits> t_This; //!< This class type

    protected:
      t_Object object; //!< Instance of t_Object defining this individual
      t_Fitness repFitness; //!< Fitness of this instance
      types::t_unsigned age; //!< Age (number of generations in existence)
      t_Quantity quantity;  //!< t_Quantity instance, eg raw fitness

    public: 
       //! Constructor, Base::age is set to 0
      Base() : age(0) {}
      //! Destructor
      ~Base() {} 
      //! Copy Constructor
      Base   (const t_This &_indiv ) 
           : object( _indiv.object ), 
             repFitness(_indiv.repFitness),
             age(_indiv.age), quantity(_indiv.quantity) {}

      //! \brief Returns true if the (polished) fitness has been set
      //! \details Keeps track of which individuals have been evaluated
      //! already, so that evaluation is not done twice, and that unevaluated
      //! individuals do not go through the cracks to become part of the
      //! population with unitialized fitness.
      //! \sa Fitness::Base::invalid, Evaluation
      bool invalid() const { return (bool) repFitness.invalid(); }
      //! \brief Invalidates current fitness
      //! \details Helps to keep track of which individuals have been evaluated
      //! already, so that evaluation is not done twice, and that unevaluated
      //! individuals do not go through the cracks to become part of the
      //! population with unitialized fitness.
      //! \sa Fitness::Base::invalid, Evaluation
      void invalidate() { repFitness.invalidate(); }

      //! \brief Copies everythin but Base::age
      //! \details Base::age is set to zero.
      //!          This function is useful if you want a new individual (age=0)
      //!          in the population, but which happens to be a clone. \sa
      //!          Taboo::History
      void clone( const t_This &_indiv );

      //! Sets the age of the indiviudal to \a _age 
      void  set_age( types::t_unsigned _age )  { age = _age; }
      //! Returns the the age of the indiviudal
      types::t_unsigned  get_age() const { return age; }

      //! \brief Ordering operator, done according to Fitness::operator<( const Fitness& )
      //! \details Watch that t_Fitness and Fitness may be different types.
      //!          Indeed, as of revision 336, and in the case of
      //!          multi-objective fitnesses, t_Fitness is a \e vectorial fitness,
      //!          whereas Fitness is a \e scalar fitness. This means that as far
      //!          as EO is concerned  only \e scalar fitnesses are compared.
      bool operator<(const t_This& _eo2) const
        { return repFitness.Fitness::operator<( _eo2.repFitness); }
      //! \brief. Ordering operator, done according to Fitness::operator>( const Fitness& )
      //! \details Watch that t_Fitness and Fitness may be different types.
      //!          Indeed, as of revision 336, and in the case of
      //!          multi-objective fitnesses, t_Fitness is a \e vectorial fitness,
      //!          whereas Fitness is a \e scalar fitness. This means that as far
      //!          as EO is concerned  only \e scalar fitnesses are compared.
      bool operator>(const t_This& _eo2) const
        { return repFitness.Fitness::operator>( _eo2.repFitness); }
      //! \brief Equalitey operator
      //! \details
      //! - If Base::invalid is false, two individuals are equal if 
      //!                              and only if, and in this order,  their
      //!                              concentrations are different, and their 
      //!                              Base::object are different.
      //! - If Base::invalid is true, two individuals are equal if and
      //!                             only if, and in this order, their
      //!                            concentration are the same, their
      //!                            (polished) fitnesses Base::repFitness are
      //!                            different, and their Base::object are different.
      //! .
      //! The concentration is obtained by Base::get_concentration(). Objects
      //! could be compared using Base::object only, but concentration and
      //! fitness should suffice in most cases, and should be much faster.
      bool operator==( const t_This &_indiv ) const;
        
      //! Returns a reference to Base::repFitness
      t_Fitness& fitness()  { return repFitness; }
      //! Returns a constant reference to Base::repFitness
      const t_Fitness& fitness() const;
      //! Returns a constant reference to Base::repFitness
      const t_Fitness& const_fitness() const { return fitness(); }
      //! Sets the fitness using a t_Fitness
      void set_fitness( const t_Fitness &_q ) { repFitness = _q; }
      //! Sets the fitness using a t_Quantity
      void set_fitness( const t_Quantity _q ) { repFitness = _q; }

      //! EO stuff
      virtual std::string className() const
        { return "Individual::Base<...>"; }
      //! More useless EO stuff
      void printOn(std::ostream &_os) const;
      //! Still more useless EO stuff
      void readFrom(std::istream &_is) {}

      //! Prints out Base::object to a stream
      void print_out( std::ostream &_stream ) const
        { _stream << (const t_Object&) object << "  "; }

      //! returns a reference to Base::object
      operator t_Object& ()  { return object; }
      //! returns a constant reference to Base::object
      operator const t_Object& () const  { return object; }
      //! returns a reference to Base::object
      t_Object& Object() { return object; }
      //! returns a constant reference to Base::object
      const t_Object& Object() const { return object; }


      //! returns the concentration of Base::object
      types::t_real get_concentration() const
        { return object.get_concentration(); }

      //! returns a reference to Base::quantity
      t_Quantity& quantities() { return quantity; }
      //! returns a constant reference to Base::quantity
      const t_Quantity& const_quantities() const { return quantity; }

      //! \brief returns a reference to the \a _n(th) component of Base::quantity
      //! \details In the case of a scalar quantity, this always returns a
      //! reference to Base::quantity, whatever the value of \a _n
      t_ScalarQuantity& quantities( types::t_unsigned _n )
        { return t_QuantityTraits::scalar(quantity,_n); }
      //! \brief returns a constant reference to the \a _n(th) component of Base::quantity
      //! \details In the case of a scalar quantity, this always returns a
      //! constant reference to Base::quantity, whatever the value of \a _n
      const t_ScalarQuantity& const_quantities( types::t_unsigned _n ) const
        { return t_QuantityTraits::scalar(quantity,_n); }

      //! \brief Save to XML
      //! \param _node  XML node to which to save this instance
      //! \param _saveop a functor capable of saving a Base::t_Object to XML
      template<class SaveOp>  bool Save( TiXmlElement &_node, SaveOp &_saveop ) const;
      //! \brief Load from XML
      //! \param _node  XML node from which to load this instance
      //! \param _loadop a functor capable of loading a Base::t_Object to XML
      template<class LoadOp> bool Load( const TiXmlElement &_node, LoadOp &_loadop );

#ifdef _MPI
      /** \ingroup MPI
       * \brief Serializes Individual::Base class for mpi purposes
       * \details It serializes Base::age, Base::repFitness, Base::object,
       * Base::quantity.
       */
      bool broadcast( mpi::BroadCast &_bc );
#endif
  };

  //! \brief Multi-Objective individual
  //! \details Multi overides a few functions of its base class \a T_BASE in
  //! order to adapt it to use with a vectorial quantity and fitness. The best
  //! way to use this class  is through Individual::Types::Vector.
  //! \param T_BASE is meant to be some instanciation of Individual::Base
  //! \warning \a T_BASE is quite strange a templatization: you should not use
  //! this class directtly, unless you know what you are doing.
  template< class T_BASE >
  class Multi : public T_BASE
  {
      typedef T_BASE t_Base; //!< Base class type
    public: 
      //! \brief Contains all traits pertaining to an individual, and then some
      //! \sa Traits::Indiv
      typedef typename t_Base::t_IndivTraits t_IndivTraits; 

    protected:
      //! Scalar fitness type
      typedef typename t_Base::t_Fitness::t_ScalarFitness t_ScalarFitness;
      //! Scalar Quantity type
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

    protected:
      using t_Base :: repFitness; 

    public:
      using t_Base :: set_fitness;

      //! \brief Sets the scalar fitness of Multi::repFitness
      //! \details For a scalar fitness, this function is exactly equivalent to
      //! Base::set_fitness( const t_Quantiy ) and hence cannot even be declared.
      //! Actually, this fact explains why Multi exists in the first place.
      void set_fitness( const t_ScalarFitnessQuantity _q )
        { (t_ScalarFitness&) repFitness = _q; }

      //! \brief Save to XML
      //! \details, we need to redeclare this function since \a _saveop cannot take
      //! an Individual::Base instance as an argument, which it would if we were
      //! to use Base::Save() directly. 
      //! \param _node  XML node to which to save this instance
      //! \param _saveop a functor capable of saving a Base::t_Object to XML
      //! \todo Find a better way to do this, such that we can use Base::Save directly
      template<class SaveOp> bool Save( TiXmlElement &_node, SaveOp &_saveop ) const;
      //! \brief Load from XML
      //! \details, we need to redeclare this function since \a _loadop cannot take
      //! a Individual::Base instance as an argument, which it would if we were
      //! to use Base::Load() directly.
      //! \param _node  XML node to which to save this instance
      //! \param _loadop a functor capable of loading a Base::t_Object from XML
      //! \todo Find a better way to do this, such that we can use Base::Load directly
      template<class LoadOp> bool Load( const TiXmlElement &_node, LoadOp &_loadop );
  };

  //! \brief Sends a Individual::Base to a stream \a _str using Base::print_out()
  template<class T_INDIVTRAITS>
  inline std::ostream& operator<< ( std::ostream &_str,
                                    const Base<T_INDIVTRAITS> &_indiv )
    { _indiv.print_out(_str); return _str; }

  //! \brief Makes it simpler to declare a scalar and vectorial individual for
  //! crystal objects 
  //! \details A crystal object is simply any physical structure. This class is
  //! meant to be used in conjunction with all the operators declared in
  //! gaoperators.h
  //! To declare a Scalar individual use as in ce.h,
  //! \code
  //!   typedef Individual::Types< SingleSite::Object, 
  //!                              SingleSite::Concentration, 
  //!                              SingleSite::Fourier        > :: Scalar t_Individual;
  //! \endcode
  //! Where it is expected that the template quantities are declared.
  //! To declare a Vectorial individual use as in molecularity.h,
  //! \code
  //!   typedef Individual::Types< Object, 
  //!                              Concentration, 
  //!                              Fourier        > :: Vector t_Individual;
  //! \endcode
  //! \param T_OBJECT the type of the object 
  //! \param T_CONCENTRATION the type of a functor capable of computing the
  //!          concentration of a \a T_OBJECT instance, as well as setting it
  //!          concentration in some cases. 
  //! \param T_FOURIER_RTOK the type of a functor capable of computing the
  //!         Fourier transform and the back Fourier transform of an object
  //! \see SingleSite::Fourier, SingleSite::Object, SingleSite::Concentration...
  template<class T_OBJECT, class T_CONCENTRATION, class T_FOURIER_RTOK>
    class Types
    {
      public:
        typedef T_OBJECT        t_Object; //!< object type 
        typedef T_CONCENTRATION t_Concentration; //!< concentration functor type
        typedef T_FOURIER_RTOK  t_FourierRtoK; //!< Fourier functor type

      protected:
        typedef t_FourierRtoK t_FourierKtoR; //!< Back Fourier transform functor
        //! Traits of scalar individual whose quantities are real
        typedef Traits::Indiv< t_Object, t_Concentration, t_FourierRtoK, t_FourierKtoR,
                               Traits::Quantity< types::t_real > > t_realtraits;
        //! Traits of vectorial individual whose quantities are real
        typedef Traits::Indiv< t_Object, t_Concentration, t_FourierRtoK, t_FourierKtoR,
                               Traits::Quantity< std::vector<types::t_real> > > t_vectortraits;
        //! Base class for templatized Multi 
        typedef Base< t_vectortraits > VectorBase;

      public:
        typedef Base< t_realtraits > Scalar; //!< A scalar individual type for real quantities
        typedef Multi< VectorBase > Vector;  //!< A vectorial individual type for real quantities
    };


} // namespace Individual
/**
 * @} */

#include "individual.impl.h"

#endif // _INDIVIDUAL_H_

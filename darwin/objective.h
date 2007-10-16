//
//  Version: $Id$
//
#ifndef _MULTIOB_OBJECTIVE_H_
#define _MULTIOB_OBJECTIVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <opt/convex_hull.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif
#include <print/xmg.h>
#include <print/manip.h>

#include "gatraits.h"
#include "loadsave.h"

/** \ingroup Genetic
 * @{ */
//! \brief Implements optimization goals by linking "raw fitness" to "polished" fitness
//! \details The raw fitness is simply whatever comes out of the functional (e.g. some
//! class derived from GA::Evaluator). The fitness is whatever quantity is used
//! to appreciate one individual over the other. The link between the two is the
//! objective. This could mean in the case of minimization a simple identity, in
//! the case of maximization a negation, or in the case of a convex-hull,
//!  the depth to that convex-hull. 
//! Note \f$\mathcal{P}(\sigma)\f$ the polished fitness (Fitness) of an
//! individual \f$\sigma\f$, and \f$q(\sigma)\f$ the raw fitness
//! (Individual::Base::quantity) of the same individual, then the objective
//! \f$\mathcal{O}\f$ is simply the %function such that 
//! \f$\mathcal{P}(\sigma) =\mathcal{Q}\left[ q(\sigma) \right]\f$.
//! 
//! All objectives are derived from a single template class, Objective::Base,
//! built to support through templatization both scalar and vectorial
//! fitnesses. In any case, are defined the member functions allowing straight up
//! evaluation, taking in a t_Quantity and returning a t_Fitness, gradient
//! evaluations, as well as standard XML read and write.
//! 
//! For ease of use a class Objective::Types is declared which contains
//! typedefs to scalar and vectorial objectives, Objective::Types::Scalar and
//! Objective::Types::Vector. It also contains a static %function  capable of
//! reading XML input and returning a pointer to one of the objectives defined below. 
//! \note Fitnesses can be made to depend upon the whole population, say for
//! niching, via objects from the Scaling namespace. Scaling are applied after
//! objectives.
//! \see Fitness, Scaling
//! \xmlinput should be all done by Objective::Types::new_from_xml()
//! At present, the allowed \b scalar objectives are the following
//! \code
//     <Objective type = "Maximize"/>
//! \endcode or, \code
//     <Objective type = "Minimize"/>
//! \endcode or, \code
//     <Objective type = "Target" target=?? />
//! \endcode where \a target expects a "target" number,  or, \code
//     <Objective type = "convexhull" />
//! \endcode
//! For scalar objectives, you can input any one of those directly whithin \<GA\> tags,
//! \code
//    <GA>
//      ... other tags
//      <Objective type="??"/>
//      ... other tags
//    </GA>
//! \endcode
//! \n\n For \a Vectorial objectives, there is at present only one option,
//! which should be inputed as follows
//! \code
//    <GA>
//      ... other tags
//        <Objective type="LinearSum" />
//           <Objective type=?? coef="14" />
//           <Objective type=?? coef="2" />
//           ... other scalar objectives
//        </Objective>
//      ... other tags
//    </GA>
//! \endcode
//! It returns the weighted linear average of scalar objectives,
//! with the weight assigned by each \a coef attribute. The (possible) scalar
//! objectives are those listed above.
//! \note It turns out that Pareto ranking is not an Objective... Its a
//! Scaling. Don't get it? doesn't matter. Just input it as given in
//! Scaling::new_from_xml();
namespace Objective
{

  //! \brief Base class for objectives.
  //! \details The goal is to define the necessary behavior for linking the
  //! "raw fitness" of Individual::Base::quantity and the "polished" fitness of
  //! Individual::Base::repFitness. As sich the following member %funtions are
  //! declared virtual:
  //!    - Base::operator()(const t_Quantity&), shouuld return a polished
  //!       fitnes computed from the "raw fitness" passed as an argument.
  //!    - Base::evaluate_gradient, should evaluate the gradient of the "polished" fitness
  //!    - Base::evaluate_with_gradient, should combine the last two behaviors
  //!    - Base::evaluate_one_gradient, should return the gradient in a specified direction
  //!    .
  //! All these behaviors are basically those required of an opt::function
  //! class. Indeed, a lamarckian scheme would expect as much. Nonetheless,
  //! there is one striking difference, namely that Base::operator()(const
  //! t_Quantity&) returns a reference to a fitness. The reason for this lies
  //! in the possibility of a vectorial objective and shall be described in
  //! more details below.
  //! To these member %functions  are added XML input/output behaviors, validity
  //! checks, and initialization routine.
  //! 
  //! As usual, the implementation is made more complex by the possiblity of a
  //! multi-objective %GA. Indeed, a multi-objective implementation would
  //! expect to have objectives designed for vectorial quantities, but also for
  //! scalar quantities. For instance, Maximization, Minimization, and Target
  //! are defined explicitely as scalar objectives which take a scalar quantity
  //! on input and return a scalar fitness (see Fitness namespace). On the
  //! other hand, LinearSum is defined explicitely for as a vectorial objective.
  //! Base has been written to be a base class for either scalar or vectorial
  //! quantities  depending on the actual template instanciation. As such, it
  //! takes as argument
  //! \param T_GA_TRAITS all %GA traits
  //! \param T_QUANTITY_TRAITS traits of the quantity
  //! \param T_VA_TRAITS the lamarckian traits
  //! 
  //! The type of fitness, whether scalar or vectorial, is deduced from \a
  //! T_QUANTITY_TRAITS using the trait class Fitness::GetScalar.
  //! More explicitely, Base::t_Fitness is defined as follows
  //! \code
  //   public: 
  //     typedef T_GA_TRAITS t_GATraits;
  //     typedef T_QUANTITY_TRAITS t_QuantityTraits;
  //   protected:
  //     typedef typename t_GATraits :: t_Fitness             t_RealFitness;
  //     typedef typename Fitness::GetScalar< t_RealFitness, 
  //                                          t_QuantityTraits::is_scalar > :: t_result t_Fitness;
  //! \endcode
  //! As can be seen, this should always lead to a scalar fitness if the
  //! quantity is scalar, and to a vectorial fitness if the quantity is
  //! vectorial.
  //! 
  //! Base contains two static members variables: Base::current_indiv and
  //! Base::fitness. The first is set to point to the whatever individual is
  //! being evaluated in (static) member %function Base::init. It allows for
  //! more complicated objectives than those relying only quanities. For
  //! instancem Objective::ConvexHull needs to know the concentration of an
  //! individual in addition to its quantity. The second static variable is an
  //! instance of a (scalar or vectorial) fitness. It is used to pass on the
  //! result of Base::operator()() without copying a possibly complex quantity
  //! (the fitness) onto the stack. This is probably not very usefull in the
  //! case of a single-objective %GA, where returning a number is sufficient, by
  //! should prove advantageous in the case of a multi-objective fitness. Note
  //! that a correctly defined fitness class should have operator t_Quantity()
  //! defined, so that on the user end, the return of Base::operator()() can be
  //! dealt with directly as a number (or vector).
  //!
  //! The class Objective::Types contains two types, Objective::Types::Scalar
  //! and Objective::Types::Vector which define a scalar and a vectorial
  //! objective. In the case of a single-objective %GA, there is no difference
  //! between Objective::Types::Scalar and Objective::Types::Vector.
  //! \xmlinput Generally, objective tags should have the following syntax
  //! \code  <Objective type="??" attribute="??" /> \endcode
  //! See each implementation for details
  template< class T_GA_TRAITS,
            class T_QUANTITY_TRAITS = typename T_GA_TRAITS :: t_QuantityTraits,
            class T_VA_TRAITS = typename T_GA_TRAITS :: t_VA_Traits >
  class Base
  {
    public: 
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
      //! \brief Traits of this type of objective's quantity
      //! \details this may be traits to a scalar quantity in the case of a
      //! scalar objective, or traits to a vectorial quantity in the case of a
      //! vectorial objective. More explicitely, \a T_QUANTITY_TRAITS defines
      //! wether this objective is scalar or vectorial.
      typedef T_QUANTITY_TRAITS t_QuantityTraits;
      typedef T_VA_TRAITS t_VA_Traits; //!< %Traits to lamarckian quanities (eg gradients)
    protected:
      //! Type of individual 
      typedef typename t_GATraits :: t_Individual          t_Individual;
      //! \brief Type of fitness used in %GA
      //! \details If the %GA is multi-objective, but this objective is meant
      //! to be scalar, then t_RealFitness is only one step into obtaining the
      //! fitness type relevant to this class.
      typedef typename t_GATraits :: t_Fitness             t_RealFitness;
      //! \brief Type of fitness relevant to this class
      //! \warning This definition expects the fitness to be defined such that
      //!          Fitness::GetScalar works. Eg, we need some standard behavior
      //!          implemented in the t_RealFitness class
      typedef typename Fitness::GetScalar< t_RealFitness, 
                                           t_QuantityTraits::is_scalar > :: t_result t_Fitness;
      //! Type of the quantity
      typedef typename t_QuantityTraits :: t_Quantity       t_Quantity;
      //! \brief Type of the scalar quantity
      //! \details In the case of scalar objectives, Base::t_Quantity and
      //!          Base::t_ScalarQuantity are equivalent
      typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
      //! Type of the gradients, for larmarckian %GA
      typedef typename t_VA_Traits :: t_QuantityGradients   t_QuantityGradients;
      //! Type of the input variables (for lamarckian %GA)
      typedef typename t_VA_Traits :: t_Type                t_VA_Type;
      //! \brief A functor for saving individuals to XML
      //! \note This functor will make a call to GA::Evaluator::Save(const
      //! t_Individual &_indiv, TiXmlElement &_node, bool _type).
      typedef GA::SaveObject<t_GATraits>                    t_SaveOp;
      //! \brief A functor for loading individuals from XML
      //! \note This functor will make a call to GA::Evaluator::Load(t_Individual
      //! &_indiv, TiXmlElement &_node, bool _type) const. 
      typedef GA::LoadObject<t_GATraits>                    t_LoadOp;
    protected:
      //! \brief points to individual being assessed
      //! \details Works in correspondence with Base::init(). This pointer is
      //! usefull when an objective needs more information on an individual than
      //! simply the quantities. See for instance Objective::Base::ConvexHull.
      static const t_Individual *current_indiv;
      //! \brief Repository for evaluation result
      //! \details In the case of complex vectorial Objective, returning a
      //! container (of results) on the stack may be quite a penalty. This static
      //! member plays the role of the stack, with the added benefit that it has
      //! all the behaviors of a fitness.
      static t_Fitness fitness;
    public:
      //! \brief true if the objective is a scalar
      bool const static is_scalar = t_QuantityTraits :: is_scalar;
      //! \brief true if the objective is a vector
      bool const static is_vector = t_QuantityTraits :: is_vector;
    public:
      //! Constructor
      Base() {}
      //! Destructor
      virtual ~Base() {}

      //! \brief initializes Base::current_indiv
      void static init(const t_Individual& _indiv);
      //! returns a polished fitness \f$\mathcal{O}(q)\f$ from a raw fitness \f$q\f$
      virtual const t_Fitness& operator()( const t_Quantity& ) = 0;
      //! \brief evaluates the gradient \f$\partial\mathcal{O}(q)\f$ from
      //! \f$q\f$ and \f$\partial q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //! it is not necessarily zero, say in the case of objectives within
      //! objectives  \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual void evaluate_gradient( const t_Quantity & _q,
                                      t_QuantityGradients & _grad,
                                      t_ScalarQuantity *_i_grad) = 0;
      //! \brief evaluates a polished fitness and its gradient \f$\partial \mathcal{O}(q)\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //! it is not necessarily zero, say in the case of objectives within
      //! objectives  \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients &_grad,
                                                       t_VA_Type *_i_grad) = 0;
      //! \brief returns the gradient in specified direction 
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _n direction of the requested gradient. This "direction" is
      //!           merely the index of the array t_QuantityGradients. 
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &_q,
                                               t_QuantityGradients &_grad,
                                               types::t_unsigned _n) = 0;
      //! \brief returns true if the objective is valid
      //! \details Moving target objectives such as Objective::ConvexHull may
      //! change during the evaluation of a population. As such, whenever the
      //! objective changes, it is marked invalid until this %function is called.
      //! See implementation of Evaluation::operator( t_Population&,
      //! t_Population& )
      virtual bool is_valid() const = 0;
      //! Returns a string describing the objective
      virtual std::string what_is() const = 0;
      //! \brief Returns a string describing the imparmanent status of a
      //! "Moving %Target" objective
      //! \details Simple objectives such as Minimization should return an 
      //! empty string. ConvexHull, which is bit more complex, returns
      //! information on the convex-hull. At present (revision 319) 
      virtual std::string print() const = 0;
      //! Save impermanent status to XML
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op) { return true; };
      //! Load imparmentant status from XML
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) { return true; };
      //! Whether status is imparmanent
      virtual bool does_store() const { return false; }
  };
 
  //! helper class to have one single %function from which to create
  //! objectives from XML
  template< class T_TYPE, bool IS_VECTORIAL > struct fork;

  //! \brief Helper class for simpler definition of scalar and vectorial
  //! objectives
  //! \details This class simply defines two types: a vectorial and a scalar
  //! objectives. It only takes one argument on instanciation and deduces those
  //! for Base automatically. Note that if T_GA_TRAITS::t_QuantityTraits defines
  //! a scalar quantity, than Objective::Types::Scalar and
  //! Objective::Types::Vector are equivalent.
  //
  //! Finally, this class also contains a %function capable of creating
  //! objectives from XML input.
  template <class T_GA_TRAITS >
    struct Types
    {
      template< class T_TYPE, bool is_vectorial > friend struct fork; 
      public:
        typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
      protected:
        //! Type of individual 
        typedef typename t_GATraits::t_Individual           t_Individual;
        //! %Traits of the quantity 
        typedef typename t_GATraits::t_QuantityTraits       t_QuantityTraits;
        //! \brief %Traits of the scalar quantity obtained from t_QuantityTraits
        //! \details if Types::t_QuantityTraits describes a scalar type, then
        //! Types::t_QuantityTraits and Types::t_ScalarQauntityTraits are
        //! equivalent
        typedef typename t_QuantityTraits::t_ScalarQuantityTraits t_ScalarQuantityTraits;

      public:
        //! Scalar objective type
        typedef Base< t_GATraits, 
                      t_ScalarQuantityTraits,
                      typename Traits::VA< std::vector<t_ScalarQuantity>,
                                           t_ScalarQuantity > >           Scalar;
        //! \brief Vector objective type
        //! \details If Types::t_QuantityTraits describes a scalar, then
        //! Types::Scalar and Types::Vector are equivalent (and scalars).
        typedef Base< t_GATraits >  Vector;



        //! \brief Creates an objective according to XML input
        //! \return A pointer to an objective if \a _node describes a valid objective, 
        //!         or returns NULL if something goes wrong.
        //! \note The return value always points to a Types::Vector. See
        //! previous note on equivalence or non-equivalence of Types::Vector and
        //! Types::Scalar.
        //! \warning It is the job of the caller (not of Base::new_from_xml() )
        //! to check the status of the return. It is also the job of the caller
        //! to deallocate the returned pointer after use.
        static Vector* new_from_xml( const TiXmlElement &_node )
          { fork<Types, t_QuantityTraits::is_vector > s; return s( _node ); }
      protected:
        //! Creates a scalar objective from XML
        static Scalar* scalar_from_xml( const TiXmlElement &_node );
        //! Creates a vector objective from XML
        static Vector* vector_from_xml( const TiXmlElement &_node );
    };

  //! \brief Implements maximization of a scalar quantity with \f$\mathcal{O}(q) = -q \f$
  //! \details In practice, this mean transforming any quantity \a q into \a -q.
  //! \note This is a <STRONG>scalar</STRONG> objective
  //! \xmlinput No attributes
  //! \code  <Objective type="maximize"/> \endcode
  template< class T_GA_TRAITS >
  class Maximize : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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

    protected:
      using t_Base :: fitness;

    public:
      Maximize() {} //!< Constructor
      Maximize( const Maximize &) {} //!< Copy Constructor
      virtual ~Maximize() {} //!< Virtual destructor
      
      //! \brief retuns a "maximization" fitness, eg \f$-q\f$
      virtual const t_Fitness& operator()(const t_Quantity& _q)
       { fitness = -_q; return fitness; }
      //! \brief evaluates a "maximization" gradient, eg \f$-\partial q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //! it is not necessarily zero, say in the case of objectives within
      //! objectives  \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual void evaluate_gradient( const t_Quantity &_q,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! \brief evaluates both \f$\mathcal{O}(q) = -q\f$ and 
      //!        \f$\partial\mathcal{O}(q) = -\partial q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //!                it is not necessarily zero, say in the case of
      //!                objectives within objectives 
      //!                \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      //! \brief evaluates gradient in one direction only,
      //! \f$\partial_n\mathcal{O}(q) = -\partial_n q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _n direction of the requested gradient. This "direction" is
      //!           merely the index of the array t_QuantityGradients. 
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &_q,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return -_grad[_n]; }
      //! The status is permanent, and this %function always returns true
      bool is_valid() const { return true; }
      //! Returns "Maximize"
      virtual std::string what_is() const { return " Maximize"; }
      //! The status is permanent, and this %function always returns an empty string
      virtual std::string print() const { return ""; }
  };
 
  //! \brief Implements minimization of a scalar quantity, \f$\mathcal{O}(q) = q\f$
  //! \details In practice, this mean doing nothing (or simple identity) since
  //!          minimization is the default behavior of Fitness.
  //! \note This is a <STRONG>scalar</STRONG> objective
  //! \xmlinput There are no attributes
  //! \code <Objective type="minimize"/> \endcode
  template< class T_GA_TRAITS >
  class Minimize : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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

    protected:
      using t_Base :: fitness;

    public:
      Minimize() {} //!< Constructor
      Minimize( const Minimize &) {} //!< Copy Constructor
      virtual ~Minimize() {} //!< Virtual Destructor
      
      //! retuns a "minimization" fitness, eg identity \f$q\f$
      virtual const t_Fitness& operator()(const t_Quantity& _q)
         { fitness = _q; return fitness; }
      //! \brief evaluates a "minimization" gradient, eg original gradient
      //!        \f$\partial q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //!                it is not necessarily zero, say in the case of
      //!                objectives within objectives 
      //!                \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual void evaluate_gradient( const t_Quantity &_q,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! \brief evaluates both \f$\mathcal{O}(q) = q\f$ and 
      //! \f$\partial\mathcal{O}(q) = \partial q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //!                it is not necessarily zero, say in the case of
      //!                objectives within objectives 
      //!                \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      //! \brief evaluates gradient in one direction only,
      //! \f$\partial_n\mathcal{O}(q) = \partial_n q\f$
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _n direction of the requested gradient. This "direction" is
      //!           merely the index of the array t_QuantityGradients. 
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &_q,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return _grad[_n]; }
      //! The status is permanent, and this %function always returns true
      bool is_valid() const { return true; }
      //! Returns "Minimize"
      virtual std::string what_is() const { return " Minimize"; }
      //! The status is permanent, and this %function always returns an empty string
      virtual std::string print() const { return ""; }
  };
  
  //! \brief Implements optimization towards a (scalar) target
  //! \f$q_0\f$ with \f$\mathcal{O}(q)= |q -q_0|\f$
  //! \note This is a <STRONG>scalar</STRONG> objective
  //! \xmlinput There is one \b required attribute \a target, which should be
  //!           an integer or a real, depending on the type of Target::t_ScalarQuantity
  //!           \code <Objective type="target" target="??" /> \endcode
  template< class T_GA_TRAITS >
  class Target : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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
    protected:
      using t_Base :: fitness;

    protected:
      t_ScalarQuantity q_0; //!< Target value of \f$\mathcal{O}(q) = |q - q_0|\f$
    public:
      //! Constructor and Initializer
      Target( t_ScalarQuantity _q_0 ) : q_0( _q_0 ) {}
      //! Copy Constructor
      Target( const Target &_c ) : q_0( _c.q_0 ) {}
      //! Destructor
      virtual ~Target() {}
      
      //! returns distance from target, eg \f$|q-q_0|\f$
      virtual const t_Fitness& operator()(const t_Quantity& _q)
        { fitness = std::abs( _q - q_0 ); return fitness; }
      //! \brief computes the gradient \f$\mathrm{sgn}(q - q_0) \partial q \f$
      //! \details It was chosen to set the gradient to zero for 
      //!          \f$q \approx q_0\f$, where \f$\approx\f$ is exact in the
      //!          case of integer quantities, and "fuzzy" in the case of
      //!          reals.
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //!                it is not necessarily zero, say in the case of
      //!                objectives within objectives 
      //!                \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual void evaluate_gradient( const t_Quantity &_q,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! \brief computes both \f$\mathcal{O}(q)=|q-q_0|\f$ and
      //!        \f$\partial\mathcal{O}(q)=\mathrm{sgn}(q - q_0) \partial q \f$
      //! \details It was chosen to set the gradient to zero for 
      //!          \f$q \approx q_0\f$, where \f$\approx\f$ is exact in the
      //!          case of integer quantities, and "fuzzy" in the case of reals.
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _i_grad the resulting gradient from this objective. Note that
      //!                it is not necessarily zero, say in the case of
      //!                objectives within objectives 
      //!                \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      //! \brief returns the gradient
      //!        \f$\partial_n\mathcal{O}(q)=\mathrm{sgn}(q - q_0) \partial_n q \f$
      //!        in specified direction \a _n
      //! \details It was chosen to set the gradient to zero for 
      //!          \f$q \approx q_0\f$, where \f$\approx\f$ is exact in the
      //!          case of integer quantities, and "fuzzy" in the case of reals.
      //! \param _q the quantity \f$q\f$
      //! \param _grad the gradient \f$\partial q\f$ as computed by some
      //!              derived class of GA::Evaluator
      //! \param _n direction of the requested gradient. This "direction" is
      //!           merely the index of the array t_QuantityGradients. 
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _q,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      //! The status is permanent, and this %function always returns true
      bool is_valid() const { return true; }
      //! The status is permanent, and this %function always returns an empty string
      virtual std::string print() const { return ""; }
      //! Returns "Target (\f$q_0\f$)"
      virtual std::string what_is() const;
  };


  /** \brief Implements the refinement of a (scalar) convex-hull, 
        \f$ \mathcal{O}(q_\sigma) = q_\sigma - C^{(n)}(x_\sigma), \f$
      \details The convex-hull is defined in opt::ConvexHull::Base. Indeed,
        that is the implementation we shall use below.\n The goal of this
        objective is to construct and refine a convex-hull. As such, the
        objective is defined as \f[ \mathcal{O}(q_\sigma) = q_\sigma -
        C^{(n)}(x_\sigma), \f] with \f$\sigma\f$ an individual, \f$q_\sigma\f$ its
        raw fitness, \f$x_\sigma\f$ its concentration, and \f$C^{(n)}(x)\f$ the
        known convex-hull at iteration (generation) \f$n\f$ of the genetic
        algorithm. This definition of the objective as the distance to the
        known convex-hull allows us to refine the convex-hull simultaneously
        throughout the concentration range.
      \note This is a <STRONG>scalar</STRONG> objective
      \xmlinput There are no special attributes
        \code <Objective type="convexhull"/> \endcode
      \xmlrestart This objective can be saved and restarted directly as a
      convex-hull. See opt::ConvexHull
  */
  template< class T_GA_TRAITS >
  class ConvexHull : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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
      //! A functor to for saving individuals
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      //! A functor to for loading individuals
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;
      //! The convex-hull object instanciated for a %GA individual
      typedef opt::ConvexHull::Base<t_Individual>         t_ConvexHull;

    protected:
      t_ConvexHull convexhull;  //!<  The convex-hull instanciation
      //! \brief is true when convex-hull.
      //! \details The convex-hull becomes invalid whenever a new break-point is found.
      //!          It is set to true again whenever ConvexHull::is_valid()
      //!          const is called. Since the latter is a constant %function,
      //!          ConvexHull::valid is mutable.
      mutable bool valid;       
      using t_Base :: fitness;
      using t_Base :: current_indiv;

    public:
      //! Constructor
      ConvexHull() : valid(true) {}
      //! Copy Constructor
      ConvexHull( const ConvexHull &_c ) : convexhull( _c.convexhull ), valid(true) {}
      //! Destructor
      virtual ~ConvexHull() {}
      
      //! \brief Returns the distance of Base::current_indiv to the
      //!     convex-hull, \f$ \mathcal{O}(q_\sigma) = q_\sigma - C^{(n)}(x_\sigma), \f$
      //! \details First, we check whether Base::current_indiv is a (new)
      //!     breaking-point. Then the subroutines returns the distance to the
      //!     (possibly) updated convex-hull.
      virtual const t_Fitness& operator()(const t_Quantity& _q);
      /** \brief Returns the gradient
              \f$ \partial\mathcal{O}(q_\sigma) = \partial q_\sigma - \partial
              C^{(n)}(x_\sigma), \f$
          \details First, we check whether Base::current_indiv is a (new)
              breaking-point. Then the subroutines returns the gradient of the
              distance to the (possibly) updated convex-hull.
          \param _q the quantity \f$q\f$
          \param _grad the gradient \f$\partial q\f$ as computed by some
                       derived class of GA::Evaluator
          \param _i_grad the resulting gradient from this objective. Note that
                         it is not necessarily zero, say in the case of
                         objectives within objectives 
                         \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      */
      virtual void evaluate_gradient( const t_Quantity &_q,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
        { evaluate_with_gradient( _q, _grad, _i_grad ); }
      /** \brief Returns both \f$ \mathcal{O}(q_\sigma) = q_\sigma -
              C^{(n)}(x_\sigma), \f$ \f$ \partial\mathcal{O}(q_\sigma) =
              \partial q_\sigma - \partial C^{(n)}(x_\sigma), \f$
          \details First, we check whether Base::current_indiv is a (new)
              breaking-point. Then the subroutines computes what it should.
          \param _q the quantity \f$q\f$
          \param _grad the gradient \f$\partial q\f$ as computed by some
                       derived class of GA::Evaluator
          \param _i_grad the resulting gradient from this objective. Note that
                         it is not necessarily zero, say in the case of
                         objectives within objectives 
                         \f$\mathcal{O}_1 \circ \mathcal{O}_2(q)\f$.
      */
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_q,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      /** \brief Returns \f$ \partial_n\mathcal{O}(q_\sigma) =
              \partial_n q_\sigma - \partial_n C^{(n)}(x_\sigma), \f$
          \details First, we check whether Base::current_indiv is a (new)
              breaking-point. Then the subroutines computes what it should.
          \param _q the quantity \f$q\f$
          \param _grad the gradient \f$\partial q\f$ as computed by some
                       derived class of GA::Evaluator
          \param _n direction of the requested gradient. This "direction" is
                    merely the index of the array t_QuantityGradients. 
      */
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _q,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      //! \brief returns true if no new breaking-point has been found since the
      //! last call to this routine
      //! \sa ConvexHull::valid, Evaluate::Base::operator(t_Population&, t_Population&)
      bool is_valid() const;
      //! Returns " Convex-Hull"
      virtual std::string what_is() const { return " Convex-Hull"; }
      //! Returns a string defining the current status of the convex-hull, as
      //! given by opt::ConvexHull::Base::print().
      virtual std::string print() const { return convexhull.print(); }
      //! \brief Saves the current status of the convex-hull
      //! \sa opt::ConvexHull::Base::Save
      //! \note This %function will generally make an indirect call to
      //!       GA::Evaluator::Save(const t_Individual&, TiXmlElement &, bool) const
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
        { return convexhull.Save( _node, _op ); };
      //! \brief Restarts using a previouly saved status of the convex-hull
      //! \sa opt::ConvexHull::Base::Load
      //! \note This %function will generally make an indirect call to
      //!       GA::Evaluator::Load(t_Individual&, const TiXmlElement &, bool)
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op)
        { return convexhull.Load( _node, _op ); };
      //! Always returns true since this is an impermanent objective
      virtual bool does_store() const { return true; }
  };


  //! \brief Base class to implement optimization towards a vector of scalar
  //!        objectives
  //! \details This class is still pure virtual and cannot be used as such. It
  //!          merely implements some general routines over a container of scalar
  //!          objectives. As suchm there is not XML input to create this class.
  //! \note This is a <STRONG>vectorial</STRONG> objective
  //! \xmlrestart This objectives saves and restart the scalar objectives in
  //!             the order given Container::objectives.
  template<class T_GA_TRAITS >
  class Container : public Types<T_GA_TRAITS> :: Vector
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Class holding all possible objective %types
      typedef Types<t_GATraits>                          t_ObjectiveType;
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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
      //! A functor to for saving individuals
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      //! A functor to for loading individuals
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;
      //! Scalar objective type
      typedef typename t_ObjectiveType :: Scalar          t_Objective;
      
    protected:
      using t_Base :: fitness;
      //! \brief Container of scalar objectives
      //! \details The objectives are owned by this class
      t_Objectives objectives; 

    public:
      //! Constructor
      Container() {}
      //! Copy Constructor
      Container( const Container &_c ) : objectives(_c.objectives) {}
      //! \brief Destructor
      //! \details Deletes the pointers held within Container::objectives
      virtual ~Container ();

      //! Returns false if one of the scalar objectives is invalid
      bool is_valid() const;
      //! Calls upon each objective to save current status
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op);
      //! Calls upon each objective to restart from a previously saved status
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op);
      //! Returns true if at least one objective is impermanent
      virtual bool does_store() const;
      //! Concatenates into a string the call to each scalar objective 
      virtual std::string print() const;
  };

  /** \brief Implements a weighted linear average of scalar objectives,
         \f$ \mathcal{O}(q) = \sum_i\omega_i\mathcal{O}_i(q)\f$.
      \note This is a \b vectorial objective
      \xmlinput A list of scalar  %Objective tags are expected.
      \code
        <Objective type="LinearSum" >
          <Objective type="minimize" coef="1"/>
          <Objective type="target"  target=0.3 coef="10"/>
        </Objective>
      \endcode
      The order of the scalar objectives should correspond to the order of the
      scalar quantities, as defined  within your class derived from
      GA::Evaluator. Refer to an actual implementation for details. 
      The \a coef  attribute is the weight of that objective in the linear sum.
      \xmlrestart This objectives saves and restart the scalar objectives in
                  the order given by Container::objectives.
  */
  template<class T_GA_TRAITS >
  class LinearSum : public Container<T_GA_TRAITS>
  {
    public:
      typedef T_GA_TRAITS t_GATraits; //!< All %GA traits
    protected:
      //! Class holding all possible objective %types
      typedef Types<t_GATraits>                          t_ObjectiveType;
      //! Type of individual in this %GA
      typedef typename t_GATraits :: t_Individual         t_Individual;
      //! Base of this class
      typedef typename Types<t_GATraits> :: Scalar        t_Base; 
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
      //! A functor to for saving individuals
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      //! A functor to for loading individuals
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;

    protected:
      using t_Base :: fitness;
      using t_Base :: objectives;
      //!  Vector containing the weights of each scalar objective
      std::vector< t_ScalarQuantity > coefs; 

    public:
      //! Constructor
      LinearSum() {}
      //! Copy Constructor
      LinearSum ( const LinearSum &_c ) : t_Base(_c), coefs(_c.coefs) {}
      //! Destructor
      virtual ~LinearSum() {}

      //! \brief Adds a scalar objective to the linear sum
      //! \details The pointer to the scalar objective is owned by this
      //!    instance of LinearSum and will be destroyed by this instance of
      //!    LinearSum.
      //! \param _objective pointer to a scalar objective
      //! \param _coef coefficient in the linear sum
      void add( t_Objective *_objective, t_ScalarQuantity _coef );
      //! Returns a weighted linear average of scalar objectives,
      //! \f$ \mathcal{O}(q) = \sum_i\omega_i\mathcal{O}_i(q)\f$.
      virtual const t_Fitness& operator()(const t_Quantity& _q);
      //! Computes both  \f$ \mathcal{O}(q) = \sum_i\omega_i\mathcal{O}_i(q)\f$ and
      //! \f$ \partial\mathcal{O}(q) = \partial \sum_i\omega_i\partial\mathcal{O}_i(q)\f$.
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *);
      //! Computes the weighted linear average of the gradients 
      //! \f$ \partial\mathcal{O}(q) = \partial \sum_i\omega_i\partial\mathcal{O}_i(q)\f$.
      virtual void evaluate_gradient( const t_Quantity &_q,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! Computes the weighted linear average of the gradients in direction \a _n,
      //! \f$ \partial_n\mathcal{O}(q) = \partial \sum_i\omega_i\partial_n\mathcal{O}_i(q)\f$.
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      //! Returns a string describing each scalar objective and its assigned weight
      virtual std::string what_is() const;
     
  };
  

   
  //! (Vectorial) Hack, don't touch
  template< class T_TYPE >
    struct fork<T_TYPE, true>
    {
      //! (Vectorial) Hack, don't touch
      typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
      {
        typename T_TYPE :: Vector* result;
        result = T_TYPE::vector_from_xml( _node );
        if( not result ) std::cerr << "Could not create multi-objective..." << std::endl;
        return result;
      }
    };
  //! (Scalar) Hack, don't touch
  template< class T_TYPE > 
    struct fork<T_TYPE, false>
    {
      //! (Scalar) Hack, don't touch
      typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
      {
        return T_TYPE::scalar_from_xml( _node );
      }
    };

}  // namespace Objective
 /* @} */
#endif //  _MULTIOB_OBJECTIVE_H_

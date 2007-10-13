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
//! \f$\mathcal{O}\f$ is simply the function such that 
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
//! Objective::Types::Vector. It also contains a static function  capable of
//! reading XML input and returning a pointer to one of the objectives defined below. 
//! \note Fitnesses can be made to depend upon the whole population, say for
//! niching, via objects from the Ranking namespace. Rankings are applied after
//! objectives.
//! \see Fitnesss, Ranking
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
  //! To these member functions  are added XML input/output behaviors, validity
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
      //! \brief returns a polished fitness from a raw fitness
      virtual const t_Fitness& operator()( const t_Quantity& ) = 0;
      //! \brief evaluates the gradient of a polished fitness
      virtual void evaluate_gradient( const t_Quantity &,
                                      t_QuantityGradients &,
                                      t_ScalarQuantity *) = 0;
      //! \brief evaluates a polished fitness and its gradient
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *) = 0;
      //! \brief returns the gradient in specified direction 
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients&,
                                               types::t_unsigned) = 0;
      //! \brief returns true if the objective is valid
      //! \details Moving target objectives such as Objective::ConvexHull may
      //! change during the evaluation of a population. As such, whenever the
      //! objective changes, it is marked invalid until this function is called.
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
 
  //! helper class to have one single function from which to create
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
  //! Finally, this class also contains a function capable of creating
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
      
      //! retuns a "maximization" fitness, eg -\a _val
      virtual const t_Fitness& operator()(const t_Quantity& _val)
       { fitness = -_val; return fitness; }
      //! evaluates a "maximization" gradient, eg negates original gradient
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! evaluates a "maximization" fitness and gradient
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      //! evaluates gradient in one direction only
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return -_grad[_n]; }
      //! The status is permanent, and this function always returns true
      bool is_valid() const { return true; }
      //! Returns "Maximize"
      virtual std::string what_is() const { return " Maximize"; }
      //! The status is permanent, and this function always returns an empty string
      virtual std::string print() const { return ""; }
  };
 
  //! \brief Implements minimization of a scalar quantity, \f$\mathcal{O}(q) = q\f$
  //! \details In practice, this mean doing nothing (or simple identity) since
  //!          minimization is the default behavior of Fitness.
  //! \note This is a <STRONG>scalar</STRONG> objective
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
      
      //! retuns a "minimization" fitness, eg identity
      virtual const t_Fitness& operator()(const t_Quantity& _val)
         { fitness = _val; return fitness; }
      //! evaluates a "minimization" gradient, eg original gradient
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      //! evaluates a "minimization" fitness and gradient
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad);
      //! evaluates gradient in one direction only
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) { return _grad[_n]; }
      //! The status is permanent, and this function always returns true
      bool is_valid() const { return true; }
      //! Returns "Minimize"
      virtual std::string what_is() const { return " Minimize"; }
      //! The status is permanent, and this function always returns an empty string
      virtual std::string print() const { return ""; }
  };
  
  //! \brief Implements optimization towards a (scalar) target
  //! \f$q_0\f$ with \f$\mathcal{O}(q)= |q -q_0|\f$
  //! \note This is a <STRONG>scalar</STRONG> objective
  template< class T_GA_TRAITS >
  class Target : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_GATraits :: t_IndivTraits       t_IndivTraits;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;

    protected:
      using t_Base :: fitness;

    protected:
      t_ScalarQuantity target; 
    public:
      Target( t_ScalarQuantity _target ) : target( _target ) {}
      Target( const Target &_c ) : target( _c.target ) {}
      virtual ~Target() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val)
        { fitness = std::abs( _val - target ); return fitness; }
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
      }
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *_i_grad)  
      {
        typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
        t_VA_Type *i_grad_result = _i_grad;
        if ( _val > target ) 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result += *i_grad;
        else 
          for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
            *i_grad_result -= *i_grad;
        return std::abs(target - _val);
      }
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n) 
      {
        return ( _val > target ) ? _grad[_n]: -_grad[_n];
      }
      bool is_valid() const { return true; }
      virtual std::string print() const { return ""; }
      virtual std::string what_is() const
      { 
        std::ostringstream sstr;
        sstr << " Target (" << target << ")";
        return sstr.str();
      }
  };
  template< class T_GA_TRAITS >
  class ConvexHull : public Types< T_GA_TRAITS > :: Scalar
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef typename Types<t_GATraits> :: Scalar       t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_Base :: t_Fitness                t_Fitness;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef opt::ConvexHull::Base<t_Individual>         t_ConvexHull;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;

    protected:
      t_ConvexHull convexhull; 
      mutable bool valid;
      using t_Base :: fitness;
      using t_Base :: current_indiv;

    public:
      ConvexHull() : valid(true) {}
      ConvexHull( const ConvexHull &_c ) : convexhull( _c.convexhull ), valid(true) {}
      virtual ~ConvexHull() {}
      
      virtual const t_Fitness& operator()(const t_Quantity& _val);
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad)
        { evaluate_with_gradient( _val, _grad, _i_grad ); }

      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients& _grad,
                                                       t_VA_Type *);
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity & _val,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      bool is_valid() const
      {
        if ( valid ) return true;
        valid = true;
        return false;
      }
      virtual std::string what_is() const { return " Convex-Hull"; }
      virtual std::string print() const { return convexhull.print(); }
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
        { return convexhull.Save( _node, _op ); };
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op)
        { return convexhull.Load( _node, _op ); };
      virtual bool does_store() const { return true; }
  };
  template< class T_GA_TRAITS >
  const typename ConvexHull<T_GA_TRAITS> :: t_Fitness&
    ConvexHull<T_GA_TRAITS> :: operator()(const t_Quantity& _val)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
    
      if ( _val >= base and convexhull.size() >= 2 )
       { fitness = _val - base;  return fitness; }
    
      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      fitness = 0.0;
      return fitness;
    }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_ScalarQuantity
    ConvexHull<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                       t_QuantityGradients &_grad,
                                                       t_VA_Type *_i_grad)
    {
      t_Quantity x = current_indiv->get_concentration();
      t_Quantity base = (t_Quantity) convexhull.evaluate( x );
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      typename t_QuantityGradients :: iterator i_grad = _grad.begin(); 
      typename t_QuantityGradients :: iterator i_grad_end = _grad.end(); 
      t_VA_Type *i_grad_result = _i_grad;
      for(; i_grad != i_grad_end; ++i_grad, ++i_grad_result )
        *i_grad_result += *i_grad - gradient;
      
      if ( _val >= base ) return _val - base;

      if ( convexhull.add( _val, *current_indiv ) )
        valid = false;

      return 0.0;
    }
  template< class T_GA_TRAITS >
  typename ConvexHull<T_GA_TRAITS> :: t_VA_Type
    ConvexHull<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity &_val,
                                                      t_QuantityGradients &_grad,
                                                      types::t_unsigned _n) 
    {
      t_Quantity x = current_indiv->get_concentration();
      types::t_real Ninv = 1.0 / ( (types::t_real ) current_indiv->Object().Container().size() );
      types::t_real gradient = convexhull.evaluate_gradient( x ) * Ninv;
      return _grad[_n] - gradient;
    }


  template<class T_GA_TRAITS >
  class Container : public Types<T_GA_TRAITS> :: Vector
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef Types<t_GATraits>                          t_ObjectiveType;
      typedef typename t_ObjectiveType :: Vector          t_Base;
      typedef typename t_GATraits :: t_Individual        t_Individual;
      typedef typename t_GATraits :: t_IndivTraits       t_IndivTraits;
      typedef typename t_Base :: t_Quantity               t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity         t_ScalarQuantity;
      typedef typename t_ObjectiveType :: Scalar          t_Objective;
      typedef std::vector< t_Objective* >                 t_Objectives;
      typedef typename t_Base :: t_VA_Traits              t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>                 t_SaveOp;
      typedef GA::LoadObject<t_GATraits>                 t_LoadOp;

    protected:
      using t_Base :: fitness;
      t_Objectives objectives;

    public:
      Container() {}
      Container( const Container &_c ) : objectives(_c.objectives) {}
      virtual ~Container () 
      {
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.end();
        for(; i_objective != i_end; ++i_objective ) delete *i_objective;
        objectives.clear();
      }

      bool is_valid() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->is_valid() ) return false;
        return true;
      }
      virtual bool Save( TiXmlElement &_node, t_SaveOp& _op)
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->Save( _node, _op ) ) return false;
        return true;
      }
      virtual bool Restart( const  TiXmlElement &_node, t_LoadOp &_op) 
      {
        typename t_Objectives::iterator i_objective = objectives.begin();
        typename t_Objectives::iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( not (*i_objective)->Restart( _node, _op ) ) return false;
        return true;
      }
      virtual bool does_store() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        for(; i_objective != i_end; ++i_objective )
          if ( (*i_objective)->does_store() ) return true;
        return false;
      }
      virtual std::string print() const
      {
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        std::ostringstream sstr;
        for(; i_objective != i_end; ++i_objective )
          sstr << (*i_objective)->print();
        return sstr.str();
      }
  };

  template<class T_GA_TRAITS >
  class LinearSum : public Container<T_GA_TRAITS>
  {
    public:
      typedef T_GA_TRAITS t_GATraits;
    protected:
      typedef Types<t_GATraits>                      t_ObjectiveType;
      typedef Container<t_GATraits>                  t_Base;
      typedef typename t_Base :: t_Fitness           t_Fitness;
      typedef typename t_Base :: t_Individual        t_Individual;
      typedef typename t_Base :: t_QuantityTraits    t_QuantityTraits;
      typedef typename t_Base :: t_Quantity          t_Quantity;
      typedef typename t_Base :: t_ScalarQuantity    t_ScalarQuantity;
      typedef typename t_Base :: t_Objective         t_Objective;
      typedef std::vector< t_Objective* >            t_Objectives;
      typedef typename t_Base :: t_VA_Traits         t_VA_Traits;
      typedef typename t_VA_Traits :: t_QuantityGradients t_QuantityGradients;
      typedef typename t_VA_Traits :: t_Type              t_VA_Type;
      typedef GA::SaveObject<t_GATraits>             t_SaveOp;
      typedef GA::LoadObject<t_GATraits>             t_LoadOp;

    protected:
      using t_Base :: fitness;
      using t_Base :: objectives;
      std::vector< t_ScalarQuantity > coefs;

    public:
      LinearSum() {}
      LinearSum ( const LinearSum &_c ) : t_Base(_c), coefs(_c.coefs) {}
      virtual ~LinearSum() {}

      void add( t_Objective *_objective, t_ScalarQuantity _coef )
      {
        if ( not _objective ) return;
        coefs.push_back( _coef );
        objectives.push_back( _objective );
      }
      virtual const t_Fitness& operator()(const t_Quantity& _val);
      virtual t_ScalarQuantity evaluate_with_gradient( const t_Quantity &,
                                                       t_QuantityGradients&,
                                                       t_VA_Type *);
      virtual void evaluate_gradient( const t_Quantity &_val,
                                      t_QuantityGradients &_grad,
                                      t_VA_Type *_i_grad);
      virtual t_VA_Type evaluate_one_gradient( const t_Quantity &,
                                               t_QuantityGradients& _grad,
                                               types::t_unsigned _n);
      virtual std::string what_is() const
      {
        std::ostringstream sstr;
        sstr << "LinearSum begin{ ";
        typename t_Objectives::const_iterator i_objective = objectives.begin();
        typename t_Objectives::const_iterator i_end = objectives.begin();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        for(; i_objective != i_end; ++i_objective, ++i_coef )
          sstr << (*i_objective)->what_is() << "[" << *i_coef << "] ";
        sstr << "} end"; 
        return  sstr.str();
      }
     
  };
  
  template< class T_GA_TRAITS >
    const typename LinearSum<T_GA_TRAITS>::t_Fitness&
      LinearSum<T_GA_TRAITS> :: operator()( const t_Quantity& _val ) 
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );

        t_ScalarQuantity inter = 0;
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();
        fitness.clear();
        for(; i_objective != i_end and i_val != i_val_end; ++i_objective, ++i_coef, ++i_val )
        {
          double r = (*i_objective)->operator()( *i_val );
          fitness.push_back( r );
          inter += ( *i_coef ) * r;
        }
          
        fitness = inter;
        return fitness;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_ScalarQuantity
      LinearSum<T_GA_TRAITS> :: evaluate_with_gradient( const t_Quantity &_val,
                                                                    t_QuantityGradients &_grad,
                                                                    t_VA_Type *_i_grad)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        t_ScalarQuantity results = 0.0;
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type *i2 = _i_grad + _val.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _val.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _val.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          results +=   (*i_objective)->evaluate_with_gradient( *i_val, *i_grad, const_grad_result )
                     * (*i_coef);
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
          
        return results;
      };
  template< class T_GA_TRAITS >
    void
      LinearSum<T_GA_TRAITS> ::  evaluate_gradient( const t_Quantity &_val,
                                                                t_QuantityGradients &_grad,
                                                                t_VA_Type *_i_grad)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type *i2 = _i_grad + _val.size();
        t_VA_Type *i1 = _i_grad;
        for( ; i1 != i2; ++i1 ) *i1 = t_VA_Type(0);

        t_VA_Type *const_grad_result = new t_VA_Type[ _val.size() ];
        t_VA_Type *const_grad_result_end = const_grad_result + _val.size(); 

        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
        {
          i1 = const_grad_result;

          for( ; i1 != const_grad_result_end; ++i1 ) *i1 = t_VA_Type(0);

          (*i_objective)->evaluate_gradient( *i_val, *i_grad, const_grad_result );
          i1 = const_grad_result;

          for(; i1 != const_grad_result_end; ++i1, ++i2, ++i_coef ) *i2 +=  (*i_coef) * (*i1);
        }

        delete[] const_grad_result;
      };
  template< class T_GA_TRAITS >
    typename LinearSum<T_GA_TRAITS>::t_VA_Type
      LinearSum<T_GA_TRAITS> :: evaluate_one_gradient( const t_Quantity & _val,
                                                                   t_QuantityGradients& _grad,
                                                                   types::t_unsigned _n)
      {
        if ( t_QuantityTraits::size(_val) != coefs.size() )
          throw std::runtime_error( "Wrong number of objective functions\n" );
        typename t_Quantity :: const_iterator i_val = _val.begin();
        typename t_Quantity :: const_iterator i_val_end = _val.end();
        typename t_QuantityGradients :: iterator i_grad = _grad.begin();
        typename t_QuantityGradients :: iterator i_grad_end = _grad.end();
        typename std::vector< t_ScalarQuantity > :: const_iterator i_coef = coefs.begin();
        typename t_Objectives :: iterator i_objective = objectives.begin();
        typename t_Objectives :: iterator i_end = objectives.begin();

        t_VA_Type result = t_VA_Type(0);
        for(; i_objective != i_end and i_val != i_val_end;
              ++i_objective, ++i_coef, ++i_val, ++i_grad )
          result +=   ( *i_coef ) 
                    * (*i_objective)->evaluate_one_gradient( *i_val, *i_grad, _n );

        return result;
      };

  template< class T_GA_TRAITS >
     typename Types<T_GA_TRAITS> :: Scalar*
      Types<T_GA_TRAITS> :: scalar_from_xml( const TiXmlElement &_node )
      {
        if ( not &_node ) return NULL;
        std::string str = "minimize"; 
        std::string name = Print::lowercase(_node.Value());
        if (    name.compare("objective") == 0 
             or name.compare("method") == 0 )
        {
          if ( _node.Attribute( "type" ) )
            str = Print::lowercase(_node.Attribute( "type" ));
        }
        else if ( _node.Attribute("objective") )
          str = _node.Attribute( "objective" );
        if ( str.compare("convexhull") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: ConvexHull" << Print::endl;
          return new ConvexHull<t_GATraits>;
        }
        else if ( str.compare("minimize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Minimize" << Print::endl;
          return new Minimize<t_GATraits>;
        }
        else if ( str.compare("maximize") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "Objective: Maximize" << Print::endl;
          return new Maximize<t_GATraits>;
        }
        else if (str.compare("target") == 0 )
        {
          if( _node.Attribute("target") )
          {
            double d; _node.Attribute("target", &d );
            Print::xmg << Print::Xmg::comment
                       << "Objective: Target (" << d << ")" << Print::endl;
            return new Target<t_GATraits>( (types::t_real) d );
          }
        }
        if ( _node.FirstChildElement( "Objective" ) )
         return scalar_from_xml( *_node.FirstChildElement( "Objective" ) ); 

        return NULL;
      }
    template< class T_GA_TRAITS >
       typename Types<T_GA_TRAITS> :: Vector*
        Types<T_GA_TRAITS> :: vector_from_xml( const TiXmlElement &_node )
        {
          if ( not &_node ) return NULL;
          std::string str = "minimize"; 
          std::string name = Print::lowercase(_node.Value());
          if (    name.compare("objective") == 0 
               or name.compare("method") == 0 )
          {
            if ( _node.Attribute( "type" ) )
              str = Print::lowercase(_node.Attribute( "type" ));
          }
          if ( Vector::t_QuantityTraits::is_vector ) // and str.compare("LinearSum") == 0 )
          {
            LinearSum<T_GA_TRAITS> *linear = new LinearSum<T_GA_TRAITS>;
            if ( not linear ) 
            {
              std::cerr << "Mememory Pb when creating LinearSum multi-objective" << std::endl;
              return NULL;
            }
            const TiXmlElement *child = _node.FirstChildElement("Objective");
            for(; child; child = child->NextSiblingElement("Objective") )
            {
              Scalar* scalar = scalar_from_xml( *child );
              if ( not scalar ) continue;
              double d = 0.0;
              if ( not child->Attribute("coef", &d) ) d = 1.0;
              linear->add( scalar, t_ScalarQuantity(d) );
            }
            return linear;
          }

          if ( _node.FirstChildElement( "Objective" ) )
           return new_from_xml( *_node.FirstChildElement( "Objective" ) ); 

          return NULL;
        }

   
    // bullshit class since there is no such thing as partial specialization of functions!!
    template< class T_TYPE >
      struct fork<T_TYPE, true>
      {
        typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
        {
          typename T_TYPE :: Vector* result;
          result = T_TYPE::vector_from_xml( _node );
          if( not result ) std::cerr << "Could not create multi-objective..." << std::endl;
          return result;
        }
      };
    template< class T_TYPE > 
      struct fork<T_TYPE, false>
      {
        typename T_TYPE :: Vector* operator()( const TiXmlElement &_node )
        {
          return T_TYPE::scalar_from_xml( _node );
        }
      };

  }  // namespace Objective
 /* @} */
#endif //  _MULTIOB_OBJECTIVE_H_

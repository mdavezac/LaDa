//
//  Version: $Id$
//
#ifndef _FITNESS_H_
#define _FITNESS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<vector>
#include<math.h>
#include<stdexcept>

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#include <mpi/mpi_object.h>
#include <print/xmg.h>
#include <print/manip.h>

#include <opt/traits.h>

/** \ingroup Genetic 
 * @{*/
//! \brief Defines ordering of individuals
//! \details These class stores the so-called "polished" fitness with which
//! individuals are rated. Indeed, each individual will contain an instance of a
//! fitness (see Traits::Indiv). A fitness class will necessarily define weak ordering
//! operators, Below, two kinds of fitnesses are implemented. A scalar fitness is
//! defined which simply holds a scalar value such as a double or an
//! integer. Around this value are accreted Load and Save to XML members, an mpi
//! broadcasting member. 
//! 
//! The multi-objective version is very similar. It defines the same functions
//! around a collection  of scalar value. It uses Pareto ordering. Note however
//! that it is defined as a derived class of the scalar fitness. As such, it is
//! possible (and sometimes interesting) to define a scalar fitness from the
//! vectorial fitnesses. For an example, see Pareto ranking with niching in
//! namespace Scaling.
//! 
//! A third class, Fitness::Types defines two typedefs, Fitness::Types::Scalar
//! and Fitness::Types::Vector. Much like Objective::Types, these two types are
//! equivalent in the case of single-objective %GA. You should use
//! Fitness::Types to obtain a fitness, rather than Fitness::Scalar or
//! Fitness::Vectorial directly.
namespace Fitness
{
  //! \brief The return value of Fitness::Scalar::compare() and derived
  //! \details Comparison depends on Fuzzy operators.
  //!          In the case of <strong> vectorial /<strong> fitnesses,
  //!          \e a fitness a is \< than a fitness \e b, if and only if all
  //!          its components but one are at least \<=  than the components
  //!          of \e b, and if at least one component is stryctly \<.
  //!          Two </strong> vectorial /<strong> fitnesses are INDIFFERENT
  //!          if they neither equal, nor stronger, nor weaker than one another.
  //! \sa Ranking::Pareto
  enum t_Comparison 
  {
    //! Fitnesses are equal
    EQUAL, 
    //! Fitnesses is weaker (see t_Comparison description )
    WEAKER, 
    //! Fitnesses is stronger (see t_Comparison description )
    STRONGER, 
    //! Fitnesses are not comparable (see t_Comparison description )
    INDIFFERENT
  };

  //! \brief %Fitness base class for <em>scalar</em> fitnesses
  //! \details Defines a fitness as a variable Fitness::Base::quantity around
  //! which are implemented weak ordering operators. This class also adds
  //! loading and saving capabilities and such niceties. We first declare the
  //! template class, with the argument below, and then specialize a scalar and
  //! vectorial fitness base class.
  //!
  //! This class also (re)implements the validity check of EO. More
  //! specifically, Fitness::Base contains a logical member,
  //! Fitness::Base::is_valid which is set to true when the quantity is set. The
  //! point is to make sure that the fitness is valid prior to using it. For this
  //! purpose, Fitness::Base::invalidate should be called whenever the fitness
  //! becomes invalid, for instance after a GA::Krossover operator is applied.
  //! 
  //! \param T_QUANTITYTRAITS traits of the quantity. Generally  a
  //!        Traits::Quantity of some kind
  //! \warning using this class with a vectorial quantity is not recommended.
  //!          Use Fitness::Vectorial instead.
  template<class T_QUANTITYTRAITS >
  class Scalar
  {
    typedef Scalar<T_QUANTITYTRAITS> t_This; //!< Type of this class

    public:
      //! \brief The traits of the fitness quantity \see Traits::Quantity
      typedef T_QUANTITYTRAITS t_QuantityTraits;
      //! \brief The type of the fitness quantity
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      //! \brief The scalar fitness type
      //! \details Used by some to access strictly scalar fitnesses. If the
      //! fitness is already a scalar than t_ScalarFitness is equivalent to the
      //! fitness itself.
      typedef t_This t_ScalarFitness;

    protected:
      t_Quantity quantity; //!< quantity against which will be judged
      bool is_valid; //!< True if Fitness::Scalar::quantity has been set

    public:
      //! Constructor
      Scalar() : is_valid( false )  {}
      //! Copy Constructor
      explicit Scalar   ( const t_This & _c )
                      : quantity( _c.quantity ), 
                        is_valid( _c.is_valid ) {}
      //! Constructor and Initializer
      void operator=( const t_Quantity _fit ) 
        { quantity = _fit; is_valid = true; }
      //! Destructor
      ~Scalar() {}


      //! \brief strict ordering operator 
      //! \details Calls a static function of Fitness::Scalar::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise). Note that since minimization is the
      //! default, the implementation calls t_QuantityTraits::greater().
      //! \see Traits::Quantity, Fuzzy
      bool operator<(const t_This & _f) const
        { return t_QuantityTraits::gt(quantity, _f.quantity); }
      //! \brief strict ordering operator 
      //! \details Calls a static function of Fitness::Scalar::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise). Note that since minimization is the
      //! default, the implementation calls t_QuantityTraits::le().
      //! \see Traits::Quantity, Fuzzy
      bool operator>(const t_This & _f) const
        { return t_QuantityTraits::le(quantity, _f.quantity); }
      //! \brief equality operator 
      //! \details Calls a static function of Fitness::Scalar::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise).
      //! \see Traits::Quantity, Fuzzy
      bool operator==(const t_This & _f) const
        { return t_QuantityTraits::eq(quantity, _f.quantity); }

      //! \brief returns true if the fitness is not valid
      bool invalid() const { return not ( (bool) is_valid ); }
      //! \brief invalidates the fitness
      void invalidate() { is_valid = false; }
      //! \brief returns a constant copy of the quantity
      operator t_Quantity() const
      { 
        __ASSERT( not is_valid, "Trying to access invalid fitness\n" )
        return quantity;
      }
      //! \brief Loads the fitenss from XML
      bool Load( const TiXmlElement & _node );
      //! \brief Saves the fitenss to XML
      bool Save( TiXmlElement & _node ) const;
      //! \brief Returns a t_Comparison giving the relationship between this
      //!        object and \a _f.
      t_Comparison compare( const t_This &_f ) const;

#ifdef _MPI
      /** \ingroup MPI
       * \brief allows the serialization of a Fitness::Scalar object.
       * \details serializes the object completely, eg both quantity and
       * validity are set.
       */
      bool broadcast( mpi::BroadCast &_bc );
#endif 
  };


  /** \brief %Base class for <em>vectorial</em> fitnesses
      \details This is the multi-objective implementation of fitness. Note that
             it is derived from a scalar fitness. Indeed, a multi-objective
             approach may at some point appreciate one individual against
             another using, say, Pareto Ranking, and hence a vectorial
             quantity, or through a single scalar number, say when using a
             linear sum of objectives. By deriving the multi-objective
             fitness from the single-objective fitness, we allow both
             approach with the same object.  Individuals can be orderer
             according to the Pareto ordering operators.
             
             As defined by the Individual::Base::Fitness type, EO only knows
             about scalar fitnesses. In other words, when EO compares two
             individuals, it will compare the \e scalar fitnesses. These have
             better be set to some logical value prior to, say, parent
             selection, or replacement. Basically, all things \e vectorial are
             used and known only by LaDa. The end-product fitness is and should
             be the \e scalar fitness.  Note however that it should not be very
             difficult to change this behavior (say add a t_EO_Fitness in
             Traits::Indiv, and set Individual::Base::Fitness to t_EO_Fitness ).

             Let \f$\mathcal{F}^{v}_t(\sigma)\f$ denote the component
             \f$t\in[0,N[\f$ of <I>N</I>-dimensional fitness 
             \f$\mathcal{F}^{v}(\sigma)\f$ of an individual \f$\sigma\f$.
             Also, let \f$\mathcal{F}(\sigma)\f$ be the \b scalar end-product
             fitness of an individual \f$\sigma\f$.  An individual
             \f$\sigma_i\f$ dominates another individual \f$\sigma\f$ if
             and only if 
             \f[
                 \mathcal{F}^v(\sigma_i) \succeq \mathcal{F}^v(\sigma_j)\quad
                 \Leftrightarrow\quad \forall t\in[0,N[,\ \mathcal{F}^v_t(\sigma_i)
                 \geq \mathcal{F}^v_t(\sigma_j).
             \f]
      
  */
  template<class T_QUANTITYTRAITS >
  class Vectorial :
        public Scalar< typename T_QUANTITYTRAITS::t_ScalarQuantityTraits >
  {
    typedef Vectorial<T_QUANTITYTRAITS> t_This;      //!< Type of this class
    //! Type of the base class
    typedef Scalar<typename T_QUANTITYTRAITS::t_ScalarQuantityTraits> t_Base; 

    public:
      //! \brief The traits of the fitness quantity \see Traits::Quantity
      typedef T_QUANTITYTRAITS t_QuantityTraits;
      //! \brief The type of the fitness quantity
      typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
      //! \brief The scalar quantity traits
      typedef typename t_QuantityTraits :: t_ScalarQuantityTraits t_ScalarQuantityTraits;
      //! \brief The scalar fitness type
      //! \details Used by some to access strictly scalar fitnesses. If the
      //! fitness is already a scalar than t_ScalarFitness is equivalent to the
      //! fitness itself.
      typedef t_Base t_ScalarFitness;


    protected:
      //! \brief True if vectorial quantity above is set
      //! \details Only applies to the vectorial Fitness::Scalar::quantity, and
      //! not to the scalar t_Base::quantity. One may be set when the other is
      //! not. Indeed, in some schemes, the scalar quantity may  never be set.
      bool vec_is_valid;
      //! \brief Quantities against which all may be judged
      //! \details This is a <strong>vectorial</strong> quantity.
      //! t_Base::quantity is a <strong>scalar</strong> quantity.
      t_Quantity vec_quantity;

    public:
      //! Constructor
      Vectorial() : vec_is_valid( false )  {}
      //! Copy Constructor
      explicit Vectorial   ( const t_This & _c )
                         : t_Base(_c), vec_is_valid( _c.vec_is_valid ),
                           vec_quantity( _c.vec_quantity ) {}
      //! Constructor and Initializer
      void operator=( const t_Quantity &_fit )
        { vec_is_valid = true; vec_quantity = _fit; }
      //! \brief Copy constructor which only sets the scalar fitness.
      //! \details Nothing is done about member variables of this class. This
      //!          is meant to be a copy operator for the scalar fitness only.
      void operator=( const t_ScalarFitness &_fit ) 
        { t_Base::operator=(_fit); }
      //! \brief Copy constructor which only sets the scalar fitness.
      //! \details Nothing is done about member variables of this class. This
      //!          is meant to be a copy operator for the scalar fitness only.
      void operator=( const typename t_ScalarFitness :: t_Quantity &_fit )
        { t_Base::operator=(_fit); }
      //! Destructor
      ~Vectorial() {}


      /** \brief Pareto ordering, \f$\mathcal{F}^v(\sigma_i) \preceq \mathcal{F}^v(\sigma_j)\f$
          \details Implements the relationship
             \f[
                 \mathcal{F}^v(\sigma_i) \preceq \mathcal{F}^v(\sigma_j)\quad
                 \Leftrightarrow\quad \forall t\in[0,N[,\ \mathcal{F}^v_t(\sigma_i)
                 \geq \mathcal{F}^v_t(\sigma_j).
             \f]
                  For generality, the static ordering operators of
                  t_ScalarQuantityTraits are use in the implementation (eg
                  fuzzy math for reals). Note however that it is Fuzzy::gt
                  which is called, since minimization is the default.
      */
      bool operator<(const t_This & _f) const;
      /** \brief Pareto ordering, \f$\mathcal{F}^v(\sigma_i) \succeq \mathcal{F}^v(\sigma_j)\f$
          \details Implements the relationship
             \f[
                 \mathcal{F}^v(\sigma_i) \succeq \mathcal{F}^v(\sigma_j)\quad
                 \Leftrightarrow\quad \forall t\in[0,N[,\ \mathcal{F}^v_t(\sigma_i)
                 \geq \mathcal{F}^v_t(\sigma_j).
             \f]
                  For generality, the static ordering operators of
                  t_ScalarQuantityTraits are use in the implementation (eg
                  fuzzy math for reals). Note however that it is Fuzzy::gt
                  which is called, since minimization is the default.
      */
      bool operator>(const t_This & _f) const;
      //! \brief Pareto ordering, eg true if both fitnesses are equally
      //! dominant and dominated
      bool operator==(const t_This & _f) const;

      //! \brief returns true if the (vectorial) fitness is not valid
      bool invalid() const { return not ( (bool) vec_is_valid ); }
      //! \brief invalidates the (vectorial) fitness
      void invalidate() { vec_is_valid = false; }
      //! \brief returns a constant reference of the vectorial quantity
      operator const t_Quantity& () const
      { 
        __ASSERT( not vec_is_valid, "Trying to access invalid fitness\n" )
        return vec_quantity;
      }
      //! \brief Loads the fitenss from XML
      bool Load( const TiXmlElement & _node );
      //! \brief Saves the fitenss to XML
      bool Save( TiXmlElement & _node ) const;

      //! \brief clears the quantity
      void clear() { vec_quantity.clear(); vec_is_valid = false; }
      //! \brief adds a component to the vectorial quantity
      void push_back( typename t_Quantity :: value_type  _var )
        { vec_quantity.push_back( _var ); vec_is_valid = true; }
      //! \brief Returns a t_Comparison giving the relationship between this
      //!        object and \a _f.
      t_Comparison compare( const t_This &_f ) const;
#ifdef _MPI
      /** \ingroup MPI
       * \brief allows the serialization of a Fitness::Vectorial object.
       * \details serializes the object completely, eg both quantity and
       * validity are set for both t_This and t_Base.
       */
      bool broadcast( mpi::BroadCast &_bc );
#endif 
  };

  //! \brief Template type returning scalar fitnesses on demand
  //! \details By default, the type is a non-scalar fitness
  template< class T_FITNESS, bool indeed = false >
  struct GetScalar{ 
    typedef T_FITNESS t_result; //!< resulting type
  };
  //! Specialization of GetScalar for a true scalar fitness type
  template< class T_FITNESS >
  struct GetScalar<T_FITNESS, true>
  {
    typedef typename T_FITNESS :: t_ScalarFitness t_result; //!< resulting type
  };



  //! \brief Helper class for defining fitnesses.
  //! \details Two fitness types are defined automatically: Types::Scalar and Types::Vector.
  //!          In the case of a multi-objective %GA (vector quantity), each is
  //!          what you would expect. In the case of single-objective %GA
  //!          (scalar quantity), both Types::Vector and Types::Scalar are
  //!          equivalent and simply typedef a Fitness::Scalar. Generally, you
  //!          will only want to use Types::Vector.
  template<class T_QUANTITYTRAITS,
           bool IS_SCALAR = T_QUANTITYTRAITS::is_scalar >
  struct Types 
  {
    typedef Vectorial<T_QUANTITYTRAITS> Vector; //!< A \e truly vectorial fitness type
    //! A scalar fitness type
    typedef typename Vectorial<T_QUANTITYTRAITS> :: t_ScalarFitness Scalar;
  };
  //! %Scalar flavor 
  template<class T_QUANTITYTRAITS >
  struct Types <T_QUANTITYTRAITS, true >
  {
    typedef Scalar<T_QUANTITYTRAITS> Vector; //!< A scalar fitness type (sic)
    typedef Scalar<T_QUANTITYTRAITS> Scalar; //!< A scalar fitness type
  };


  //! Dumps fitness to a stream
  template<class T_QUANTITYTRAITS>
  std::istream & operator>>( std::istream &_is,
                             Scalar<T_QUANTITYTRAITS> &_fit );
  
  //! Retrieves fitness from a stream
  template<class T_QUANTITYTRAITS>
  std::ostream & operator<<( std::ostream &_os,
                             const Scalar<T_QUANTITYTRAITS> &_fit );


  //! Dumps \e vectorial fitness to a stream
  template<class T_QUANTITYTRAITS>
    std::istream & operator>>( std::istream &_is,
                               Vectorial<T_QUANTITYTRAITS> &_fit );

  //! Retrieves \e vectorial fitness from a stream
  template<class T_QUANTITYTRAITS>
    std::ostream & operator<<( std::ostream &_os,
                               const Vectorial<T_QUANTITYTRAITS> &_fit );
}

#include "fitness.impl.h"

/*@}*/

#endif // _FITNESS_H_

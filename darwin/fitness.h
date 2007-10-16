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
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif
#include <print/xmg.h>
#include <print/manip.h>

#include<opt/traits.h>

/** \ingroup Genetic 
 * @{*/
//! \brief Defines ordering of individuals
//! \details This class stores the so-called "polished" fitness with which
//! individuals are rated. Indee, each individual will contain an instance of a
//! fitness. A fitness class will necessarily define weak ordering
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
namespace Fitness
{
  //! \brief %Base class for fitnesses
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
  //! purpose, Fitness::Base::invalidate should be called whenever the fintess
  //! becomes invalid, for instance after a GA::Krossover operator is applied.
  //! 
  //! \param T_QUANTITYTRAITS traits of the quantity. Generally  a
  //!        Traits::Quantity of some kind
  //! \param IS_SCALAR should be set to true for scalar fitnesses. Mainly, this
  //! allows us to define a multi-objective fitness as a derived type, capable of
  //! being both a  vectorial quantity (say for Pareto ranking) and a scalar
  //! quantity (say for a linear sum of objectives).
  template<class T_QUANTITYTRAITS,
           bool IS_SCALAR = T_QUANTITYTRAITS::is_scalar >
  class Base {};

  //! %Ftiness %base class for <em>scalar</em> fitnesses
  template<class T_QUANTITYTRAITS >
  class Base<T_QUANTITYTRAITS, true>
  {
    typedef Base<T_QUANTITYTRAITS, true> t_This; //!< Type of this class

    //! \brief Dumps fitness to a stream
    template<class TQUANTITYTRAITS>
      friend std::istream & operator>>( std::istream &_is,
                                        const Base<TQUANTITYTRAITS, true> &_fit );
    
    //! \brief Retrieves fitness from a stream
    template<class TQUANTITYTRAITS>
      friend std::ostream & operator<<( std::ostream &_os,
                                        Base<TQUANTITYTRAITS, true> &_fit );
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
      bool is_valid; //!< True if Fitness::Base::quantity has been set

    public:
      //! Constructor
      Base() : is_valid( false )  {}
      //! Copy Constructor
      Base( const Base & _c ) : quantity( _c.quantity ), is_valid( _c.is_valid ) {}
      //! Constructor and Initializer
      Base( const types::t_real _fit ) : quantity( _fit ), is_valid( true ) {}
      //! Destructor
      ~Base() {}


      //! \brief strict ordering operator 
      //! \details Calls a static function of Fitness::Base::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise).
      //! \see Traits::Quantity, Traits::Fuzzy
      bool operator<(const Base & _f) const
        { return t_QuantityTraits::less(quantity, _f.quantity); }
      //! \brief strict ordering operator 
      //! \details Calls a static function of Fitness::Base::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise).
      //! \see Traits::Quantity, Traits::Fuzzy
      bool operator>(const Base & _f) const
        { return t_QuantityTraits::greater(quantity, _f.quantity); }
      //! \brief equality operator 
      //! \details Calls a static function of Fitness::Base::t_QuantityTraits.
      //! This allows type specific implementation, such as fuzzy math for
      //! reals (to avoid numerical noise).
      //! \see Traits::Quantity, Traits::Fuzzy
      bool operator==(const Base & _f) const
        { return t_QuantityTraits::equal(quantity, _f.quantity); }

      //! \brief returns true if the fitness is not valid
      bool invalid() const { return not is_valid; }
      //! \brief invalidates the fitness
      void invalidate() { is_valid = false; }
      //! \brief returns a constant copy of the quantity
      operator t_Quantity() const;
      //! \brief Loads the fitenss from XML
      bool Load( const TiXmlElement & _node );
      //! \brief Saves the fitenss to XML
      bool Save( TiXmlElement & _node ) const;

#ifdef _MPI
      /** \ingroup MPI
       * \brief allows the serialization of a Fitness::Base object.
       * \details serializes the object completely, eg both quantity and
       * validity are set.
       */
      bool broadcast( mpi::BroadCast &_bc );
#endif 
  };




  //! \brief %Base class for <em>vectorial</em> fitnesses
  //! \details This is the multi-objective implementation of fitness. Note that
  //! it is derived from a scalar fitness. Indeed, a multi-objective
  //! approach may at some point appreciate one individual against another using,
  //! say, Pareto Ranking, and hence a vectorial quantity, or through a single
  //! scalar number, say when using a linear sum of objectives. By deriving the
  //! multi-objective fitness from the single-objective fitness, we allow both
  //! approach with the same object.
  template<class T_QUANTITYTRAITS >
  class Base<T_QUANTITYTRAITS, false> :
        public Base< typename T_QUANTITYTRAITS::t_ScalarQuantityTraits, true >
  {
    typedef Base<T_QUANTITYTRAITS, false> t_This;      //!< Type of this class
    //! Type of the base class
    typedef Base<typename T_QUANTITYTRAITS::t_ScalarQuantityTraits, true> t_Base; 

    //! \brief Dumps fitness to a stream
    template<class TQUANTITYTRAITS>
      friend std::istream & operator>>( std::istream &_is,
                                        Base<TQUANTITYTRAITS, false> &_fit );

    //! \brief Retrieves fitness from a stream
    template<class TQUANTITYTRAITS>
      friend std::ostream & operator<<( std::ostream &_os,
                                        const Base<TQUANTITYTRAITS, false> &_fit );

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
      //! \brief Quantities against which all may be judged
      //! \details This is a <strong>vectorial</strong> quantity.
      //! t_Base::quantity is a <strong>scalar</strong> quantity.
      t_Quantity quantity;
      //! \brief True if vectorial quantity above is set
      //! \details Only applies to the vectorial Fitness::Base::quantity, and
      //! not to the scalar t_Base::quantity. One may be set when the other is
      //! not. Indeed, in some schemes, the scalar quantity may  never be set.
      bool is_valid;

    public:
      //! Constructor
      Base() : is_valid( false )  {}
      //! Copy Constructor
      Base( const Base & _c ) : t_Base(_c), quantity( _c.quantity ), is_valid( _c.is_valid ) {}
      //! Constructor and Initializer
      Base( const t_Quantity &_fit ) : t_Base(), quantity( _fit ), is_valid( true ) {}
      //! Destructor
      ~Base() {}


      //! \brief initializes t_Base with a scalar fitness.
      void operator=( const t_ScalarFitness &_fit ) { (t_Base&) *this = _fit; }
      //! \brief initializes t_Base with a scalar quantity
      void operator=( const typename t_ScalarFitness::t_Quantity &_fit )
        { (t_Base&) *this = _fit; }

      //! \brief Pareto ordering 
      bool operator<(const Base & _f) const;
      //! \brief Pareto ordering 
      bool operator>(const Base & _f) const;
      //! \brief Pareto ordering, eg true if both fitnesses are equally
      //! dominant and dominated
      bool operator==(const Base & _f) const;

      //! \brief returns true if the (vectorial) fitness is not valid
      bool invalid() const { return not is_valid; }
      //! \brief invalidates the (vectorial) fitness
      void invalidate() { is_valid = false; }
      //! \brief returns a constant reference of the vectorial quantity
      operator const t_Quantity& () const;
      //! \brief Loads the fitenss from XML
      bool Load( const TiXmlElement & _node );
      //! \brief Saves the fitenss to XML
      bool Save( TiXmlElement & _node ) const;

      //! \brief clears the quantity
      void clear() { quantity.clear(); }
      //! \brief adds a component to the vectorial quantity
      void push_back( typename t_Quantity :: value_type  _var )
        { quantity.push_back( _var ); }
#ifdef _MPI
      /** \ingroup MPI
       * \brief allows the serialization of a Fitness::Base object.
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
}


#include "fitness.impl.h"

/*@}*/

#endif // _FITNESS_H_

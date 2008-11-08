//
//  Version: $Id$
//
#ifndef _DARWIN_RESULTS_H_
#define _DARWIN_RESULTS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <list>
#include <algorithm>
#include <iostream>

#include <eo/eoPop.h>

#include <tinyxml/tinyxml.h>

#include <opt/function_base.h>
#include <opt/fitness_function.h>
#include <opt/convex_hull.h>
#include <print/xmg.h>
#include <mpi/mpi_object.h>

#include "objective.h"
#include "loadsave.h"
#include "gatraits.h"
#include "loadsave.h"

namespace LaDa
{
  namespace GA
  {
  /** \ingroup Genetic 
   * @{*/
  //! \brief Implements storage capacity for use with Evaluation classes
  //! \details Storage is performed <em>via</em> a simple call to
  //! Store::Base::operator(const t_Indiviudal& _indiv). Whether \a _indiv is
  //! stored will depend upon the particular class  used.
  //! Two storage classes have defined at this point:
  //! - Store::Base, an abstract base class
  //! - Store::Conditional, stores an individual depending a template T_CONDITION functor.
  //! A number of T_CONDITION have defined in namespace Store::Condition
  //! \sa Evaluation
  namespace Store
  { 
    //!  Abstract base class for results and storage
    template<class T_GATRAITS>
    class Base
    {
      public:
        typedef T_GATRAITS t_GATraits; //!< defines all GA classes \sa Traits::GA

      private:
        typedef typename t_GATraits :: t_Evaluator  t_Evaluator; //!< Evaluator type
        typedef typename t_GATraits :: t_Individual t_Individual;//!< Individual type

      protected:
        //! \brief reference to functional evaluator \sa GA::Darwin::evaluator
        //! \details used only for load and save stuff. Such functions should probably become
        //! functors to be set in Traits::GA in the future.
        t_Evaluator &evaluator;
        //! \brief for print_out purposes
        //! \details Set to true each time a new individual is stores. Should be
        //! set to false upon call to print(). It is mutable so that it can be
        //! changed upon a call to a Base::print_out(), even though print_out
        //! functions are generally decalred const.
        mutable bool new_results;

      public:
        //! Constructor and basic Initializer
        Base (t_Evaluator &_eval) : evaluator(_eval), new_results(false) {};
        //! Copy Constructor
        Base (const Base<t_GATraits> &_c) : evaluator(_c.eval), new_results(_c.new_results) {};
        //! Destructor
        virtual ~Base() {};

        //! \brief Should reload stored individuals from XML input
        //! \details Can use Base::evaluator reference for this purpose
        virtual bool Restart( const TiXmlElement &_node ) = 0;
        //! \brief Should save stored individuals to XML ouput
        //! \details Can use Base::evaluator reference for this purpose
        virtual bool Save( TiXmlElement &_node ) const = 0;

        //! \brief Stores \a _indiv 
        //! \details Can potentially store \a _indiv if \a _indiv statisfies
        //! whatever condition the base class cares about. This operator will be
        //! called by Evaluation classes.
        //! \sa Evaluation::Base::evaluate, Evaluation::WithHistory::Evaluate
        virtual void operator()( const t_Individual &_indiv ) = 0;
        //! \brief returns true if new results have been found
        bool newresults() const { return new_results; }
        //! \brief prints new results
        //! \param _age curreng GA iteration
        //! \param is_comment Whether print_out should be in Print::Xmg::comment format
        //! \sa Print::Xmg
        virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
          { new_results = false; }
        //! Returns a string defining this derived class
        virtual std::string what_is() const = 0;
        //! Other print method
        virtual std::string print() const = 0; 
        //! Applies \a _op to all stored individuals.
        virtual void apply_all( eoMonOp<const t_Individual> *_op ) const = 0; 
        //! Applies \a _op to best stored individuals.
        virtual void apply_best( eoMonOp<const t_Individual> *_op ) const = 0; 
    };

    //! \brief stores individuals depending upon the return value of a
    //! conditional functor of template T_CONDITON
    //! \details T_CONDITION must return true when Conditional should
    //! <em>not</em> store the individual under consideration. This choice may
    //! seem odd, but it works better when using std::remove_if to remove
    //! individuals which once fulfilled the condition, but for whatever reason, do
    //! not anymore.
    //! For examples of T_CONDITION, see namespace Store::Condition below.
    template<class T_CONDITION, class T_GATRAITS>
    class Conditional : public Base<T_GATRAITS>
    {
      public:
        typedef T_CONDITION t_Condition; //!< type of functor usd to determine storage
        typedef T_GATRAITS t_GATraits; //!< all GA classes \sa Traits::GA

      private:
        typedef typename t_GATraits :: t_Evaluator  t_Evaluator; //!< Evaluator type
        typedef Base<t_GATraits>                    t_Base;      //!< %Base class type
        typedef typename t_GATraits :: t_Individual t_Individual; //!< Individual type
        typedef std::list<t_Individual>             t_Container; //!< Type of storage container

      protected:
        using t_Base :: new_results; 
        using t_Base :: evaluator;
        t_Condition condition; //!< Condition which determines storage \sa Store::Condition
        //! \brief whether to print condition, results, or all in Conditional::print_results() const
        types::t_unsigned print_what; 
        const static types::t_unsigned PRINT_RESULTS = 1; //!< print_what says print results 
        const static types::t_unsigned PRINT_CONDITION = 2;  //!< print_what says print condition 
        t_Container results; //!< Container where to store "good" individuals

      public:
        //! \brief Constructor and Initializor
        //! \details Conditional::condition is initialized entirely from XML input
        //! Conditional itself only reads wether there exists a print attribute, and sets 
        //! Conditional::print_what accordingly
        Conditional   (t_Evaluator &_eval, const TiXmlElement &_node)
                    : t_Base( _eval), condition( _node ),
                      print_what(PRINT_CONDITION)  { Load( _node ); }
        //! \brief Constructor and Initializor
        //! \details Conditional::condition is initialized both from XML input and \a _type
        template< class T_TYPE >
        Conditional   (t_Evaluator &_eval, T_TYPE _type, const TiXmlElement &_node)
                    : t_Base( _eval), condition( _type, _node ), 
                      print_what( PRINT_CONDITION ) { Load( _node ); }
        //! Destructor
        virtual ~Conditional() {}

        //! \brief simply stores results for which T_CONDTION::operator()( const t_Individual& )
        //!    returns false
        virtual void operator()( const t_Individual &_indiv );

        //! \brief Reloads stored individuals from XML input
        //! \details Uses Base::evaluator reference for this purpose.
        //! \sa darwin::LoadObject
        bool Restart( const TiXmlElement &_node );
        //! \brief Saves stored individuals into XML input
        //! \details Uses Base::evaluator reference for this purpose.
        //! \sa darwin::SaveObject
        bool Save( TiXmlElement &_node ) const;

        //! \brief Prints results to Print::Xmg
        //! \sa Base::print_results, Base::new_results, GA::IslandsContinuator
        virtual void print_results(types::t_unsigned _age, bool is_comment = false) const;
        //! \brief Prints condition only
        virtual std::string print() const { return condition.print(); }
        //! \brief Returns a string characterizing Conditional (and its Conditional::condition)
        virtual std::string what_is() const;
       
        //! Applies \a _op to all stored individuals.
        virtual void apply_all( eoMonOp<const t_Individual> *_op ) const;
        //! Applies \a _op to best stored individuals.
        virtual void apply_best( eoMonOp<const t_Individual> *_op ) const
          { condition.apply2optimum(_op); }

      private:
        //! Checks from XML input what to print out in print_results 
        bool Load( const TiXmlElement &_node );
    };

    //! \brief Defines some functors for conditional storage
    //! \details Implements the functors used in storage class Store::Conditional. Only
    //! a few behaviors are require:
    //! - a constructor which takes an xml element only
    //! - a constructor which takes an xml element only and some other type
    //! - the functor member, operator()(const t_Individual&) which return false
    //! - Restart for reloading previous state in XML format
    //! - Save for saving current state in XML format
    //! - print for printing current state to a string
    //! - what_is returns a string defining this condition
    //! when an individual should <em> not </em> be stored.
    //! At this point, Three classes are defined
    //! - Condition :: BaseOptima a base class which basic implements for "best-of" behaviors
    //! - Condition :: FromObjective implements \e Best-Of behavior using Objective classes
    //! - Condition :: Optima implements \e Best-Of behavior using ordering
    //! operators acting upon indiividuals
    //! \sa Store::Conditional 
    namespace Condition
    {
      //! \brief %Base class for \e Best-Of Behavior
      //! \details Only implements BaseOptima::Restart, BaseOptima::Save,
      //! BaseOptima::print_results for the one "best" individual stored as
      //! BaseOptima::optimum. BaseOptima::optimum will be used in derived classes
      //! for comparison
      template< class T_GATRAITS >
      class BaseOptima 
      {
        public:
          typedef T_GATRAITS t_GATraits; //!< All GA classes \sa Traits::GA

        protected:
          typedef typename t_GATraits :: t_Evaluator  t_Evaluator; //!< Evaluator type
          typedef typename t_GATraits :: t_Individual t_Individual;//!< Individual type
          typedef GA::SaveObject<t_GATraits>           t_SaveOp; //!< Functor for saving individuals
          typedef GA::LoadObject<t_GATraits>           t_LoadOp; //!< Functor for loading individuals

        protected:
          t_Individual *optimum; //!< "best" individual
          //! Whether the memory pointed to by optimum is owned by this object
          bool owns_optimum;

        public:
          //! Constructor
          BaseOptima   ( const TiXmlElement &_node )
                     : optimum(NULL), owns_optimum(false) {};
          //! Copy Constructor
          BaseOptima   ( const BaseOptima &_c )
                     : optimum(_c.optimum), owns_optimum(false) {};
          //! Destructor
          ~BaseOptima();

          //! \brief Reloads stored BaseOptima::optimum from XML input
          bool Restart( const TiXmlElement &_node, t_LoadOp & _op);
          //! \brief Saves stored BaseOptima::optimum from XML input
          bool Save( TiXmlElement &_node, t_SaveOp & _op) const;
          //! \brief Returns a string characterizing BaseOptima
          std::string what_is() const { return " BaseOptima "; } 
          //! \brief prints out BaseOptima::optimum characteristics
          std::string print() const;
          //! Applies \a _op to best stored individuals.
          void apply2optimum( eoMonOp<const t_Individual> *_op ) const { (*_op)( *optimum ); }
      };

      //! \brief Implements \e Best-Of behavior using an individual's fitness and an Objective
      //! \details An individual is deemed good for storage depending on the FromObjective::objective,
      //! FromObjective::objective can be the provided as a pointer upon
      //! constructing this object, or it can be created from an XML input. From
      //! there on, and individual is stored if is better than the inherited member
      //! BaseOptima::optimum, or when it is within FromObjective::delta of
      //! BaseOptima::optimum. If is is better that BaseOptima::optimum, it
      //! replaces BaseOptima::optimum.
      template< class T_GATRAITS >
      class FromObjective : public BaseOptima<T_GATRAITS>
      {
        public:
          typedef T_GATRAITS t_GATraits; //!< All GA %types \sa Traits::GA

        protected:
          typedef typename t_GATraits :: t_Evaluator           t_Evaluator; //!< Evaluator type
          typedef BaseOptima<t_GATraits>                       t_Base;     //!< base class type
          typedef typename t_GATraits :: t_Individual          t_Individual; //!< individual type
          typedef typename t_GATraits :: t_QuantityTraits      t_QuantityTraits;//!< quantity traits
          typedef typename t_QuantityTraits :: t_Quantity      t_Quantity;  //!< quantity type
           //! scalar quantity type
          typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity; 
          typedef typename Objective::Types<t_GATraits>        t_Objective; //!< objective typedef
          typedef GA::SaveObject<t_GATraits> t_SaveOp; //!< Save functor
          typedef GA::LoadObject<t_GATraits> t_LoadOp; //!< load functor

        protected:
          using t_Base :: optimum;
          using t_Base :: owns_optimum;
          typename t_Objective :: t_Vector *objective;  //!< objective type
          t_ScalarQuantity val; //!< optimal value
          t_ScalarQuantity end_val; //!< end of storage interval
          t_ScalarQuantity delta; //!< length of storage interval
          bool owns_objective; //!< wether FomOjbective::objective pointer is owned by this object
          //! \brief True if condition should be printed, as well as objective
          //! \details. More specifically, this variable modulates the behavior
          //! of FromObjective::print so that when it is true, BaseOptima::print() is
          //! called in addition to Objective::Base::print()
          bool print_condition;

        public:
          using t_Base :: print; 
          using t_Base :: Restart; 
          using t_Base :: Save; 

        public:
          //! \brief XML constructor 
          //! \sa BaseOptima::BaseOptima(const TiXmlElement&), 
          //! Store and Store::Conditional overview
          FromObjective   ( const TiXmlElement &_node );
          //! \brief XML and objective constructor 
          //! \sa BaseOptima::BaseOptima(const TiXmlElement&), 
          //! Store and Store::Conditional overview
          FromObjective   ( typename t_Objective::t_Vector* _type, const TiXmlElement &_node );
          //! \brief Copy Constructor
          FromObjective   (const FromObjective &_c )
                        : t_Base( _c ), objective( _c.objective ),
                          val( _c.val ), end_val( _c.end_val ), delta( _c.delta ),
                          owns_objective( false ), print_condition( _c.print_condition ) {}
          //! \brief Destructor
          //! \details deleletes FromObjective::objective only if objective is
          //! non-null and if FromObjective::objective is owned according to
          //! FromObjective::owns_objective.
          ~FromObjective();

          //! \brief returns false if \a _indiv should <em>not</em> be stored
          bool operator()( const t_Individual &_indiv );

          //! \brief returns a string with stuff that FomObjective::objective
          //!        store, eg convexhull
          //! \details A bit complicated... If the objective does \e not store
          //!          anything itself (eg as returned by
          //!          Objective::Base::does_store() const), then this function returns
          //!          with a call to BaseOptima::print() const. If on the other
          //!          hand Objective::Base::does_store() const returns true,
          //!          then two further cases appear. If
          //!          FromObjective::print_condition is set, both the print
          //!          routines of the condition and objective are called.
          //!          Otherwise, only the print routine of the objective
          //!          function is called.
          std::string print() const;
          //! Return a string characterizing FromOjbective
          std::string what_is() const;
          //! \brief Reloads previously saved state from XML input
          //! \param _node XML node from which to reload
          //! \param _op Darwin::LoadObject with which to reload t_Individuals
          bool Restart( const TiXmlElement &_node, t_LoadOp & _op)
            { return t_Base::Restart( _node, _op) and objective->Restart( _node, _op); }
          //! Save current state to XML output
          //! \param _node XML node to which current state should be saved
          //! \param _op Darwin::SaveObject with which to save t_Individuals
          bool Save( TiXmlElement &_node, t_SaveOp & _op) const;
      };

      //! \brief Functor which returns true if an individual is not as good as
      //!        the current optimum
      //! \details Comparisons are done using the ordering operators defined in Individual
      //! (eg depending on individuals fitness)
      template< class T_GATRAITS >
      class Optima : public BaseOptima<T_GATRAITS>
      {
        public:
          typedef T_GATRAITS t_GATraits; //!< All GA objects \sa Traits::GA
        protected:
          typedef BaseOptima<t_GATraits>                       t_Base; //!< %Base class type
          typedef typename t_GATraits :: t_Evaluator           t_Evaluator; //!< Evaluator type
          typedef typename t_GATraits :: t_Individual          t_Individual; //!< Individual type

        protected:
          using t_Base :: optimum;

        public:
          //! \brief XML Constuctor 
          //! \sa BaseOptima::BaseOptima(const TiXmlElement), 
          //! Store and Store::Condition overviews.
          Optima ( const TiXmlElement &_node ) : t_Base(_node) {}
          //! Destructor
          ~Optima () {}

          //! \brief Returns true if _indiv is <em>not</em> better than
          //!        inherited member BaseOptima::optimum
          //! \details Set to work with std::remove_if. BaseOptima::optimum is
          //! set to _indiv if _indiv is better than BaseOptimat::optimum. By
          //! better, we mean the ordering
          //! operator(const t_Individual&, const t_Individual) 
          //! which should be defined.
          bool operator()( const t_Individual &_indiv )
          {
            if (  (not optimum) or optimum->invalid() ) return true;
            if (  _indiv > *optimum ) return true; 
            *optimum = _indiv;
            return false; 
          }
          //! Returns a string characteristic of Optima.
          std::string what_is() const { return " optimum "; }
      };
    } // namespace Condition

    //! Typedef collection of a few standard Objective %types 
    template<class T_GATRAITS>
      struct Types
      {
        //! \brief A conditional storage type which stores only best individuals
        //! \details Best is defined by  bool operator(const t_Individual&, const t_Individual) 
        typedef Store::Conditional< Store::Condition::Optima< T_GATRAITS >,
                                    T_GATRAITS > Optima;
        //! \brief A conditional storage type which stores only best individuals
        //!  according to objective
        typedef Store::Conditional< Store::Condition::FromObjective< T_GATRAITS >,
                                    T_GATRAITS >  FromObjective;

        //! \brief The type of the abstract base class. 
        typedef Store::Base< T_GATRAITS > Base;
      };
  } // namespace Store

  } //namespace GA
  /*@}*/
} // namespace LaDa

#include "store.impl.h"


#endif // _RESULTS_H_

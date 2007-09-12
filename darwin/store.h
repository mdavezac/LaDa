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

#include "opt/opt_function_base.h"
#include "opt/fitness_function.h"
#include "opt/convex_hull.h"
#include "print/xmg.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "objective.h"
#include "taboos.h"
#include "loadsave.h"
#include "gatraits.h"
#include "loadsave.h"

namespace Store
{ 
  // abstract base class for results and storage
  template<class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
  class Base
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_GA_TRAITS t_GA_Traits;

    private:
      typedef typename t_GA_Traits :: t_Individual t_Individual;

    protected:
      t_Evaluator &evaluator;
      mutable bool new_results;

    public:
      Base (t_Evaluator &_eval) : evaluator(_eval), new_results(false) {};
      Base (const Base<t_Evaluator> &_c) : evaluator(_c.eval), new_results(_c.new_results) {};
      virtual ~Base() {};

      virtual bool Restart( const TiXmlElement &_node ) = 0;
      virtual bool Save( TiXmlElement &_node ) const = 0;

      virtual void operator()( const t_Individual &_indiv ) = 0;
      bool newresults() const { return new_results; }
      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const = 0;
      virtual std::string what_is() const = 0;
      virtual std::string print() const = 0;

#ifdef _MPI
      virtual void synchronize() = 0;
#endif
  };

  // Store depends on T_CONDTION object
  // T_CONDITION must return true when Conditional should NOT store
  // For examples of T_CONDITION, see namespace Condition below
  // easier to use are "templated typedef" in struct Condition below
  template<class T_EVALUATOR, class T_CONDITION, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
  class Conditional : public Base<T_EVALUATOR, T_GA_TRAITS>
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef T_CONDITION t_Condition;
      typedef T_GA_TRAITS t_GA_Traits;

    private:
      typedef Base<t_Evaluator, t_GA_Traits> t_Base;
      typedef typename t_GA_Traits :: t_Individual t_Individual;
      typedef std::list<t_Individual> t_Container;

    protected:
      using t_Base :: new_results;
      using t_Base :: evaluator;
      t_Condition condition;
      t_Container results;
#ifdef _MPI
      t_Container new_optima;
#endif 

    public:
      Conditional   (t_Evaluator &_eval, const TiXmlElement &_node)
                  : t_Base( _eval), condition( _node ) {}
      template< class T_TYPE >
      Conditional   (t_Evaluator &_eval, T_TYPE _type, const TiXmlElement &_node)
                  : t_Base( _eval), condition( _type, _node ) {}
      virtual ~Conditional() {}

      // non-mpi version simply stores best N results 
      // MPI version is a bit of a hack
      // the object is for MPI version to be able to store in "unsynchronized" new_optima by default
      // and then into "synchronized" result upon call to synchronize
      virtual void operator()( const t_Individual &_indiv )
      {
        // if first evaluation, adds it directly
        if( condition( _indiv ) ) return;

        // checks wether result should be stored
        typename t_Container :: iterator i_found;
#ifdef _MPI
        if ( not new_optima.empty() )
        {
          i_found = std::find( new_optima.begin(), new_optima.end(), _indiv );
          if ( i_found != new_optima.end() ) return;
        }
#endif
        if ( not results.empty() )
        {
          i_found = std::find( results.begin(), results.end(), _indiv );
          if ( i_found != results.end() ) return;
        }
        new_results = true;
#ifdef _MPI
        if ( not new_optima.empty() ) new_optima.remove_if( condition );
        new_optima.push_back( _indiv );
#else
        if ( not results.empty() ) results.remove_if( condition );
        results.push_back( _indiv );
#endif
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        darwin::LoadObject<t_Evaluator> loadop( evaluator, &t_Evaluator::Load,
                                                darwin::LOADSAVE_LONG);
        if ( not condition.Restart( *xmlresults, loadop ) )
        {
          Print::xmg << Print::Xmg::comment << "Could not load condition" << Print::endl
                     << Print::Xmg::comment << "Aborting Load of Result"  << Print::endl;
          return false;
        }

        results.clear();
        darwin::LoadIndividuals( *xmlresults, loadop, results );
        // This particular line acts as Restart for condition if necessary
        std::remove_if(results.begin(), results.end(), condition); 
        Print::xmg << Print::Xmg::comment << "Reloaded Optimum and "
                   << results.size() << " target results" << Print::endl;
        print_results(0, true);
        return true;
      }
      bool Save( TiXmlElement &_node ) const
      {
        darwin::SaveObject<t_Evaluator> saveop( evaluator, &t_Evaluator::Save,
                                                darwin::LOADSAVE_LONG);
        TiXmlElement *parent = new TiXmlElement("Results");
        if ( not parent ) 
        {
          std::cerr << "Memory error while saving results" << std::endl;
          return false;
        }

        if (      condition.Save(*parent, saveop) 
             and  SaveIndividuals( *parent, saveop, results.begin(), results.end() ) ) 
        {
          _node.LinkEndChild(parent);
          return true;
        }

        std::cerr << "Could not save resuls as requested" << std::endl;
        delete parent;
        return false;
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
      {
        typename t_Container :: const_iterator i_indiv = results.begin();
        typename t_Container :: const_iterator i_end = results.end();
        for (; i_indiv != i_end; ++i_indiv )
          Print::xmg << ( is_comment ? Print::Xmg::comment: Print::Xmg::clear )
                     << std::setw(12) << std::setprecision(7)
                     << _age << " "
                     << i_indiv->get_concentration() << " "
                     << i_indiv->fitness() << Print::endl;
        new_results = false;
      }
      virtual std::string print() const { return condition.print(); }
      virtual std::string what_is() const
      { 
        std::ostringstream sstr; sstr << " Conditional on"  << condition.what_is(); 
        return sstr.str();
      }
#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not condition.broadcast(_bc) ) return false;
        types::t_int n = results.size();
        if ( not _bc.serialize(n) ) return false;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          results.resize(n);
        typename t_Container :: iterator i_res = results.begin();
        typename t_Container :: iterator i_res_end = results.end();
        for(; i_res != i_res_end; ++i_res )
          if( not i_res->broadcast( _bc ) ) return false;
        return true;
      }
      virtual void synchronize()
      {
        mpi::AllGather allgather( mpi::main );
        typename t_Container :: iterator i_res = new_optima.begin();
        typename t_Container :: iterator i_res_end = new_optima.end();
        for(; i_res != i_res_end; ++i_res )
          i_res->broadcast( allgather );
        allgather.allocate_buffers();
        i_res = new_optima.begin();
        for(; i_res != i_res_end; ++i_res )
          i_res->broadcast( allgather );
        allgather();
        new_optima.clear();
        t_Individual indiv;
        while( indiv.broadcast( allgather ) )
        {
          if( condition( indiv ) ) continue;
          
          // checks wether result should be stored
          typename t_Container :: iterator i_found;
          i_found = std::find( results.begin(), results.end(), indiv );
          if ( i_found != results.end() )
            return;
          new_results = true;
          results.push_back( indiv );
        }
        results.remove_if( condition );
      }
#endif
  };

  namespace Condition
  {
    // Condition checks fitness only
    template< class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    class BaseOptima 
    {
      public:
        typedef T_EVALUATOR t_Evaluator;
        typedef T_GA_TRAITS t_GA_Traits;

      protected:
        typedef typename t_GA_Traits :: t_Individual t_Individual;
        typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
        typedef darwin::SaveObject<t_Evaluator> t_SaveOp;
        typedef darwin::LoadObject<t_Evaluator> t_LoadOp;

      protected:
        t_Individual optimum;

      public:
        BaseOptima   ( const TiXmlElement &_node ) {};
        BaseOptima   ( const BaseOptima &_c ) : optimum(_c.optimum) {};
        ~BaseOptima() {};

        bool Restart( const TiXmlElement &_node, t_LoadOp & _op)
        {
          const TiXmlElement *child = _node.FirstChildElement("optimum");
          if( not child ) return false;
          return optimum.Load( *child, _op );
        }
        bool Save( TiXmlElement &_node, t_SaveOp & _op) const
        {
          if ( optimum.invalid() )
          {
            std::cerr << "Optimum is invalid!! when trying to save" << std::endl;
            Print::out << "Optimum is invalid!! when trying to save\n";
            return false;
          }

          TiXmlElement *child = new TiXmlElement("optimum");
          if( not child )
          { 
            std::cerr << "Memory allocation error while save results" << std::endl;
            return false;
          }
          if ( not optimum.Save( *child, _op ) )
          {
            std::cerr << "Could not save optimum" << std::endl;
            return false;
          }

          _node.LinkEndChild(child);
          return true;
        }
#ifdef _MPI
        bool broadcast( mpi::BroadCast &_bc )
        {
          return optimum.broadcast(_bc);
        }
#endif
        std::string what_is() const { return " BaseOptima "; } 
        std::string print() const
        {
          std::ostringstream sstr;
          std::string bitstring; bitstring <<  (const typename t_IndivTraits :: t_Object&) optimum;
          sstr << std::setw(12) << std::setprecision(7) 
               << bitstring << " "
               << optimum.get_concentration() << " "
               << optimum.fitness() << std::endl;
          return sstr.str(); 
        }
    };
    // Condition checks objective 
    template< class T_EVALUATOR, class T_GA_TRAITS = Traits :: GA<T_EVALUATOR> >
    class FromObjective : public BaseOptima<T_EVALUATOR>
    {
      public:
        typedef T_EVALUATOR t_Evaluator;
        typedef T_GA_TRAITS t_GA_Traits;

      protected:
        typedef BaseOptima<t_Evaluator, t_GA_Traits> t_Base;
        typedef typename t_GA_Traits :: t_Individual t_Individual;
        typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
        typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
        typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
        typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;
        typedef typename Objective::Types<t_Evaluator, t_GA_Traits> t_Objective;
        typedef darwin::SaveObject<t_Evaluator> t_SaveOp;
        typedef darwin::LoadObject<t_Evaluator> t_LoadOp;

      protected:
        using t_Base :: optimum;
        typename t_Objective :: Vector *objective;
        t_ScalarQuantity val;
        t_ScalarQuantity end_val;
        t_ScalarQuantity delta;
        bool owns_objective;

      public:
        FromObjective   () : objective(NULL), owns_objective(false) {}
        FromObjective   ( const TiXmlElement &_node ) 
                      : t_Base(_node), objective(NULL),
                        val(0), end_val(0), delta(0), owns_objective(false)
        {
          std::string obj_type;
          double d=0;
          if ( _node.Attribute("delta") )
            _node.Attribute("delta", &d);
          delta = (t_ScalarQuantity) std::abs(d);

          objective = t_Objective :: new_from_xml( _node );
          if ( not objective )
            throw std::runtime_error("Could not Load objective from input in conditional store\n"); 
          owns_objective = true;
        }
        FromObjective   ( typename t_Objective::Vector* _type, const TiXmlElement &_node )
                      : t_Base(_node), objective(_type),
                        val(0), end_val(0), delta(0), owns_objective(false)
        {
          std::string obj_type;
          double d=0;
          if ( _node.Attribute("delta") )
            _node.Attribute("delta", &d);
          delta = (t_ScalarQuantity) std::abs(d);
 
          if ( not objective )
            throw std::runtime_error("Could not Load objective from input in conditional store\n"); 
        }
        FromObjective   ( const FromObjective &_c )
                      : t_Base(_c), objective(_c.objective),
                        val(_c.val), end_val(_c.end_val), delta(_c.delta),
                        owns_objective(false) {}
        ~FromObjective()  { if ( objective and owns_objective ) delete objective; }

        bool operator()( const t_Individual &_indiv )
        {
          objective->init( _indiv );
          if ( optimum.invalid() )
          {
            optimum = _indiv;
            val = (*objective)( optimum.quantities() );
            end_val = val + delta;
            return false;
          }
          t_ScalarQuantity indiv_val = (*objective)(_indiv.quantities());
          if ( t_QuantityTraits::less(indiv_val, val) )
          {
            optimum = _indiv;
            val = indiv_val;
            end_val = val + delta;
            return false;
          }

          return t_QuantityTraits::greater(indiv_val, end_val); 
        }
        std::string print() const
          { return objective->print();  }
        std::string what_is() const
        { 
          std::ostringstream sstr;
          sstr << objective->what_is() << " Objective, with delta = " << delta;
          return sstr.str();
        } 
        bool Restart( const TiXmlElement &_node, t_LoadOp & _op)
          { return t_Base::Restart( _node, _op) and objective->Restart( _node, _op); }
        bool Save( TiXmlElement &_node, t_SaveOp & _op) const
        {
          if (     t_Base::Save( _node, _op)  
               and objective->Save( _node, _op) )  return true;
          
          std::cerr << "Could not save objective" << std::endl;
          return false;
        }
    };

    // returns true if should remove
    // ie false if shouldn't store;
    template< class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    class Optima : public BaseOptima<T_EVALUATOR, T_GA_TRAITS>
    {
      public:
        typedef T_EVALUATOR t_Evaluator;
        typedef T_GA_TRAITS t_GA_Traits;
      protected:
        typedef BaseOptima<t_Evaluator, t_GA_Traits> t_Base;
        typedef typename t_GA_Traits :: t_Individual t_Individual;
        typedef typename t_GA_Traits :: t_IndivTraits t_IndivTraits;
        typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
        typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
        typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;

      protected:
        using t_Base :: optimum;

      public:
        Optima ( const TiXmlElement &_node ) : t_Base(_node) {}
        ~Optima () {}

        bool operator()( const t_Individual &_indiv )
        {
          if (  (not optimum.invalid()) and _indiv > optimum ) return true; 
          optimum = _indiv;
          return false; 
        }
        std::string what_is() const { return " optimum "; }
    };
  } // namespace Condition

  template<class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    struct Types
    {
      typedef Store::Conditional< T_EVALUATOR,
                                  Store::Condition::Optima< T_EVALUATOR, T_GA_TRAITS >,
                                  T_GA_TRAITS > Optima;
      typedef Store::Conditional< T_EVALUATOR,
                                  Store::Condition::FromObjective< T_EVALUATOR, T_GA_TRAITS >,
                                  T_GA_TRAITS >  FromObjective;
      typedef Store::Base< T_EVALUATOR, T_GA_TRAITS > Base;
    };
} // namespace Store

#endif // _RESULTS_H_

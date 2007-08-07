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

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

#include "objective.h"
#include "taboos.h"
#include "print_xmgrace.h"
#include "loadsave.h"
#include "gatraits.h"

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
      bool new_results;

    public:
      Base (t_Evaluator &_eval) : evaluator(_eval), new_results(false) {};
      Base (const Base<t_Evaluator> &_c) : evaluator(_c.eval), new_results(_c.new_results) {};
      virtual ~Base() {};

      virtual bool Restart( const TiXmlElement &_node ) = 0;
      virtual void Save( TiXmlElement &_node ) const = 0;

      virtual void operator()( t_Individual &_indiv ) = 0;
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
      virtual ~Conditional() {}

      // non-mpi version simply stores best N results 
      // MPI version is a bit of a hack
      // the object is for MPI version to be able to store in "unsynchronized" new_optima by default
      // and then into "synchronized" result upon call to synchronize
      virtual void operator()( t_Individual &_indiv )
      {
        // if first evaluation, adds it directly
        if( condition( _indiv ) )
          return;

        // checks wether result should be stored
        typename t_Container :: iterator i_found;
#ifdef _MPI
        i_found = std::find( new_optima.begin(), new_optima.end(), _indiv );
        if ( i_found != new_optima.end() )
          return;
#endif
        i_found = std::find( results.begin(), results.end(), _indiv );
        if ( i_found != results.end() )
          return;
        new_results = true;
#ifdef _MPI
        new_optima.push_front( _indiv );
#else
        results.push_front( _indiv );
#endif
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        darwin::LoadObject<t_Evaluator> loadworst( evaluator, &t_Evaluator::Load,
                                                   darwin::LOADSAVE_SHORT);
        if ( not condition.Restart( *xmlresults, loadworst ) )
        {
          darwin::printxmg.add_comment("Could not load condition");
          darwin::printxmg.add_comment("Aborting load of results");
          return false;
        }

        results.clear();
        darwin::LoadObject<t_Evaluator> loadop( evaluator, &t_Evaluator::Load,
                                                darwin::LOADSAVE_LONG);
        darwin::LoadIndividuals( *xmlresults, loadop, results );
        std::ostringstream sstr;
        sstr << "Reloaded Optimum ";
        sstr << " and " << results.size() << " target results";
        darwin::printxmg.add_comment( sstr.str() );
        print_results(0, true);
        return true;
      }
      void Save( TiXmlElement &_node ) const
      {
        darwin::SaveObject<t_Evaluator> saveworst( evaluator, &t_Evaluator::Save,
                                                   darwin::LOADSAVE_SHORT);
        darwin::SaveObject<t_Evaluator> saveop( evaluator, &t_Evaluator::Save,
                                                darwin::LOADSAVE_LONG);
        TiXmlElement *parent = new TiXmlElement("Results");
        if ( not parent ) return;
        condition.Save(*parent, saveworst);
        SaveIndividuals( *parent, saveop, results.begin(), results.end() );
        _node.LinkEndChild(parent);
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
      {
        typename t_Container :: const_iterator i_indiv = results.begin();
        typename t_Container :: const_iterator i_end = results.end();
        for (; i_indiv != i_end; ++i_indiv )
        {
          std::ostringstream sstr; 
          sstr << std::setw(12) << std::setprecision(7)
               << _age << " "
               << i_indiv->get_concentration() << " "
               << i_indiv->fitness();
          is_comment ?   darwin::printxmg.add_comment( sstr.str() )
                       : darwin::printxmg.add_line( sstr.str() );
        }
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
          if( condition( indiv ) )
            continue;
          
          // checks wether result should be stored
          typename t_Container :: iterator i_found;
          i_found = std::find( results.begin(), results.end(), indiv );
          if ( i_found != results.end() )
            return;
          new_results = true;
          results.push_front( indiv );
        }
      }
#endif
  };

  namespace Condition
  {
    // Condition checks objective 
    template< class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
    class BaseOptima 
    {
      public:
        typedef T_INDIVIDUAL t_Individual;
        typedef T_INDIV_TRAITS t_IndivTraits;

      protected:
        t_Individual optimum;

      public:
        BaseOptima   ( const TiXmlElement &_node ) {};
        ~BaseOptima() {};

        template< class T_OP > bool Restart( const TiXmlElement &_node, T_OP & _op)
        {
          const TiXmlElement *child = _node.FirstChildElement("optimum");
          if( not child ) return false;
          return optimum.Load( *child, _op );
        }
        template< class T_OP > bool Save( TiXmlElement &_node, T_OP & _op) const
        {
          TiXmlElement *child = new TiXmlElement("optimum");
          if( not child ) return false;
          if ( not optimum.Save( *child, _op ) ) return false;
          _node.LinkEndChild(child);
          return true;
        }
#ifdef _MPI
        bool broadcast( mpi::BroadCast &_bc )
        {
          return optimum.broadcast(_bc);
        }
#endif
    };
    // Condition checks objective 
    template< class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits :: Indiv<T_INDIVIDUAL> >
    class FromObjective : public BaseOptima<T_INDIVIDUAL>
    {
      public:
        typedef T_INDIVIDUAL t_Individual;
        typedef T_INDIV_TRAITS t_IndivTraits;

      protected:
        typedef BaseOptima<t_Individual, t_IndivTraits> t_Base;
        typedef typename t_IndivTraits :: t_QuantityTraits t_QuantityTraits;
        typedef typename t_QuantityTraits :: t_Quantity t_Quantity;
        typedef typename t_QuantityTraits :: t_ScalarQuantity t_ScalarQuantity;

      protected:
        using t_Base :: optimum;
        ::Objective::Base<t_ScalarQuantity> *objective;
        types::t_unsigned n;
        types::t_real val;
        types::t_real end_val;
        types::t_real delta;

      public:
        FromObjective   ( const TiXmlElement &_node ) 
                      : t_Base(_node), objective(NULL), n(0), 
                        val(0), end_val(0), delta(0)
        {
          std::string obj_type;
          int i = 0; double d=0;
          if ( _node.Attribute("n") )
            _node.Attribute("n", &i);
          n = (types::t_unsigned ) std::abs(i);
          if ( _node.Attribute("delta") )
            _node.Attribute("delta", &d);
          delta = (types::t_real) std::abs(d);

          objective = Objective::new_from_xml( _node );
          if ( not objective )
            throw std::runtime_error("Could not Load objective from input in conditional store\n"); 

        }
        ~FromObjective()  { if ( objective ) delete objective; }

        bool operator()( const t_Individual &_indiv )
        {
          if ( optimum.invalid() )
          {
            optimum = _indiv;
            val = (*objective)( optimum.quantities(n) );
            end_val = val + delta;
            return false;
          }
          types::t_real indiv_val = (*objective)(_indiv.quantities(n));
          if ( indiv_val < val )
          {
            optimum = _indiv;
            val = indiv_val;
            end_val = val + delta;
            return false;
          }

          return indiv_val >= end_val; 
        }
    };

    // returns true if should remove
    // ie false if shouldn't store;
    template< class T_INDIVIDUAL, class T_INDIV_TRAITS = Traits::Indiv<T_INDIVIDUAL> >
    class Optima : public BaseOptima<T_INDIVIDUAL, T_INDIV_TRAITS>
    {
      public:
        typedef T_INDIVIDUAL t_Individual;
        typedef T_INDIV_TRAITS t_IndivTraits;
      protected:
        typedef BaseOptima<t_Individual, t_IndivTraits> t_Base;
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
    };
  } // namespace Condition

  template<class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    struct Optima
    {
      typedef Store::Conditional< T_EVALUATOR,
                                  Store::Condition::Optima< typename T_GA_TRAITS :: t_Individual,
                                                            typename T_GA_TRAITS :: t_IndivTraits>,
                                  T_GA_TRAITS > auto_template;
    };
  template<class T_EVALUATOR, class T_GA_TRAITS = Traits::GA<T_EVALUATOR> >
    struct FromObjective
    {
      typedef Store::Conditional< T_EVALUATOR,
                                  Store::Condition::FromObjective< typename T_GA_TRAITS :: t_Individual,
                                                                   typename T_GA_TRAITS :: t_IndivTraits>,
                                  T_GA_TRAITS >  auto_template;
    };
} // namespace Store

#endif // _RESULTS_H_

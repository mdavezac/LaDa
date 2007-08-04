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

#include "taboos.h"
#include "print_xmgrace.h"
#include "loadsave.h"

namespace Store
{ 
  // abstract base class for results and storage
  template<class T_EVALUATOR>
  class Base
  {
    public:
      typedef T_EVALUATOR t_Evaluator;

    private:
      typedef typename t_Evaluator::t_Individual t_Evaluator;
      typedef typename t_Individual::t_Object t_Object;

    protected:
      t_Evaluator &evaluator;

    public:
      Store (const t_Evaluator &_eval) : evaluator(_eval) {};
      virtual ~Store() {};

      virtual bool Restart( const TiXmlElement &_node ) = 0;
      virtual void Save( TiXmlElement &_node ) const = 0;

      virtual void operator()( t_Individual &_indiv ) = 0;
#ifdef _MPI
      virtual void synchronize();
#endif
  };

  // optima = best result
  template<class T_EVALUATOR>
  class Optima : public Base<T_EVALUATOR>
  {
    public:
      typedef T_Evaluator t_Evaluator;
    private:
      typedef typename t_Evaluator::t_Individual t_Individual;
      typedef typename t_Individual::t_Object t_Object;

    protected:
      // returns false if shouldn't store
    class Conditional
    {
      protected:
        const t_Individual optimum;
      public:
        Conditional () {}

        bool operator()( const t_Individual &_indiv ) const
        {
          if ( not optimum.invalid() )
            return _indiv > optimum; 
        }
    };

    protected:
      Conditional conditional;
      std::vector<t_Individual> results;
#ifdef _MPI
      std::vector<t_Individual> new_optima;
#endif 

    public:
      Optima   (t_Evaluator &_eval)
             : Base<t_Evaluator>( _eval)
      {
        printxmg.add_comment("Store Optimum");
      };
      virtual ~BestN() {}

      // non-mpi version simply stores best N results 
      // MPI version is a bit of a hack
      // the object is for MPI version to be able to store in "unsynchronized" new_optima by default
      // and then into "synchronized" result upon call to synchronize
#ifndef _MPI
      virtual void operator()( t_Individual &_indiv )
#else
      virtual void operator()( t_Individual &_indiv,
                               std::vector<t_Individual> &_container = new_optima )
#endif
      {
#ifndef _MPI
        std::vector<t_Individual> &_container = results;
#endif
        // if first evaluation, adds it directly
        if( conditional( _indiv ) )
          return;

        // checks wether result should be stored
        typename std::list<t_Individual> :: iterator i_found;
        i_found = std::find( _container.begin(), _container.end(), _indiv );
        if ( i_found != results.end() )
          return;
        new_results = true;
        _container.push_front( _indiv );
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        const TiXmlElement *xmloptimum = xmlresults->FirstChildElement("Optimum");
        if ( not xmloptimum ) return false;
        LoadObject<t_Individual, t_Evaluator> loadworst( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
        if ( not conditional.optimum.Load( *xmloptimum, loadworst ) )
        {
          conditional.optimum.invalidate();
          printxmg.add_comment("Could not load worst individual");
          printxmg.add_comment("Aborting load of results");
          return false;
        }

        results.clear();
        LoadObject<t_Individual, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_LONG);
        LoadIndividuals( *xmlresults, loadop, results );
        std::ostringstream sstr;
        sstr << "Reloaded Optimum ";
        sstr << " and " << results.size() << " target results";
        printxmg.add_comment( sstr.str() );
        print_results(0, true);
        return true;
      }
      void Save( TiXmlElement &_node ) const
      {
        SaveObject<t_Individual, t_Evaluator> saveworst( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
        SaveObject<t_Individual, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_LONG);
        TiXmlElement *parent = new TiXmlElement("Results");
        TiXmlElement *child = new TiXmlElement("Optimum");
        conditional.optimum.Save(*child, saveworst);
        parent->LinkEndChild(child);
        SaveIndividuals( *parent, saveop, results.begin(), results.end() );
        _node.LinkEndChild(parent);
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
      {
        typename std::list< t_Individual > :: const_iterator i_indiv = results.begin();
        typename std::list< t_Individual > :: const_iterator i_end = results.end();
        for (; i_indiv != i_end; ++i_indiv )
        {
          std::ostringstream sstr; 
          sstr << std::setw(12) << std::setprecision(7)
               << _age << " "
               << i_indiv->get_concentration() << " "
               << i_indiv->get_quantity();
          is_comment ? printxmg.add_comment( sstr.str() ) : printxmg.add_line( sstr.str() );
        }
      }
#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not optimum.broadcast(_bc) ) return false;
        types::t_int n = results.size();
        if ( not _bc.serialize(n) ) return false;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          results.resize(n);
        typename std::list<t_Individual> :: iterator i_res = results.begin();
        typename std::list<t_Individual> :: iterator i_res_end = results.end();
        for(; i_res != i_res_end; ++i_res )
          if( not i_res->broadcast( _bc ) ) return false;
        return true;
      }
      virtual void synchronize( t_Population &_pop )
      {
        mpi::AllGather allgather( mpi::main );
        typename std::list<t_Individual> :: iterator i_res = new_optima.begin();
        typename std::list<t_Individual> :: iterator i_res_end = new_optima.end();
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
          operator()( indiv, results );
      }
#endif
  };

  //  All N best results
  template<class T_EVALUATOR>
  class FromObjective : public Base<T_EVALUATOR>
  {
    public:
      typedef T_EVALUATOR  t_Evaluator;
    private:
      typedef typename t_Evaluator::t_Individual t_Individual;
      typedef typename t_Individual::t_Object t_Object;

    protected:
      using Optimum<t_Evaluator> :: optimum;

    // Conditional returns false if indiv should not be stored
    class Conditional
    {
      protected:
        const objective::Base<typename t_Individual::t_Quantity> *obj;
        types::t_unsigned n;
        t_Individual optimum;
        types::t_real val;
        types::t_real end_val;
        types::t_real delta

      public:
        Conditional  ( const objective *obj, types::t_unsigned _i, types::t_real _delta )
                    : obj( _obj ), n(_i), delta(_delta) {}

        bool operator()( const t_Individual &_indiv )
        {
          if ( optimum.invalid() )
          {
            optimum = _indiv;
            val = (*obj)( optimum.quantities(n) );
            end_val = val + delta;
            return false;
          }
          types::t_real indiv_val = (*obj)(_indiv.quantities(n));
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

    protected:
      Conditional conditional;
      std::vector<t_Individual> results;
#ifdef _MPI
      std::vector<t_Individual> new_optima;
#endif 

    public:
      FromObjective (t_Evaluator &_eval, t_Objective *objective,
                     types::t_unsigned _n = 0, types::t_real _delta = 0)
             : Base<t_Evaluator>( _eval), conditional(objective, _n, _delta )
      {
        std::ostringstream sstr;
        sstr << "Store results at most " << _delta << " from objective " << _n;
        printxmg.add_comment(sstr.str());
      };
      virtual ~BestN() {}

      // non-mpi version simply stores best N results 
      // MPI version is a bit of a hack
      // the object is for MPI version to be able to store in "unsynchronized" new_optima by default
      // and then into "synchronized" result upon call to synchronize
#ifndef _MPI
      virtual void operator()( t_Individual &_indiv )
#else
      virtual void operator()( t_Individual &_indiv,
                               std::vector<t_Individual> &_container = new_optima )
#endif
      {
#ifndef _MPI
        std::vector<t_Individual> &_container = results;
#endif
        // if first evaluation, adds it directly
        if ( conditional( _indiv ) )
          return; 
        typename std::list<t_Individual> :: iterator i_found;
        i_found = std::find( results.begin(), results.end(), _indiv );
        if ( i_found != results.end() )
          return;
        new_results = true;
        _container.push_front( _indiv );
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        const TiXmlElement *xmloptimum = xmlresults->FirstChildElement("worst");
        if ( not xmloptimum ) return false;
        LoadObject<t_Individual, t_Evaluator> loadworst( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
        if ( not conditioanl.optimum.Load( *xmloptimum, loadworst ) )
        {
          conditional.optimum.invalidate();
          printxmg.add_comment("Could not load worst individual");
          printxmg.add_comment("Aborting load of results");
          return false;
        }

        results.clear();
        LoadObject<t_Individual, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_LONG);
        LoadIndividuals( *xmlresults, loadop, results );
        std::ostringstream sstr;
        sstr << "Reloaded Optimum ";
        sstr << " and " << results.size() << " target results";
        printxmg.add_comment( sstr.str() );
        print_results(0, true);
        return true;
      }
      void Save( TiXmlElement &_node ) const
      {
        SaveObject<t_Individual, t_Evaluator> saveworst( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
        SaveObject<t_Individual, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_LONG);
        TiXmlElement *parent = new TiXmlElement("Results");
        TiXmlElement *child = new TiXmlElement("optimum");
        conditional.optimum.Save(*child, saveworst);
        parent->LinkEndChild(child);
        parent->LinkEndChild(child);
        SaveIndividuals( *parent, saveop, results.begin(), results.end() );
        _node.LinkEndChild(parent);
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
      {
        typename std::list< t_Individual > :: const_iterator i_indiv = results.begin();
        typename std::list< t_Individual > :: const_iterator i_end = results.end();
        for (; i_indiv != i_end; ++i_indiv )
        {
          std::ostringstream sstr; 
          sstr << std::setw(12) << std::setprecision(7)
               << _age << " "
               << i_indiv->get_concentration() << " "
               << i_indiv->get_quantity();
          is_comment ? printxmg.add_comment( sstr.str() ) : printxmg.add_line( sstr.str() );
        }
      }
#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not optimum.broadcast(_bc) ) return false;
        types::t_int n = results.size();
        if ( not _bc.serialize(n) ) return false;
        if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          results.resize(n);
        typename std::list<t_Individual> :: iterator i_res = results.begin();
        typename std::list<t_Individual> :: iterator i_res_end = results.end();
        for(; i_res != i_res_end; ++i_res )
          if( not i_res->broadcast( _bc ) ) return false;
        return true;
      }
      virtual void synchronize( t_Population &_pop )
      {
        mpi::AllGather allgather( mpi::main );
        typename std::list<t_Individual> :: iterator i_res = new_optima.begin();
        typename std::list<t_Individual> :: iterator i_res_end = new_optima.end();
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
          operator()( indiv, results );
      }
#endif
  };
} // namespace Store

#endif // _RESULTS_H_

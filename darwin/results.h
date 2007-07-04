#ifndef _DARWIN_RESULTS_H_
#define _DARWIN_RESULTS_H_

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

namespace darwin 
{ 
  // abstract base class for results and storage
  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class Results
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;
    private:
      typedef typename t_Individual::t_Object t_Object;
      typedef typename t_Evaluator::t_Functional t_Functional;
      typedef typename t_Functional::t_Type t_Type;
      typedef Taboo_Base<t_Individual> t_Taboo;

    protected:
      History< t_Individual, std::list<t_Individual> > *history;
      t_Evaluator &evaluator;
      t_Functional *objective;
      mutable bool new_results;
      t_Individual *current_individual;
      const t_Taboo *const taboo;

    public:
      types::t_unsigned nb_eval;
      types::t_unsigned nb_grad;

    public:
      Results   ( t_Evaluator &_eval, const t_Taboo * const _t)
              : history(NULL), evaluator(_eval), objective(NULL), 
                new_results(false), current_individual(NULL),
                taboo(_t), nb_eval(0), nb_grad(0) { };
      virtual ~Results() {};

      void set_evaluator ( t_Evaluator *_e ) { evaluator = _e; }
      void set_history( History< t_Individual, std::list<t_Individual> > *_h ) { history = _h; }
      bool is_new_results() const 
        { return new_results; }
      void set_new_results(bool _n = false) const 
        { new_results = _n; }
      void set_object( t_Individual _indiv )
        { evaluator.set_object( _indiv, objective ); }
      typename t_Functional::t_Container* const get_objective_variables()
        { return objective->get_variables(); }
      void invalidate()
        { current_individual->invalidate(); }
      
      virtual void evaluate( t_Individual &_indiv ) = 0;
      virtual void evaluate_with_gradient( t_Type *_i_grad, t_Individual &_indiv ) = 0;
      virtual void print_results(types::t_unsigned _age, bool is_comment = false ) const = 0;

      virtual bool Restart( const TiXmlElement &_node ) = 0;
      virtual void Save( TiXmlElement &_node ) const = 0;
      virtual void printout( std::ostream &_stream ) {};
      virtual void printout( const std::string &_filename) {};

      typename t_Functional::t_Type evaluate()
      {
        evaluator.set_object( *current_individual, objective );
        evaluate(*current_individual);
        return current_individual->fitness();
      }
      typename t_Functional::t_Type evaluate_with_gradient( typename t_Functional::t_Type *_i_grad )
      {
        evaluator.set_object( *current_individual, objective );
        evaluate_with_gradient(_i_grad, *current_individual);
        return current_individual->fitness();
      }
      void evaluate_gradient( typename t_Functional::t_Type *_i_grad )
      {
        nb_grad += objective->size();
        evaluator.set_object( *current_individual, objective );
        objective->evaluate_gradient(_i_grad); 
      }
      typename t_Functional::t_Type evaluate_one_gradient( types::t_unsigned _pos )
      {
        ++nb_grad; 
        evaluator.set_object( *current_individual, objective );
        return objective->evaluate_one_gradient(_pos); 
      }
      bool is_taboo() const
      {
        if ( not taboo )
          return false;
        evaluator.set_object( *current_individual, objective );
        return (*taboo)(*current_individual);
      }
      bool init() { return true; }; // done below
      // this function should fully set-up the objective function
      virtual void init( t_Individual &_indiv ) = 0;
        
      virtual void evaluate( t_Population &_parents, t_Population &_offsprings ) 
      {
        typename t_Population :: iterator i_indiv = _offsprings.begin();
        typename t_Population :: iterator i_last = _offsprings.end();

        for ( ; i_indiv != i_last; ++i_indiv )
        {
          init( *i_indiv );
          evaluate(*i_indiv);
        }

#ifdef _MPI
        synchronize( _offsprings ); // for mpi purposes
#endif

        i_indiv = _offsprings.begin();
        std::cout << "New Individuals:" << std::endl; 
        for (types::t_int i = 0 ; i_indiv != i_last; ++i, ++i_indiv )
          std::cout << " Offspring " << i 
                    << " Fitness: " << i_indiv->fitness() << std::endl;
        std::cout << std::endl; 
      }
#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc ) = 0;
      virtual void synchronize( t_Population &_pop ) {};
#endif

    protected:
      virtual bool is_known( t_Individual &_indiv )
      {
        if ( not history ) 
          return false;
        return history->set_quantity(_indiv);
      }

  };

  //  only one result = optimum
  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class Optimum : public Results<T_INDIVIDUAL, T_EVALUATOR, T_POPULATION>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;
    private:
      typedef typename t_Evaluator :: t_Functional t_Functional;
      typedef typename t_Functional::t_Type t_Type;
      typedef Taboo_Base<t_Individual> t_Taboo;
      typedef typename t_Individual :: t_Object t_Object;
    protected:
      using Results<t_Individual, t_Evaluator, t_Population> :: evaluator;
      using Results<t_Individual, t_Evaluator, t_Population> :: objective;
      using Results<t_Individual, t_Evaluator, t_Population> :: history;
      using Results<t_Individual, t_Evaluator, t_Population> :: new_results;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_eval;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_grad;
      using Results<t_Individual, t_Evaluator, t_Population> :: current_individual;

    protected:
      t_Individual optimum;

    public:
      Optimum   (const TiXmlElement& _parent, t_Evaluator &_eval, const t_Taboo* const  _t) 
              : Results<t_Individual, t_Evaluator, t_Population>( _eval, _t )
        { printxmg.add_comment("Objective: Find Optimum"); }
      virtual ~Optimum() {};

      // sets up objective function fully
      void init( t_Individual &_indiv )
      {
        current_individual = &_indiv;
        objective = (t_Functional*) evaluator.init(*current_individual); 
      }

      virtual void evaluate( t_Individual &_indiv )
      {
        if ( not _indiv.invalid() )
          return;
        
        // returns true if fitness can be set
        if ( is_known(_indiv) ) 
        {
          types::t_real quantity = _indiv.get_quantity();
          _indiv.set_fitness( quantity );
          return;
        }

        ++nb_eval;
        types::t_real quantity = objective->evaluate();
        _indiv.set_quantity( quantity );
        _indiv.set_fitness( quantity );
        if ( history )
          history->add(_indiv);
        if ( optimum.invalid() or _indiv > optimum ) 
        {
          new_results = true;
          optimum = _indiv;
        }
        return;
      }
      virtual void evaluate_with_gradient( t_Type *_i_grad, t_Individual &_indiv )
      {
        if ( not _indiv.invalid() )
          return;
        
        // returns true if fitness can be set
        if ( is_known(_indiv) ) 
        {
          types::t_real quantity = _indiv.get_quantity();
          _indiv.set_fitness( quantity );
          return;
        }

        ++nb_eval;
        nb_grad += objective->size();
        types::t_real quantity = objective->evaluate_with_gradient( _i_grad );
        _indiv.set_quantity( quantity );
        _indiv.set_fitness( quantity );
        if ( history )
          history->add(_indiv);
        if ( optimum.invalid() or _indiv > optimum ) 
        {
          new_results = true;
          optimum = _indiv;
        }
        return;
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false ) const
      {
        std::ostringstream sstr; 
        sstr << std::setw(12) << std::setprecision(7)
             << _age << " "
             << optimum.get_concentration() << " "
             << optimum.get_quantity();
        is_comment ? printxmg.add_comment( sstr.str() ) : printxmg.add_line( sstr.str() );
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        const TiXmlElement *xmloptimum = xmlresults->FirstChildElement("Optimum");
        if ( not xmloptimum )
          return false;
        LoadObject<t_Object, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_LONG);
        if( not optimum.Load( *xmloptimum, loadop ) );  
        printxmg.add_comment("Loading Optimum from Input");
        print_results(0, true);
        return true;
      }
      void Save( TiXmlElement &_node ) const
      {
        TiXmlElement *parent = new TiXmlElement("Results");
        TiXmlElement *child = new TiXmlElement("Optimum");
        SaveObject<t_Object, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_LONG);
        optimum.Save(*child, saveop);
        parent->LinkEndChild(child);
        _node.LinkEndChild(parent);
      }

#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc )
      {
        if ( not optimum.broadcast(_bc) ) return false;
        return true;
      }
      virtual void synchronize( t_Population &_pop )
      {
        mpi::BroadCast bc( mpi::main );
        for( types::t_int i = 0; i < mpi::main.size(); ++i )
        {
          if ( i == mpi::main.rank() and new_results )
            optimum.broadcast(bc);
          if ( not bc.allocate_buffers() ) continue; // nothing to broadcast!!
          if ( i == mpi::main.rank() and new_results )
            optimum.broadcast(bc);
          bc();
          t_Individual indiv;
          indiv.broadcast(bc);
          bc.reset();
          if ( optimum.invalid() or  indiv > optimum )
          {
            std::cout << "Synchronize: Optimum Invalid " << std::endl;
            optimum = indiv;
          }
        }
      }
#endif
  };

  //  All results less than delta from optimum
  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class BestOf : public Results<T_INDIVIDUAL, T_EVALUATOR, T_POPULATION>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;
    private:
      typedef typename t_Evaluator::t_Functional t_Functional;
      typedef typename t_Functional::t_Type t_Type;
      typedef Taboo_Base<t_Individual> t_Taboo;
      typedef typename t_Individual :: t_Object t_Object;
    protected:
      using Results<t_Individual, t_Evaluator, t_Population> :: evaluator;
      using Results<t_Individual, t_Evaluator, t_Population> :: objective;
      using Results<t_Individual, t_Evaluator, t_Population> :: history;
      using Results<t_Individual, t_Evaluator, t_Population> :: new_results;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_eval;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_grad;
      using Results<t_Individual, t_Evaluator, t_Population> :: current_individual;

    class Remove_If
    {
      protected:
        const types::t_real &value;
        const types::t_real &delta;
      public:
        Remove_If   ( const types::t_real &_val, const types::t_real &_delta )
                  : value(_val), delta(_delta) {};

        bool operator()( const t_Individual &_obj ) const
         { return std::abs( _obj.get_quantity() - value ) > delta; }
    };

    protected:
      types::t_real delta;
      Remove_If remove_if;
      std::list<t_Individual> results;
#ifdef _MPI
      std::list<t_Individual> new_optima;
#endif 
      T_INDIVIDUAL optimum;

    public:
      BestOf   (const TiXmlElement &_parent, t_Evaluator &_eval, const t_Taboo * const _t)
             : Results<t_Individual, t_Evaluator, t_Population>( _eval, _t ), 
               delta(types::tolerance), remove_if(optimum.get_quantity(), delta)
      {
        if ( _parent.Attribute("delta") )
          _parent.Attribute("delta", &delta);
        std::ostringstream sstr;
        sstr << "Objective: Find Optimum and All Structures ";
        sstr << delta << " from Optimum";
        printxmg.add_comment(sstr.str());
      };
      virtual ~BestOf() {}

      void init( t_Individual &_indiv )
      {
        current_individual = &_indiv;
        objective = (t_Functional*) evaluator.init(*current_individual); 
      }

      // evaluator.init( _indiv ) MUST be called prior to 
      // this function!!
      virtual void evaluate( t_Individual &_indiv )
      {
        if ( not _indiv.invalid() )
          return;
        
        // returns true if fitness can be set
        if ( is_known(_indiv) ) 
        {
          types::t_real quantity = _indiv.get_quantity();
          _indiv.set_fitness( quantity );
          return;
        }

        ++nb_eval; 
        types::t_real quantity = objective->evaluate();
        _indiv.set_quantity( quantity );
        _indiv.set_fitness( quantity );
        do_store( _indiv );
      }

      virtual void evaluate_with_gradient( t_Type *_i_grad, t_Individual &_indiv )
      {
        if ( not _indiv.invalid() )
          return;
        
        // returns true if fitness can be set
        if ( is_known(_indiv) ) 
        {
          types::t_real quantity = _indiv.get_quantity();
          _indiv.set_fitness( quantity );
          return;
        }

        ++nb_eval; 
        nb_grad += objective->size(); 
        types::t_real quantity = objective->evaluate_with_gradient( _i_grad );
        _indiv.set_quantity( quantity );
        _indiv.set_fitness( quantity );
        do_store( _indiv );
      }

      void do_store( t_Individual &_indiv )
      {
        if ( history )
          history->add(_indiv);

        // if first evaluation, adds it directly
        if ( optimum.invalid() )
        {
          new_results = true;
          optimum = _indiv;
#ifdef _MPI
          new_optima.push_front( _indiv );
#else
          results.push_front( _indiv );
#endif

          return;
        }

        // checks wether result should be stored
        if ( _indiv > optimum )
        {
          optimum = _indiv;
          results.remove_if( remove_if );
          new_results = true;
#ifdef _MPI
          new_optima.push_front( _indiv );
#else
          results.push_front( _indiv );
#endif
          return;
        }

        if ( not remove_if(_indiv) )
        {
          typename std::list<t_Individual> :: iterator i_found;
          i_found = std::find( results.begin(), results.end(), _indiv );
          if ( i_found != results.end() )
            return;
          new_results = true;
#ifdef _MPI
          new_optima.push_front( _indiv );
#else
          results.push_front( _indiv );
#endif
        }  // if _indiv > optimum ...
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        const TiXmlElement *xmloptimum = xmlresults->FirstChildElement("Optimum");
        if ( not xmloptimum ) return false;
        LoadObject<t_Object, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_LONG);
        if ( not optimum.Load( *xmloptimum, loadop ) )
        {
          optimum.invalidate();
          printxmg.add_comment("Could not load optimum");
          printxmg.add_comment("Aborting load of results");
          return false;
        }

        results.clear();
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
        SaveObject<t_Object, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_LONG);
        TiXmlElement *parent = new TiXmlElement("Results");
        TiXmlElement *child = new TiXmlElement("Optimum");
        optimum.Save(*child, saveop);
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
          if ( optimum.invalid() or indiv > optimum )
          {
            optimum = indiv;
            results.remove_if( remove_if );
          }
          else if ( not remove_if(indiv) )
          {
            typename std::list<t_Individual> :: iterator i_found;
            i_found = std::find( results.begin(), results.end(), indiv );
            if ( i_found == results.end() )
              results.push_front( indiv );
          }  // if _indiv > optimum ...
      }
#endif
  };

  //  All results less than delta from target 
  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class Target : public BestOf<T_INDIVIDUAL, T_EVALUATOR, T_POPULATION>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;
    private:
      typedef typename t_Evaluator::t_Functional t_Functional;
      typedef Taboo_Base<t_Individual> t_Taboo;
    protected:
      using Results<t_Individual, t_Evaluator, t_Population> :: objective;
      using Results<t_Individual, t_Evaluator, t_Population> :: evaluator;
      using BestOf<t_Individual, t_Evaluator, t_Population> :: delta;
      using BestOf<t_Individual, t_Evaluator, t_Population> :: remove_if;
      using BestOf<t_Individual, t_Evaluator, t_Population> :: evaluate;
      using Results<t_Individual, t_Evaluator, t_Population> :: current_individual;

    protected:
      types::t_real target;

    public:
      Target   (const TiXmlElement &_parent, t_Evaluator &_eval, const t_Taboo* const _t)
             : BestOf<t_Individual, t_Evaluator, t_Population>( _parent, _eval, _t ), 
               target(0)
      {
        if ( _parent.Attribute("target") )
          _parent.Attribute("target", &target);
        std::ostringstream sstr;
        sstr << "Objective: Find all Structures within "
             << delta << " of target value " << target;
        printxmg.remove_last();
        printxmg.add_comment(sstr.str());

        delta *= delta;
        objective = new function::Quadratic<typename t_Evaluator::t_Functional>( target );
      };
      virtual ~Target() { delete objective; };

      void init( t_Individual &_indiv )
      {
        current_individual = &_indiv;
        static_cast< function::Quadratic<t_Functional>* >(objective)->
            set_functional( (t_Functional*) evaluator.init(*current_individual) );
      }
  };
  
  //  Convex hull results
  template<class T_INDIVIDUAL, class T_EVALUATOR, class T_POPULATION = eoPop<T_INDIVIDUAL> >
  class ConvexHull : public Results<T_INDIVIDUAL, T_EVALUATOR, T_POPULATION>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_EVALUATOR  t_Evaluator;
      typedef T_POPULATION t_Population;
    private:
      typedef typename t_Individual::t_Object t_Object;
      typedef typename t_Evaluator::t_Functional t_Functional;
      typedef typename t_Functional::t_Type t_Type;
      typedef ch::Base<t_Object> t_ConvexHull;
      typedef function::Fitness<t_Functional, t_ConvexHull > t_Fitness;
      typedef Taboo_Base<t_Individual> t_Taboo;
    protected:
      using Results<t_Individual, t_Evaluator, t_Population> :: evaluator;
      using Results<t_Individual, t_Evaluator, t_Population> :: objective;
      using Results<t_Individual, t_Evaluator, t_Population> :: history;
      using Results<t_Individual, t_Evaluator, t_Population> :: new_results;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_eval;
      using Results<t_Individual, t_Evaluator, t_Population> :: nb_grad;
      using Results<t_Individual, t_Evaluator, t_Population> :: current_individual;

    protected:
      t_ConvexHull convexhull;
      bool valid_ch;

    public:
      ConvexHull   (const TiXmlElement &_parent, t_Evaluator &_eval, const t_Taboo* const _t)
                 : Results<t_Individual, t_Evaluator, t_Population>( _eval, _t ), 
                   convexhull(), valid_ch(true)
      {
        printxmg.add_comment("Objective: Find Convex-Hull");
        objective = new t_Fitness;
        static_cast< t_Fitness* >(objective)->set_baseline( &convexhull );
        t_Individual pure;
        evaluator.initialize(pure);
        init(pure);
        typename t_Functional :: t_Container :: iterator i_var = objective->begin();
        typename t_Functional :: t_Container :: iterator i_end = objective->end();
        for( ; i_var != i_end; ++i_var )
          *i_var = -1.0;
        types::t_real Epure = static_cast< t_Fitness* >(objective)->get_quantity()->evaluate();
        evaluator.set_object( pure, objective );
        convexhull.add( Epure, pure );
        std::ostringstream sstr;
        sstr << "Adding points (" << pure.get_concentration() << "," << Epure << ") and (";
        for(i_var = objective->begin() ; i_var != i_end; ++i_var )
          *i_var = 1.0;
        Epure = static_cast< t_Fitness* >(objective)->get_quantity()->evaluate();
        evaluator.set_object( pure, objective );
        convexhull.add( Epure, pure );
        sstr << pure.get_concentration() << "," << Epure << ") to convexhull";
        printxmg.add_comment( sstr.str() );
      };
      virtual ~ConvexHull() { delete objective; };

      void init( t_Individual &_indiv )
      {
        current_individual = &_indiv;
        static_cast< t_Fitness* >(objective)->
            set_quantity( static_cast<t_Functional*>(evaluator.init(*current_individual)) );
      }
      virtual void evaluate( t_Individual &_indiv )
      {
        if ( _indiv.invalid() )
        {
          // returns true if fitness can be set
          if ( is_known(_indiv) ) 
          {
            types::t_real quantity = _indiv.get_quantity();
            types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
            _indiv.set_fitness( quantity - base );
            return;
          }

          ++nb_eval;
          types::t_real quantity = objective->evaluate();
          _indiv.set_fitness( quantity );
          types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
          quantity += base;
          _indiv.set_quantity( quantity );
          if ( history )
            history->add(_indiv);
          if ( quantity < base ) 
            if( convexhull.add( quantity, _indiv ) )
            {
              new_results = true;
              valid_ch = false;
            }
          return;
        }

        if ( valid_ch )
          return;
        
        types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
        types::t_real quantity = _indiv.get_quantity();
        _indiv.set_fitness( quantity - base );
        if ( quantity < base ) 
          if( convexhull.add( quantity, _indiv ) )
          {
            new_results = true;
            valid_ch = false;
          }
      }
      virtual void evaluate_with_gradient( t_Type *_i_grad, t_Individual &_indiv )
      {
        if ( _indiv.invalid() )
        {
          // returns true if fitness can be set
          if ( is_known(_indiv) ) 
          {
            types::t_real quantity = _indiv.get_quantity();
            types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
            _indiv.set_fitness( quantity - base );
            objective->evaluate_gradient( _i_grad );
            return;
          }

          ++nb_eval;
          types::t_real quantity = objective->evaluate_with_gradient( _i_grad );
          _indiv.set_fitness( quantity );
          types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
          quantity += base;
          _indiv.set_quantity( quantity );
          if ( history )
            history->add(_indiv);
          if ( quantity < base ) 
            if( convexhull.add( quantity, _indiv ) )
            {
              new_results = true;
              valid_ch = false;
            }
          return;
        }
        
        types::t_real base = convexhull.evaluate( _indiv.get_concentration() );
        types::t_real quantity = _indiv.get_quantity();
        _indiv.set_fitness( quantity - base );
        objective->evaluate_gradient( _i_grad );
        if ( quantity < base ) 
          if( convexhull.add( quantity, _indiv ) )
          {
            new_results = true;
            valid_ch = false;
          }
      }
      virtual void evaluate( t_Population &_parents, t_Population &_offsprings ) 
      {
        typename t_Population :: iterator i_indiv = _offsprings.begin();
        typename t_Population :: iterator i_last  = _offsprings.end();

        for ( ; i_indiv != i_last; ++i_indiv )
        {
          init(*i_indiv);
          evaluate(*i_indiv);
        }


        // convex hull has changed => reevaluate
#ifdef _MPI
        if ( not mpi::main.all_sum_all( valid_ch ) )
        {
          synchronize( _offsprings );
          i_last = _offsprings.end();
          if( mpi::main.rank() == mpi::ROOT_NODE )
#else
        if ( not valid_ch )
        {
#endif
          std::cout << std::endl << "Base line changed" << std::endl << std::flush; 
          i_indiv = _offsprings.begin();
          for ( ; i_indiv != i_last; ++i_indiv )
          {
            init(*i_indiv);
            evaluate(*i_indiv);
          }
          i_indiv = _parents.begin();
          i_last = _parents.end();
          for ( ; i_indiv != i_last; ++i_indiv )
          {
            init(*i_indiv);
            evaluate(*i_indiv);
          }
          valid_ch = true;
        } 
        i_indiv = _offsprings.begin();
        i_last = _offsprings.end();
        std::cout << "New Individuals:" << std::endl; 
        for (types::t_int i = 0 ; i_indiv != i_last; ++i, ++i_indiv )
          std::cout << " Offspring " << i 
                    << " Fitness: " << i_indiv->fitness() 
                    << " Quantity: " << i_indiv->get_quantity() 
                    << " Baseline: " << i_indiv->get_quantity() - i_indiv->fitness() << std::endl; 
        std::cout << std::endl << std::flush; 
      }

      virtual void print_results(types::t_unsigned _age, bool is_comment = false) const
      {
        typename t_ConvexHull :: t_Vertices :: const_iterator i_vertex = convexhull.begin_vertex();
        typename t_ConvexHull :: t_Vertices :: const_iterator i_end = convexhull.end_vertex();
        for (; i_vertex != i_end; ++i_vertex )
        {
          std::ostringstream sstr; 
          sstr << std::setw(12) << std::setprecision(7)
               << i_vertex->x << " "
               << i_vertex->y;
          is_comment ? printxmg.add_comment(sstr.str()) : printxmg.add_line( sstr.str() );
        }
        is_comment ? printxmg.add_comment( "&" ) : printxmg.add_line( "&" );
      }

      bool Restart( const TiXmlElement &_node )
      {
        const TiXmlElement *xmlresults = &_node;
        std::string name = _node.Value();
        if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
        if ( not xmlresults ) return false;
        LoadObject<t_Object, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_LONG);
        if( not convexhull.Load( *xmlresults, loadop ) ) return false;
        printxmg.add_comment("Loading ConvexHull from Input");
        typename t_ConvexHull :: t_Vertices :: const_iterator i_vertex = convexhull.begin_vertex();
        typename t_ConvexHull :: t_Vertices :: const_iterator i_end = convexhull.end_vertex();
        for (; i_vertex != i_end; ++i_vertex )
        {
          std::ostringstream sstr; 
          sstr << std::setw(12) << std::setprecision(7)
               << i_vertex->x << " "
               << i_vertex->y;
          printxmg.add_comment( sstr.str() );
        }
        printxmg.add_comment( "&" );
        return true;
      }
      void Save( TiXmlElement &_node ) const
      {
        TiXmlElement *parent = new TiXmlElement("Results");
        SaveObject<t_Object, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_LONG);
        convexhull.Save( *parent, saveop );
        _node.LinkEndChild(parent);
      }
#ifdef _MPI
      virtual bool broadcast( mpi::BroadCast &_bc )
      {
        if (    ( not convexhull.broadcast(_bc) ) 
             or ( not _bc.serialize(valid_ch) )    ) return false;
        return true;
      }

      // for each process,
      //  _ broadcasts individuals with fitness = 0
      //  _ add each individual to convex_hull
      virtual void synchronize( t_Population &_pop )
      {
        mpi::AllGather allgather( mpi::main );
        typename t_Population :: iterator i_indiv = _pop.begin();
        typename t_Population :: iterator i_end = _pop.end();
        for(; i_indiv != i_end; ++i_indiv )
          if( ( (types::t_real)i_indiv->fitness() ) < atat::zero_tolerance )
            i_indiv->broadcast(allgather);
        allgather.allocate_buffers();
        i_indiv = _pop.begin();
        for(; i_indiv != i_end; ++i_indiv )
          if( ( (types::t_real)i_indiv->fitness() ) < atat::zero_tolerance )
            i_indiv->broadcast(allgather);
        allgather();
        t_Individual indiv;
        while( indiv.broadcast(allgather) )
        {
          convexhull.add( indiv.get_quantity(), indiv );
          new_results = true;
        }
      }
#endif
  };
} // namespace darwin

#endif // _RESULTS_H_

//
//  Version: $Id$
//
#include<iostream>
#include<sstream>

#include <eo/eoDetTournamentSelect.h>
#include <eo/eoReduceMerge.h>
#include <eo/utils/eoRNG.h>

#include <pescan_interface/interface.h>
#include <print/stdout.h>
#include <print/xmg.h>
#include <print/manip.h>
#include <opt/debug.h>
#include <opt/tinyxml.h>

#include "functors.h"
#include "statistics.h"
#include "taboos/populations.h"

#include "debug.h"

namespace LaDa
{

  namespace GA
  {
#   define OPENXMLINPUT(name) \
      TiXmlDocument doc( name.string() ); \
      TiXmlHandle docHandle( &doc ); \
      __DOASSERT( not doc.LoadFile(), \
                     doc.ErrorDesc() << "\n" \
                  << "Could not load input file " << name \
                  << " in  Darwin<T_EVALUATOR>.\nAborting.\n" ) 
 
    template<class T_EVALUATOR>
    Darwin<T_EVALUATOR> :: ~Darwin()
    {
      if ( store ) delete store;
      if ( objective ) delete objective;
      if ( evaluation ) delete evaluation;
      if ( scaling ) delete scaling;
    }
    // reads in different parameters
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Load_Parameters(const TiXmlElement &_parent)
    {
      opt::explore_attributes( att_factory, _parent );
      // some checking
      if ( std::floor( pop_size * replacement_rate ) == 0 )
      {
        Print::xmg << Print::Xmg::comment
                   << "Error: replacement_rate is too small." << Print::endl 
                   << Print::Xmg::comment 
                   << "Error: setting replacement_rate too 1.0 / pop_size ."
                   << Print::endl;
        Print::out << "Error: replacement_rate is too small.\n"
                   << "Error: setting replacement_rate too 1.0 / pop_size .\n";
        replacement_rate = 1.00 / pop_size + 10*types::tolerance;
      }
 
      return true;
    }
 
    // creates history object if required
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: make_History(const TiXmlElement &_parent)
    {
      if( not topology.history() ) return;
      // checks if there are more than one taboo list
      const TiXmlElement *child = _parent.FirstChildElement("History");
      if ( not child ) return;
      Print::xmg << Print::Xmg::comment << "Track History" << Print::endl;
      history.is_on();
    }
    
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_Method(const TiXmlElement &_parent)
    {
      __TRYBEGIN
        if( topology.objective() ) 
          objective = topology.objective<t_GATraits>( _parent );
        __DOASSERT( topology.objective() and not objective, "Could not create multi-objective.\n" )
        
        evaluation = topology.evaluation< t_GATraits,
                                          Evaluation::WithHistory >( evaluator );
        {
          Evaluation::WithHistory<t_GATraits>*
            withhistory( static_cast< Evaluation::WithHistory<t_GATraits>* >( evaluation ) );
          withhistory->set( history );
        }
        Load_Storage( _parent );
        Evaluation::Base<t_GATraits>*
          base( static_cast< Evaluation::Base<t_GATraits>* >( evaluation ) );
        base->set( objective );
        base->set( store );
        
        if( not topology.scaling() ) return;
        scaling = Scaling::new_from_xml<t_GATraits>( _parent, &evaluator );
      __TRYEND(, "Could not create method.\n" )
    }
    
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_Storage(const TiXmlElement &_parent)
    {
      __TRYBEGIN
      if( not topology.store() ) return;
      const TiXmlElement *child = NULL, *store_xml = NULL;
      
      store = topology.special_store<t_GATraits>( evaluator );
      if( store ) goto endstorage;
      
      
      store_xml = _parent.FirstChildElement("Store");
      if ( store_xml )
      {
        store = new typename t_Store::FromObjective( evaluator, *store_xml );
        goto endstorage;
      }
 
      // Note that one of these tags have already been found and should exist
      child = _parent.FirstChildElement("Objective");
      if ( not child ) child = _parent.FirstChildElement("Objectives");
      if ( not child ) child = _parent.FirstChildElement("Method");
 
      store = new typename t_Store::FromObjective( evaluator, objective, *child );
      if ( store ) goto endstorage;
 
      store = new typename t_Store::Optima( evaluator, *child );
 
      endstorage:
        Print::xmg << Print::Xmg::comment << "Store: "
                   << store->what_is() << Print::endl;
      
      __TRYEND(, "Could not create storage interface.\n" )
    }
    
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: cleanup()
    {
      if( objective ) delete objective;
      if( store ) delete store;
      if( evaluation ) delete evaluation;
      if( scaling ) delete scaling;
      eostates.~eoState();
    }
 
    // create Taboos
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_Taboos( const TiXmlElement &_node )
    {
      if( not topology.taboos() ) return;
      if( not _node.FirstChildElement("Taboos") ) return;
      Taboo::Factory::create_container( taboo_factory, *_node.FirstChildElement("Taboos"), taboos );
      topology.special_taboos<t_GATraits>( taboos );
    }
 
    // Loads mating operations
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Load_Mating (const TiXmlElement &_parent)
    {
      if( not topology.mating() ) return true;
      const TiXmlElement *child = _parent.FirstChildElement("Breeding");
      if ( not child )
      {
        std::cerr << "No Breeding operator found in input. " << std::endl;
        std::cerr << "Giving Up" << std::endl;
        return false;
      }
      Print::xmg << Print::Xmg::comment << "Breeding Operator Begin" << Print::endl;
      breeder_ops = make_genetic_op( *child, breeder_ops );
      Print::xmg << Print::Xmg::comment << "Breeding Operator End" << Print::endl;
      __DOASSERT( not breeder_ops,
                 "Error while creating breeding operators")
      return true;
    }
    
    // Loads CheckPoints
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_CheckPoints (const TiXmlElement &_parent)
    {
      __MPICODE( if( not topology.continuators() ) return; )
      LaDa::Factory::visit_xml( checkpoint_factory, _parent, checkpoints );
    } // end of make_check_point
    
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Restart()
    {
      if( not topology.restart() ) return true;
      TiXmlDocument doc( restart_filename.string() );
      TiXmlHandle docHandle( &doc ); 
      if  ( !doc.LoadFile() ) 
      { 
        std::cerr << doc.ErrorDesc() << std::endl; 
        return false;
      } 
      TiXmlElement *parent = docHandle.FirstChild("Job").FirstChild("Restart").Element();
      if ( not parent ) return false;
 
      return Restart( *parent );
    }
 
    // Loads restarts operations
    // Should be in a Restart XML Section
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Restart (const TiXmlElement &_node)
    {
      if( not topology.restart() ) return true;
      const TiXmlElement *parent = &_node;
      std::string name = _node.Value();
      if ( name.compare("Restart") )
        parent = _node.FirstChildElement("Restart");
      if ( not parent ) return false;
     
      if ( do_restart & SAVE_POPULATION ) islands.clear();
      const TiXmlElement *child = parent->FirstChildElement();
      for(; child; child = child->NextSiblingElement() )
      {
        std::string name = child->Value();
        if (     name.compare("Results") == 0
             and ( do_restart & SAVE_RESULTS ) )
        {
          store->Restart(*child);
        }
        else if (     name.compare("History") == 0
                  and ( do_restart & SAVE_HISTORY ) and history.is_on() )
        {
          LoadObject<t_GATraits> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
          history.off(); // clears history.
          history.on();
          if ( LoadIndividuals( *child, loadop, history) )
            Print::xmg << Print::Xmg::comment << "Reloaded " << history.size() 
                       << " individuals into history" << Print::endl;
          
        } // end of "Population" Restart
        else if (     name.compare("Population") == 0
                  and ( do_restart & SAVE_POPULATION ) )
        {
          if( islands.size() > nb_islands ) continue;
          LoadObject<t_GATraits> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
          t_Population pop;
          LoadIndividuals( *child, loadop, pop, pop_size);
          if( pop.size() == 0 ) continue;
          islands.push_back( pop );
          Print::xmg << Print::Xmg::comment << "Loaded " << pop.size() 
                     << " individuals into island " << islands.size()
                     << Print::endl;
          Print::out << "Loaded " << pop.size() 
                     << " individuals into island " << islands.size()
                     << Print::endl;
          
          
        } // end of "Population" Restart
      }
      return true;
    }
 
    template<class T_EVALUATOR>
    eoReplacement<typename T_EVALUATOR::t_GATraits::t_Individual>*
      Darwin<T_EVALUATOR> :: make_replacement()
    {
      if (not topology.replacement() ) return NULL;
      eoTruncate<t_Individual>* truncate       = new  eoTruncate<t_Individual>;
      eoMerge<t_Individual>* merge             = new  eoPlus<t_Individual>;
      eoReduceMerge<t_Individual>* reducemerge
        = new eoReduceMerge<t_Individual>( *truncate, *merge );
      eostates.storeFunctor(truncate);
      eostates.storeFunctor(merge);
      eostates.storeFunctor(reducemerge);
 
      return reducemerge;
    }
 
    // makes genetic operators
    template<class T_EVALUATOR>
    eoGenOp<typename T_EVALUATOR::t_GATraits::t_Individual>*
      Darwin<T_EVALUATOR> :: make_genetic_op( const TiXmlElement &el,
                                              eoGenOp<t_Individual> *current_op )
    {
      return Operator::create_eo_operator( el, operator_factory, eostates );
    }
 
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: make_breeder()
    {
      eoSelectOne<t_Individual> *select;
 
      breeder = topology.breeder<t_GATraits>( evaluator ); 
      if( not breeder ) return;
      select = new eoDetTournamentSelect<t_Individual>(tournament_size);
      breeder->set(replacement_rate);
      breeder->set(counter);
      breeder->set(select);
      breeder->set(breeder_ops);
      topology.set<t_GATraits>( breeder, &taboos );
 
      eostates.storeFunctor(breeder);
      eostates.storeFunctor(select);
    }
    
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: populate ()
    {
      if( not topology.populate() ) 
      {
        islands.resize( nb_islands );
        return;
      }
      islands.resize( nb_islands );
      typename t_Islands :: iterator i_pop = islands.begin();
      typename t_Islands :: iterator i_end = islands.end();
      t_Object&(t_Individual::*ptr_func)( void ) = &t_Individual::Object;
      for(; i_pop != i_end; ++i_pop)
        population_creator( *i_pop, pop_size );
    }
 
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: random_populate ( t_Population &_pop,
                                                  types::t_unsigned _size)
    {
      types::t_unsigned i = 0, j = 0;
      types::t_unsigned maxiter = _size * 50;
      while ( _pop.size() < _size and i < maxiter)
      {
        t_Individual indiv;
        evaluator.initialize(indiv);
        if ( not ( taboos and (*taboos)(indiv) ) )
        {
          _pop.push_back(indiv);
          __DODEBUGCODE( Print::out << "Generated Individual " 
                                    << _pop.size() << " of " 
                                    << _size << " " << indiv << ".\n"; )
            ++j;
        }
        ++i;
      } // while ( i_pop->size() < target )
      __DOASSERT( j < _size,
                     "Created " << j << " individuals in " << i << " iterations.\n"
                  << "Are taboos/concentration constraints to restrictive?\n" )
      _pop.resize( _size );
    }
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: partition_populate ( t_Population &_pop,
                                                     types::t_unsigned _size)
    {
      while ( _pop.size() < _size )
      {
        // first random object
        t_Individual indiv;
        evaluator.initialize(indiv);
        types::t_unsigned start = 0;  
        types::t_unsigned end = indiv.Object().size();  
      
        // indiv: ----
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
        
        indiv.Object().mask( start, end );         // indiv: ****
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
        
        indiv.Object().mask( start, end/2 );       // indiv: --**
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end );         // indiv: **--
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end/4 );       // indiv: -*--
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end );         // indiv: *-**
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( end/2, end*3/4 );     // indiv: *--*
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end );         // indiv: -**-
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( end/4, end/2 );       // indiv: --*-
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end );         // indiv: **-*
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end/2 );       // indiv: ---*
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
        if ( _pop.size() >= _size ) break;
 
        indiv.Object().mask( start, end );         // indiv: *---
        if ( not ( taboos and (*taboos)(indiv) ) ) _pop.push_back(indiv);
      }
      _pop.resize( _size );
    }
 
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: presubmit()
    {
      // GA::Evaluator does not know about Traits::GA, and even less about
      // Traits::GA::t_Population, hence this bit of hackery.
      std::list< t_Individual > dummy;
      eoPop<t_Individual> dummy2;
      evaluator.presubmit(dummy);
      if ( dummy.empty() ) return;
      offspring.resize( dummy.size() );
      std::copy( dummy.begin(), dummy.end(), offspring.begin() );
      (*evaluation)(dummy2, offspring);
      offspring.clear();
      evaluation->nb_eval = 0;
      evaluation->nb_grad = 0;
    }
 
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: run()
    {
      presubmit();
 
      Print::xmg << Print::flush;
      populate();
 
      offspring.clear();
      typename t_Islands :: iterator i_island_begin = islands.begin();
      typename t_Islands :: iterator i_island_end = islands.end();
      typename t_Islands :: iterator i_island;
      __MPICODE
      ( 
        boost::mpi::communicator world;
        world.barrier(); 
      )
      // A first eval of pop.
      for ( i_island = i_island_begin; i_island != i_island_end; ++i_island )
        __TRYDEBUGCODE( (*evaluation)(offspring, *i_island);,
                        "Error while evaluating starting population\n" )
      if( do_starting_population_only ) 
      { 
        checkpoints(islands);
        Print::xmg << "Performing evaluation of starting population only.\n" << Print::flush;
        return;
      }
      types::t_unsigned n = 0;
 
      Print::out << "\nEntering Generational Loop" << Print::endl;
      do
      {
        Print::out << "Iteration " << n << Print::endl;
        Print::xmg << Print::flush;
        ++n;
        i_island = i_island_begin;
        for (; i_island != i_island_end; ++i_island )
        {
          topology.syncpop( *i_island );
          types::t_unsigned pSize = i_island->size();
          offspring.clear(); // new offspring
          
          __DODEBUGCODE( Print::out << "Scaling prior to breeding" << Print::endl; )
          if( scaling ) (*scaling)( *i_island );
 
          __DODEBUGCODE( Print::out << "Breeding" << Print::endl; )
 
          if( breeder ) (*breeder)(*i_island, offspring);
          
          __DODEBUGCODE( Print::out << "Evaluating" << Print::endl; )
 
          // eval of parents + offspring if necessary
          __TRYDEBUGCODE( (*evaluation)(*i_island, offspring);,
                          "Error while Evaluating generation " << n << "\n" )
 
          __DODEBUGCODE(Print::out << "Scaling prior to replacing" << Print::endl;)
          // Does scaling simultaneously over both populations
          if( scaling ) (*scaling)( *i_island, offspring );
 
          __DODEBUGCODE(Print::out << "Replacing" << Print::endl;)
          // after replace, the new pop. is in population
          if( replacement ) (*replacement)(*i_island, offspring); 
          
          __ASSERT( pSize != i_island->size(),
                       "Population "
                    << ( pSize > i_island->size() ? "shrinking": "growing" )
                    << " from " << pSize << " to " << i_island->size() << "\n"
                  )
 
          __DODEBUGCODE(Print::out << "Continuing" << Print::endl;)
        }
      } while ( checkpoints( islands ) ); 
 
      Print::xmg << Print::flush;
    }
 
 
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Load(const t_Path &_filename) 
    {
      // Initializes each proc differently
      filename = _filename;
      TiXmlDocument doc;
      TiXmlHandle docHandle( &doc ); 
      opt::read_xmlfile( filename, doc );
      const TiXmlElement *node;

      // Loads topology and assigns comms to evaluators
      node = docHandle.FirstChild("Job").Element();
      topology.Load( *node, evaluator );

      // Loads evaluator first 
      __TRYASSERT( not evaluator.Load(*docHandle.FirstChild("Job").Element() ),
                      "Could not load functional input from "
                   << filename << "\n" )
 
      // Load checkpoints and filenames.
      node = docHandle.FirstChild("Job")
                      .FirstChild("GA").Element();
      Load_CheckPoints( *node );
          
      // finds <GA> ... </GA> block 
      __TRYASSERT( not Load_Parameters( *node ), 
                   "Error while reading GA attributes\n" )
 
 
 
      make_History( *node );
      Load_Taboos( *node );
      Load_Method( *node );
      __DOASSERT( not  Load_Mating( *node ),
                  "Error while loading mating operators.\n" )
      
      make_breeder();
      replacement = make_replacement();
 
      const TiXmlElement *restart_xml = opt::find_node( *node, "Filenames", "restart" );
      if ( restart_xml )
      {
        std::string str = "all";
        if ( restart_xml->Attribute("what") ) str = restart_xml->Attribute("what");
        if ( str.find("all") != std::string::npos )
          do_restart |= SAVE_POPULATION | SAVE_HISTORY | SAVE_RESULTS;
        else
        {
          if ( str.find("pop") != std::string::npos )
            do_restart |= SAVE_POPULATION;
          if ( str.find("history") != std::string::npos and history.is_on() )
            do_restart |= SAVE_HISTORY;
          if ( str.find("results") != std::string::npos )
            do_restart |= SAVE_RESULTS;
        }
      }
      if ( do_restart )
      {
        if ( restart_filename != filename )
        { 
          if ( not Restart() )
          {
            std::cerr << "Could not load restart from file" << restart_filename;
            Print::xmg << Print::Xmg :: comment << "Starting from scratch" << Print::endl;
            Print::out << "Could not load restart from file" << restart_filename << "\n"
                       << "Starting from scratch\n";
          }
        } // for following, we know docHandle.FirstChild("Job") exists
        else if ( not Restart( *docHandle.FirstChild("Job").Element() ) )
        {
          std::cerr << "Could not load restart from file" << restart_filename;
          Print::xmg << Print::Xmg :: comment << "Starting from scratch" << Print::endl;
          Print::out << "Could not load restart from file" << restart_filename << "\n"
                     << "Starting from scratch\n";
        }
      }
      Print::xmg << Print::endl;
 
      // reseeds accordint to mpi/serial topology
      topology.reseed();
 
      Print::xmg << Print::Xmg::comment << "Mating Tournament Size: "
                                        << tournament_size << Print::endl
                 << Print::Xmg::comment << "Offspring Replacement Rate: "
                                        << replacement_rate << Print::endl
                 << Print::Xmg::comment << "Population Size: "
                                        << pop_size << Print::endl
                 << Print::Xmg::comment << topology.print_seeds() << Print::endl 
                                        << Print::endl
                 << Print::Xmg::comment << topology.print() << Print::endl;
      if ( max_generations )
        Print::xmg << Print::Xmg::comment << "Maximum Number of Generations: "
                                          << max_generations << Print::endl;
      else
        Print::xmg << Print::Xmg::comment << "Unlimited Number of Generations"
                                          << Print::endl;
      Print::xmg << Print::Xmg::comment << "Number of Islands: "
                                        << nb_islands << Print::endl;
      if( scaling ) Print::xmg << Print::Xmg::comment << scaling->what_is() << Print::endl;
      std::string evalparams = evaluator.print();
      Print::xmg << Print::make_commented_string( evalparams ) << Print::endl;
      Print::xmg << Print::flush;
      return true;
    }
  } // namespace GA
} // namespace LaDa

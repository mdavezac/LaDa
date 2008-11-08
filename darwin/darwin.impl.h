//
//  Version: $Id$
//
#include<iostream>
#include<sstream>

#include <eo/eoDetTournamentSelect.h>
#include <eo/eoReduceMerge.h>
#include <eo/utils/eoRNG.h>
#include <tinyxml/tinystr.h>

#include <pescan_interface/interface.h>
#include <print/stdout.h>
#include <print/xmg.h>
#include <print/manip.h>
#include <opt/debug.h>

#include "functors.h"
#include "statistics.h"
#include "minimizergenop.h"

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
      // tournament size when selecting parents
      const TiXmlAttribute *att = _parent.FirstAttribute();
      for(; att; att = att->Next() )
      {
        std::string str = att->Name();
        
        if ( str.compare("tournament")==0 )
        {
          tournament_size = std::abs(att->IntValue());
          if ( tournament_size < 2 ) tournament_size = 2;
        }
        else if ( str.compare("rate")==0 ) // offspring rate
        {
          replacement_rate = std::abs(att->DoubleValue());
          if ( replacement_rate > 1 ) replacement_rate = 0.5;
        }
        else if ( str.compare("popsize")==0 ) // population size
        {
          pop_size = std::abs(att->IntValue());
          if ( pop_size < 1 ) pop_size = 1;
          
        }
        else if ( str.compare("maxgen")==0 ) // maximum number of generations
          max_generations = std::abs(att->IntValue());
        else if ( str.compare("islands")==0 ) // number of islands
        {
          nb_islands = std::abs(att->IntValue());
          if ( not nb_islands ) nb_islands = 1;
        }
        else if ( str.compare("print")==0 ) // print at each generation
          do_print_each_call = true;
        else if ( str.compare("populate")==0 ) 
        {
          std::string string = att->Value();
          populate_style = RANDOM_POPULATE;
          if ( string.compare("partition") == 0 )
            populate_style = PARTITION_POPULATE;
        }
        else if ( not topology.LoadSeeds( *att ) )
          evaluator.LoadAttribute( *att ); 
      }
 
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
      history = topology.history<t_GATraits>( eostates );
    }
    
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_Method(const TiXmlElement &_parent)
    {
      __TRYBEGIN
        if( topology.objective() ) 
          objective = topology.objective<t_GATraits>( _parent );
        __DOASSERT( not objective, "Could not create multi-objective.\n" )
        
        if( topology.history() and history )
        {
          evaluation = topology.evaluation< t_GATraits,
                                            Evaluation::WithHistory >( evaluator );
          Evaluation::WithHistory<t_GATraits>*
            withhistory( static_cast< Evaluation::WithHistory<t_GATraits>* >( evaluation ) );
          withhistory->set( history );
        }
        if ( not evaluation )
          evaluation = topology.evaluation<t_GATraits, Evaluation::Base >( evaluator );
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
      // checks if there are more than one taboo list
      if ( not _node.FirstChildElement("Taboos") ) return;
      const TiXmlElement &parent = *_node.FirstChildElement("Taboos");
      const TiXmlElement *child = parent.FirstChildElement();
 
      // Checks for topology specifics
      if( not topology.taboos() ) return;
      taboos = topology.special_taboos<t_GATraits>( eostates );
      if( taboos ) return;
 
      // creates Taboo container if there are more than one taboo list
      taboos = new Taboos<t_Individual>;
 
      // creates pop taboo
      child = parent.FirstChildElement("Population");
      if ( not child ) child = parent.FirstChildElement("PopTaboo");
      if (child)
      {
        Print::xmg << Print::Xmg::comment << "Population Taboo" << Print::endl;
        IslandsTaboos<t_GATraits> *poptaboo
           = new IslandsTaboos<t_GATraits>( islands );
        eostates.storeFunctor(poptaboo);
        static_cast< Taboos<t_Individual>* >(taboos)->add( poptaboo );
      }
      
      // creates offspring taboo
      child = parent.FirstChildElement("Offspring");
      if ( not child ) child = parent.FirstChildElement("OffspringTaboo");
      if (child)
      {
        Print::xmg << Print::Xmg::comment << "Offspring Taboo" << Print::endl;
        OffspringTaboo<t_GATraits> *offspringtaboo 
           = new OffspringTaboo<t_GATraits>( &offspring );
        eostates.storeFunctor(offspringtaboo);
        static_cast< Taboos<t_Individual>* >(taboos)->add( offspringtaboo );
      }
 
      // makes history a taboo list if it exists
      child = parent.FirstChildElement("History");
      if (child and history)
      {
        Print::xmg << Print::Xmg::comment << "History Taboo" << Print::endl;
        static_cast< Taboos<t_Individual>* >(taboos)->add( history );
      }
      else if (child)
        std::cerr << "HistoryTaboo found in Taboos tags, but not History tag found!!" << std::endl
                  << "Include History tag if you want HistoryTaboo" << std::endl;
 
 
      // then evaluator specific taboos
      child = parent.FirstChildElement();
      for( ; child ; child = child->NextSiblingElement() )
      {
        Taboo_Base<t_Individual> *func = evaluator.LoadTaboo( *child );
        eostates.storeFunctor(func);
        static_cast< Taboos<t_Individual>* >(taboos)->add( func );
      }
 
      // now some cleaning up
      Taboos<t_Individual>* save_taboos = static_cast< Taboos<t_Individual>* >(taboos);
      types::t_unsigned n = save_taboos->size();
      switch( n )
      {
        case 0: delete taboos; taboos = NULL; break; // no taboos! 
        case 1: taboos = save_taboos->front(); delete save_taboos; break; // move lone taboo to front
        default: eostates.storeFunctor( taboos ); break; // multiple taboos, good as is
      }
        
      return;
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
      breeder_ops = new SequentialOp<t_GATraits>;
      eostates.storeFunctor( breeder_ops );
      breeder_ops = make_genetic_op(*child->FirstChildElement(), breeder_ops);
      Print::xmg << Print::Xmg::comment << "Breeding Operator End" << Print::endl;
      __DOASSERT( not breeder_ops,
                 "Error while creating breeding operators")
      return true;
    }
    
    // Loads CheckPoints
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: Load_CheckPoints (const TiXmlElement &_parent)
    {
      // contains all checkpoints
      std::string str = "stop"; 
      if (      _parent.FirstChildElement("Filenames") 
           and  _parent.FirstChildElement("Filenames")->Attribute("stop") )
        str = _parent.FirstChildElement("Filenames")->Attribute("stop");
 
      str = Print::reformat_home(str);
      continuator = new IslandsContinuator<t_GATraits>(max_generations, str );
      eostates.storeFunctor( continuator );
 
      __MPICODE( if( not topology.continuators() ) return; )
 
      GenCount &generation_counter = continuator->get_generation_counter();
 
      // Creates SaveEvery object if necessary
      if (      _parent.FirstChildElement("Save") 
           and  _parent.FirstChildElement("Save")->Attribute("every") )
      {
        types::t_int n;
        _parent.FirstChildElement("Save")->Attribute("every", &n);
        
        if ( n >= 0 and do_save )
        {
          Print::xmg << Print::Xmg::comment << "Will Save Every " << n
                     << " Generations " << Print::endl;
          SaveEvery<t_This> *save = new SaveEvery<t_This>( *this,
                                                           &Darwin::Save, 
                                                           std::abs(n) );
          eostates.storeFunctor( save );
          continuator->add( *save );
        }
      }
   
      // Creates Statistics -- only one type available...
      const TiXmlElement *child = _parent.FirstChildElement("Statistics");
      for(; child; child = child->NextSiblingElement("Statistics") )
      {
        if( not child->Attribute("type") ) continue;
        std::string name = child->Attribute("type");
        if( name.compare("census") == 0 )
        {
          continuator->add( eostates.storeFunctor( new TrueCensus< t_GATraits >() ) );
          Print::xmg << Print::Xmg::comment
                     << "Statistics: True population size, discounting twins"
                     << Print::endl;
        }
        else if( name.compare("AverageFitness") == 0 )
        {
          continuator->add( eostates.storeFunctor( new AverageFitness< t_GATraits >() ) );
          Print::xmg << Print::Xmg::comment
                     << "Statistics: Average Fitness" << Print::endl;
        }
        else if( name.compare("AverageQuantity") == 0 )
        {
          continuator->add(
              eostates.storeFunctor( new AverageQuantities< t_GATraits >() ) );
          Print::xmg << Print::Xmg::comment
                     << "Statistics: Average Quantity " << Print::endl;
        }
      }
 
      // Creates Terminators
      child = _parent.FirstChildElement("Terminator");
      for( ; child; child = child->NextSiblingElement("Terminator") )
      {
        eoContinue<t_Individual> *terminator = NULL;
        std::string type = "<";
        std::string ref = "";
        if ( child->Attribute("ref" ) )
          ref = child->Attribute( "ref" );
        int max = 0;
        child->Attribute("value", &max);
 
        if ( max <= 0 ) continue;
        if ( type.compare("<") != 0 ) continue;
 
        if ( ref.compare("evaluation") != 0 ) continue;
        
        terminator = new Terminator< types::t_unsigned,
                                     std::less<types::t_unsigned>, t_GATraits >
                                   ( evaluation->nb_eval, (types::t_unsigned) abs(max),
                                     std::less<types::t_unsigned>(), "nb_eval < term" );
        eostates.storeFunctor( terminator );
        continuator->add( *terminator );
        Print::xmg << Print::Xmg::comment << "Terminating after " << max
                   << " evaluations" << Print::endl;
        
        // end if max
      }
      
      // Load object specific continuator
      eoF<bool> *specific = evaluator.LoadContinue( _parent );
      if ( specific )
      {
        eostates.storeFunctor(specific); 
        continuator->add(
            eostates.storeFunctor(new Continuator<t_GATraits>(*specific)) );
      }
 
      // Creates Print object
      {
        typedef PrintGA< Store::Base<t_GATraits>,
                         Evaluation::Abstract<t_Population> > t_PrintGA;
        t_PrintGA* printga = new t_PrintGA( *store, *evaluation,
                                            generation_counter, do_print_each_call);
        eostates.storeFunctor(printga);
        continuator->add( *printga );
      }
 
      // Print Offsprings
      child = _parent.FirstChildElement("Print"); 
      for(; child; child = child->NextSiblingElement("Print") )
      {
        if( not child->Attribute("type") ) continue;
        std::string name = child->Attribute("type");
        if( name.compare("offspring") == 0 )
        {
          continuator->add(
              eostates.storeFunctor( new PrintOffspring< t_GATraits >
                                                       ( generation_counter ) ) );
          Print::xmg << Print::Xmg::comment
                     << "Print: offspring" << Print::endl;
        }
        else if(    name.compare("pop") == 0 
                 or name.compare("population") == 0 )
        {
          continuator->add( eostates.storeFunctor( new PrintPop< t_GATraits >() ) );
          Print::xmg << Print::Xmg::comment
                     << "Print: current population" << Print::endl;
        }
        else
        {
          eoMonOp<const t_Individual> *op =  evaluator.LoadPrintBest( *child );
          if ( not op ) continue;
          eostates.storeFunctor( op );
          Apply2Best<t_GATraits> *printbest = new Apply2Best<t_GATraits>( *store );
          try{ printbest = new Apply2Best<t_GATraits>( *store ); }
          catch(...) { cleanup(); __THROW_ERROR( "Memory allocation error.\n" ) }
          __DOASSERT(not printbest,  "Memory allocation error.\n");
          printbest->set_functor( op );
          continuator->add(*printbest);
        }
      }
 
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
      TiXmlElement *parent = docHandle.FirstChild("Restart").Element();
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
      if ( not parent )
        return false;
     
      const TiXmlElement *child = parent->FirstChildElement();
      for(; child; child = child->NextSiblingElement() )
      {
        std::string name = child->Value();
        if (     name.compare("Results") == 0
             and ( do_restart & SAVE_RESULTS ) )
          store->Restart(*child);
        else if (     name.compare("History") == 0
                  and ( do_restart & SAVE_HISTORY ) and history )
        {
          LoadObject<t_GATraits> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
          history->clear();
          if ( LoadIndividuals( *child, loadop, *history) )
            Print::xmg << Print::Xmg::comment << "Reloaded " << history->size() 
                       << " individuals into history" << Print::endl;
          
        } // end of "Population" Restart
        else if (     name.compare("Population") == 0
                  and ( do_restart & SAVE_POPULATION ) )
        {
          LoadObject<t_GATraits> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
          islands.clear();
          const TiXmlElement *island_xml = child->FirstChildElement("Island");
          for(; island_xml; island_xml = island_xml->NextSiblingElement("Island") )
          {
            if( islands.size() > nb_islands )
              break;
            t_Population pop;
            LoadIndividuals( *island_xml, loadop, pop, pop_size);
            islands.push_back( pop );
            Print::xmg << Print::Xmg::comment << "Loaded " << pop.size() 
                       << " individuals into island " << islands.size()
                       << Print::endl;
          }
        } // end of "Population" Restart
      }
      return true;
    }
 
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Save()
    {
      if ( not do_save ) return true;
      if ( not topology.save() ) return true;
 
      TiXmlDocument doc;
      TiXmlElement *node;
      if ( save_filename == filename )
      {
        doc.LoadFile(save_filename.string());
        TiXmlHandle docHandle( &doc );
        node = docHandle.FirstChild("Job").Element();
        TiXmlNode* child = docHandle.FirstChild("Job")
                                    .FirstChild("Restart").Node();
        for(; child; child = docHandle.FirstChild("Job").FirstChild("Restart").Node() )
          node->RemoveChild(child);
        node = new TiXmlElement("Restart");
        if ( not node )
          return false;
        child = docHandle.FirstChild("Job").Element();
        child->LinkEndChild( node );
      }
      else
      {
        doc.SetTabSize(1);
        doc.LinkEndChild( new TiXmlDeclaration("1.0", "", "") );
        node = new TiXmlElement("Restart");
        if ( not node )
          return false;
        doc.LinkEndChild( node );
      }
 
      Print::xmg << Print::Xmg::comment <<  "Saving ";
      int is_saving = 0;
      if ( do_save & SAVE_RESULTS )
      {
        ++is_saving;
        Print::xmg << "Results";
        store->Save( *node );
      }
      if ( do_save & SAVE_HISTORY and history)
      {
        // Printout stuff
        ++is_saving;
        if ( is_saving > 1 )
          Print::xmg << ", ";
        Print::xmg << " History";
 
        TiXmlElement *xmlhistory = new TiXmlElement("History");
        SaveObject<t_GATraits> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
        SaveIndividuals( *xmlhistory, saveop, history->begin(), history->end());
        node->LinkEndChild( xmlhistory );
      }
      if ( do_save & SAVE_POPULATION )
      {
        // Printout stuff
        ++is_saving;
        if ( is_saving == 2 )
          Print::xmg << " and";
        if ( is_saving > 2 )
          Print::xmg << ", and ";
        Print::xmg << " Population";
 
        SaveObject<t_GATraits> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
        typename t_Islands :: const_iterator i_island = islands.begin();
        typename t_Islands :: const_iterator i_island_end = islands.end();
        TiXmlElement *xmlpop = new TiXmlElement("Population");
        for(; i_island != i_island_end; ++i_island )
        {
          TiXmlElement *xmlisland = new TiXmlElement("Island");
          SaveIndividuals( *xmlisland, saveop, i_island->begin(), i_island->end());
          xmlpop->LinkEndChild( xmlisland );
        }
        node->LinkEndChild( xmlpop );
      }
 
      if ( not doc.SaveFile(save_filename.string() ) )
      {
        std::cerr << "Could not save results in " << save_filename << std::endl;
        Print::xmg << Print::Xmg::clear;
        Print::out << "********** Save failed!! \n" << Print::endl;
        return false;
      }
 
      Print::xmg << " in " << save_filename << Print::endl;
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
      eoOp<t_Individual>* this_op;
      const TiXmlElement *sibling = &el;
      if (not sibling) return NULL;
   
      Print::xmg << Print::Xmg::indent;
 
      for ( ; sibling; sibling = sibling->NextSiblingElement() )
      {
        std::string str = sibling->Value();
        double prob = 0.0;
        int period = 0;
        this_op = NULL;
        bool is_gen_op = false;
        
        // then creates sibling
        if ( str.compare("TabooOp") == 0  and taboos )
        {
          Print::xmg << Print::Xmg::comment << "TabooOp Begin" << Print::endl;
          eoGenOp<t_Individual> *taboo_op
            = make_genetic_op( *sibling->FirstChildElement(), NULL);
          if ( not taboo_op )
            Print::xmg << Print::Xmg::removelast;
          else
          {
            eoMonOp<t_Individual> *utterrandom = NULL;
            try
            {
              Print::xmg << Print::Xmg::comment << "TabooOp End" << Print::endl;
              utterrandom = new mem_monop_t<t_GATraits>
                                  ( evaluator, &t_Evaluator::initialize,
                                    std::string( "Initialize" ) );
              __DOASSERT( not utterrandom, "Memory Allocation error.\n" )
              eostates.storeFunctor(utterrandom);
              this_op = new TabooOp<t_Individual> ( *taboo_op, *taboos, 
                                                    (types::t_unsigned)(pop_size+1),
                                                    *utterrandom );
              __DOASSERT( not this_op, "Memory Allocation error.\n" )
              eostates.storeFunctor( static_cast< TabooOp<t_Individual> *>(this_op) );
              is_gen_op = true;
            }
            __CATCHCODE( if( utterrandom) delete utterrandom; cleanup(),
                         "Error encountered while creating Taboo genetic operator.\n")
          }
        }
        else if ( str.compare("TabooOp") == 0 )
        {
          Print::xmg << Print::Xmg :: unindent;
          this_op = make_genetic_op( *sibling->FirstChildElement(), NULL);
          Print::xmg << Print::Xmg :: indent;
        }
        else if ( str.compare("UtterRandom") == 0 )
        {
          Print::xmg << Print::Xmg::comment << "UtterRandom" << Print::endl;
          this_op = new mem_monop_t<t_GATraits>
                           ( evaluator, &t_Evaluator::initialize,
                             std::string( "UtterRandom" ) );
          eostates.storeFunctor( static_cast< TabooOp<t_Individual> *>(this_op) );
        }
        else if ( str.compare("Minimizer") == 0 )
        {
          Minimizer_Functional<t_GATraits> *func = NULL;
          MinimizerGenOp<t_GATraits> *mingenop = NULL;
          try
          {
            func =  new Minimizer_Functional<t_GATraits>( *evaluation, *taboos ); 
            __DOASSERT( not func, "Memory Allocation Error.\n")
          
            mingenop = new MinimizerGenOp<t_GATraits>( *func ); 
            __DOASSERT( not mingenop, "Memory Allocation Error.\n")
 
            if ( not mingenop->Load( *sibling ) ) 
              { delete mingenop; this_op = NULL; }
            else
            {
              eostates.storeFunctor( mingenop );
              this_op = mingenop;
            }
          }
          catch( std::exception &_e )
          {
            cleanup(); 
            __THROW_ERROR( "Could  not load minimizer genetic operator.\n" << _e.what() )
          }
          catch(...)
          { 
            cleanup();
            this_op = NULL;
            __THROW_ERROR( "Memory Allocation Error.\n" ) 
          }
        }
        else if ( str.compare("Operators") == 0 )
        {
          if ( sibling->Attribute("type") )
          {
            std::string sstr = sibling->Attribute("type");
            if ( sstr.compare("and") == 0 ) 
            {
              Print::xmg << Print::Xmg::comment << "And Begin" << Print::endl;
              SequentialOp<t_GATraits> *new_branch = new SequentialOp<t_GATraits>;
              this_op = make_genetic_op( *sibling->FirstChildElement(), new_branch);
              if ( not this_op )
              {
                Print::xmg << Print::Xmg::removelast;
                Print::out << " Failure, or No Operators Found"
                           << " in \"And\" operator \n";
                delete new_branch;
              }
              else
              {
                Print::xmg << Print::Xmg::comment << "And End" << Print::endl;
                eostates.storeFunctor( new_branch );
              }
            }
          }
          if ( not this_op )
          {
            Print::xmg << Print::Xmg::comment << "Or Begin" << Print::endl;
            ProportionalOp<t_GATraits> *new_branch = new ProportionalOp<t_GATraits>;
            this_op = make_genetic_op( *sibling->FirstChildElement(), new_branch);
            if ( not this_op )
            {
              Print::xmg << Print::Xmg::removelast;
              delete new_branch;
            }
            else
            {
              Print::xmg << Print::Xmg::comment << "Or End" << Print::endl;
              eostates.storeFunctor( new_branch );
            }
          }
          is_gen_op = true;
        }
        else // operators specific to t_Evaluator 
        {
          eoGenOp<t_Individual> *eoop = evaluator.LoadGaOp( *sibling );
          eostates.storeFunctor( eoop );
          if ( eoop ) this_op = eoop;
          is_gen_op = true;
        }
        if ( this_op and sibling->Attribute("period", &period) )
        {
          if (     (    ( not max_generations )
                     or (types::t_unsigned)std::abs(period) < max_generations ) 
               and period > 0 )
          {
            Print::xmg << Print::Xmg::addtolast << " period = "
                       << period << Print::endl;
            this_op = new PeriodicOp<t_GATraits>( *this_op,
                                                  (types::t_unsigned) abs(period),
                                                   continuator 
                                                      ->get_generation_counter(),
                                                    eostates );
            eostates.storeFunctor( static_cast< PeriodicOp<t_GATraits> *>(this_op) );
            is_gen_op = true;
          }
        }
        if ( this_op and current_op )
        {
          if (not sibling->Attribute("prob", &prob) )
            prob = 1.0;
          Print::xmg << Print::Xmg::addtolast << " prob = " << prob << Print::endl;
          if ( current_op->className().compare("GA::SequentialOp") == 0 )
            static_cast< SequentialOp<t_GATraits>* >(current_op)
               ->add( *this_op, static_cast<double>(prob) );
          else if ( current_op->className().compare("GA::ProportionalOp") == 0 )
            static_cast< ProportionalOp<t_GATraits>* >(current_op)
              ->add( *this_op, static_cast<double>(prob) );
        }
        else if ( this_op )
        {
          if ( is_gen_op )
            current_op = static_cast<eoGenOp<t_Individual>*> (this_op);
          else 
            current_op = &wrap_op<t_Individual>(*this_op, eostates);
        }
      }
      
      if ( not current_op  )
      {
        std::cerr << " Error while creating genetic operators " << std::endl;
        throw;
      } 
 
      Print::xmg << Print::Xmg :: unindent;
      
      return current_op;
    }
 
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: make_breeder()
    {
      eoSelectOne<t_Individual> *select;
 
      breeder = topology.breeder<t_GATraits>( evaluator ); 
      if( not breeder ) return;
      select = new eoDetTournamentSelect<t_Individual>(tournament_size);
      breeder->set(replacement_rate);
      breeder->set(&continuator->get_generation_counter());
      breeder->set(select);
      breeder->set(breeder_ops);
      topology.set<t_GATraits>( breeder, taboos );
 
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
      types::t_unsigned target = pop_size;
      for(; i_pop != i_end; ++i_pop)
      {
        switch ( populate_style )
        {
          case RANDOM_POPULATE:
            random_populate(*i_pop, target); break;
          case PARTITION_POPULATE:
            partition_populate(*i_pop, target); break;
        }
      }
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
      __MPICODE( ::LaDa::mpi::main->barrier(); )
      // A first eval of pop.
      for ( i_island = i_island_begin; i_island != i_island_end; ++i_island )
        __TRYDEBUGCODE( (*evaluation)(offspring, *i_island);,
                        "Error while evaluating starting population\n" )
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
      } while ( continuator->apply( i_island_begin, i_island_end ) );
 
      Print::xmg << Print::flush;
      Save();
    }
 
 
#  ifdef _MPI
    template<class T_EVALUATOR>
    void Darwin<T_EVALUATOR> :: LoadAllInputFiles(std::string &_input, 
                                                  std::string &_restart, 
                                                  std::string &_evaluator ) 
    {
      __NOTMPIROOT( (*::LaDa::mpi::main),
        boost::mpi::broadcast( *::LaDa::mpi::main, _input, 0 );
        boost::mpi::broadcast( *::LaDa::mpi::main, _restart, 0 );
        boost::mpi::broadcast( *::LaDa::mpi::main, _evaluator, 0 );
        return;
      )
 
      OPENXMLINPUT(filename)
      std::ostringstream stream;
      const TiXmlElement *parent = doc.RootElement();
      stream << *parent;
      _input = stream.str();
      _restart = _input; 
      _evaluator = _input;
 
      if ( restart_filename == filename )
      {
        doc.LoadFile( restart_filename.string() );
        if  ( !doc.LoadFile() ) 
        { 
          std::cerr << __SPOT_ERROR
                    << doc.ErrorDesc() << "\n"
                    << "Could not load restart file.\n"
                    << "Will ignore restart input and start run from scratch"
                    << std::endl;
          _restart = "";
          goto nextfilename;
        } 
 
        parent = doc.RootElement();
        std::ostringstream stream;
        stream << *parent;
        _restart = stream.str();
      }
      nextfilename:
        if ( evaluator_filename != filename )
        {
          doc.LoadFile( evaluator_filename.string() );
          __DOASSERT( not doc.LoadFile(), " Could not load restart file\n" )
          parent = doc.RootElement();
          std::ostringstream stream;
          stream << *parent;
          _evaluator = stream.str();
        }
       
        boost::mpi::broadcast( *::LaDa::mpi::main, _input, 0 );
        boost::mpi::broadcast( *::LaDa::mpi::main, _restart, 0 );
        boost::mpi::broadcast( *::LaDa::mpi::main, _evaluator, 0 );
    }
#  endif
 
    template<class T_EVALUATOR>
    bool Darwin<T_EVALUATOR> :: Load(const t_Path &_filename) 
    {
      // Initializes each proc differently
      filename = _filename;
      evaluator_filename = _filename;
      restart_filename = _filename;
      save_filename = _filename;
 
      t_Path xmg_filename = "out.xmg";
      t_Path out_filename = "out";
 
      TiXmlDocument doc( filename.string() ); 
      TiXmlHandle docHandle( &doc ); 
      const TiXmlElement *parent, *child;
 
      __NOTMPIROOT( (*::LaDa::mpi::main), goto syncfilenames; )
 
      // first loads all inputfiles 
      __DOASSERT( not doc.LoadFile(), 
                     doc.ErrorDesc() << "\n" 
                  << "Could not load input file " << filename 
                  << " in  Darwin<T_EVALUATOR>.\nAborting.\n" ) 
 
      parent = docHandle.FirstChild("Job").Element();
      __DOASSERT( not parent, 
                  "Could not find <Job> tag in  " << filename << "\n"; )
      parent = docHandle.FirstChild("Job")
                        .FirstChild("GA").Element();
      __DOASSERT( not parent, 
                  "Could not find <GA> tag in  " << filename << "\n"; )
      child = parent->FirstChildElement("Filenames");
      for( ; child; child = child->NextSiblingElement("Filenames") )
      {
        if (     child->Attribute("evaluator") 
             and evaluator_filename == filename )
          evaluator_filename = Print::reformat_home(child->Attribute("evaluator"));
        if (     child->Attribute("restart") 
                  and restart_filename == filename )
          restart_filename = Print::reformat_home(child->Attribute("restart"));
        if (     child->Attribute("save") 
                  and save_filename == filename )
          save_filename = Print::reformat_home(child->Attribute("save"));
        if ( child->Attribute("out")  )
          out_filename = child->Attribute("out");
        if ( child->Attribute("xmgrace") )
          xmg_filename = Print::reformat_home(child->Attribute("xmgrace"));
      }
 
 
      syncfilenames:
      __MPISERIALCODE( 
        // MPI code
        __TRYCODE(    Print::out.sync_filename( out_filename );
                      Print::xmg.sync_filename( xmg_filename );,
                   "Caught error while synchronizing output filenames\n" 
        )
        Print::out << "Starting genetic algorithm run on processor "
                   << ( 1 + ::LaDa::mpi::main->rank() ) << " of " 
                   << ::LaDa::mpi::main->size() << ".\n\n";,
        // Serial code
        Print::xmg.init( xmg_filename );
        Print::out.init( out_filename );
        Print::out << "Starting (serial) genetic algorithm run\n\n";
      )
 
      __ROOTCODE( (*::LaDa::mpi::main),
        Print::out << "GA Input file is located at " << evaluator_filename << "\n"
                   << "Functional Input file is located at "
                     << evaluator_filename << "\n"
                   << "Restart Input file is located at " << restart_filename << "\n"
                   << "Xmgrace output file is located at "
                     << Print::xmg.get_filename() << "\n"
                   << "Will Save to file located at "
                     << save_filename << "\n" << Print::endl
                   << Print::flush;
 
        Print::xmg << Print::Xmg :: comment << "new GA run" << Print::endl
                   << Print::Xmg::comment << "GA Input file is located at "
                                          << evaluator_filename << Print::endl
                   << Print::Xmg::comment << "Functional Input file is located at "
                                          << evaluator_filename << Print::endl
                   << Print::Xmg::comment << "Restart Input file is located at "
                                          << restart_filename << Print::endl
                   << Print::Xmg::comment << "Standard output file is located at "
                                          << Print::out.get_filename() << Print::endl
                   << Print::Xmg::comment << "Will Save to file located at "
                                          << save_filename << Print::endl;
      )
      __NOTMPIROOT( (*::LaDa::mpi::main), 
        Print::out << "Xmgrace output file is located at "
                     << Print::xmg.get_filename() << "\n"
                   << "Will Save to file located at "
                     << save_filename << "\n" << Print::endl
                   << Print::flush;
        Print::xmg << Print::Xmg :: comment << "new GA run" << Print::endl;
      )
 
#  ifdef _MPI
        // broadcasts all input files and filenames
        std::string input_str, restart_str, evaluator_str;
        LoadAllInputFiles(input_str, restart_str, evaluator_str);
        
        if ( evaluator_filename == filename ) // works for all nodes!!
          doc.Parse( evaluator_str.c_str() );
#  endif
 
      // Loads evaluator first 
      if( evaluator_filename != filename )
        __TRYASSERT( not evaluator.Load(evaluator_filename.string()),
                        "Could not load functional input from "
                     << evaluator_filename << "\n" )
      else
      {
        __TRYASSERT( not evaluator.Load(*docHandle.FirstChild("Job").Element() ),
                        "Could not load functional input from "
                     << filename << "\n" )
      }
 
      // Loads topology and assigns comms to evaluators
      __MPICODE(  doc.Parse( input_str.c_str() ); )
      parent = docHandle.FirstChild("Job").Element();
      topology.Load( *parent, evaluator );
          
      // finds <GA> ... </GA> block 
      parent = docHandle.FirstChild("Job")
                        .FirstChild("GA").Element();
      __TRYASSERT( not Load_Parameters( *parent ), 
                   "Error while reading GA attributes\n" )
 
 
 
      make_History( *parent );
      Load_Taboos( *parent );
      Load_Method( *parent );
      __DOASSERT( not  Load_Mating( *parent ),
                  "Error while loading mating operators.\n" )
      Load_CheckPoints( *parent );
      
      make_breeder();
      replacement = make_replacement();
 
      const TiXmlElement *restart_xml = parent->FirstChildElement("Restart");
      if ( restart_xml and restart_xml->Attribute("what") )
      {
        std::string str = restart_xml->Attribute("what");
        if ( str.find("all") != std::string::npos )
          do_restart |= SAVE_POPULATION | SAVE_HISTORY | SAVE_RESULTS;
        else
        {
          if ( str.find("pop") != std::string::npos )
            do_restart |= SAVE_POPULATION;
          if ( str.find("history") != std::string::npos and history )
            do_restart |= SAVE_HISTORY;
          if ( str.find("results") != std::string::npos )
            do_restart |= SAVE_RESULTS;
        }
      }
      if ( do_restart )
      {
        __NOTMPIROOT( (*::LaDa::mpi::main), doc.Parse( restart_str.c_str() ); )
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
      do_save = 0;
 
      __NOTMPIROOT( (*::LaDa::mpi::main), goto out; )
 
      do_save = SAVE_RESULTS;
      restart_xml = parent->FirstChildElement("Save");
      Print::xmg << Print::Xmg::comment << "Will Save Results";
      if ( restart_xml and restart_xml->Attribute("what") )
      {
        std::string str = restart_xml->Attribute("what");
        if ( str.find("all") != std::string::npos )
        {
          do_save |= SAVE_POPULATION | SAVE_HISTORY;
          Print::xmg << ", Population";
          if ( history )
            Print::xmg << ", and History";
          goto out;
        }
        if ( str.find("pop") != std::string::npos )
        {
          do_save |= SAVE_POPULATION;
          Print::xmg << ", Population";
        }
        if ( str.find("history") != std::string::npos and history)
        {
          do_save |= SAVE_HISTORY;
          Print::xmg << ", History";
        }
      }
    out:
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
                 << Print::Xmg::comment << (  populate_style == RANDOM_POPULATE ?
                                             "Random Starting Population":
                                             "Partitionned Starting Population" )
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
      __ROOTCODE( (*::LaDa::mpi::main), Print::xmg << Print::Xmg::comment << "Will save to "
                                             << save_filename << Print::endl; )
      std::string evalparams = evaluator.print();
      Print::xmg << Print::make_commented_string( evalparams ) << Print::endl;
      Print::xmg << Print::flush;
      return true;
    }
  } // namespace GA
} // namespace LaDa


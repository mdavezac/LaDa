//
//  Version: $Id$
//
#include<iostream>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <eo/eoDetTournamentSelect.h>
#include <eo/eoReduceMerge.h>
#include <tinyxml/tinystr.h>

#include "pescan_interface/interface.h"
#include "print/stdout.h"
#include "print/xmg.h"
#include "print/manip.h"

#include "functors.h"
#include "statistics.h"
#include "minimizergenop.h"

void keycheck( types::t_int _s=20 )
{
  std::string result; result.resize(_s);
  for (types::t_int i=0 ; i < _s; ++i)
    result[i] = rng.flip() ? '1' : '0';
  Print::out << "KeyCheck: " << result << Print::endl;
}


namespace GA
{
# define OPENXMLINPUT \
    TiXmlDocument doc( filename.c_str() ); \
    TiXmlHandle docHandle( &doc ); \
    if  ( !doc.LoadFile() ) \
    { \
      std::cerr << doc.ErrorDesc() << std::endl; \
      throw "Could not load input file in  Darwin<T_EVALUATOR> "; \
    } 

  template<class T_EVALUATOR>
  Darwin<T_EVALUATOR> :: ~Darwin()
  {
    if ( store ) delete store;
    if ( objective ) delete objective;
    if ( evaluation ) delete evaluation;
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
      else if ( str.compare("seed")==0 ) // seed from given number
      { 
        types::t_int d = att->IntValue();
        Print::xmg << Print::Xmg::comment <<  "Seed: " 
                   << Print::fixed << Print::setw(8) << d << Print::endl;
#ifdef _MPI
        rng.reseed( std::abs(d) * (mpi::main.rank()+1) );
#else
        rng.reseed( std::abs(d) );
#endif
      }
      else if ( str.compare("populate")==0 ) // seed from given number
      {
        std::string string = att->Value();
        Print::xmg << Print::Xmg::comment << "Populate Style: ";
        if ( string.compare("partition") == 0 )
        {
          populate_style = PARTITION_POPULATE;
          Print::xmg << "Partition Populate";
        }
        else
          Print::xmg << "Random Populate";
        Print::xmg << Print::endl;
      }
      else
        evaluator.LoadAttribute( *att );
    }

    // some checking
    if ( std::floor( pop_size * replacement_rate ) == 0 )
    {
      Print::xmg << Print::Xmg::comment << "Error: replacement_rate is too small." << Print::endl 
                 << Print::Xmg::comment << "Error: setting replacement_rate too 1.0 / pop_size ." << Print::endl;
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
    // checks if there are more than one taboo list
    const TiXmlElement *child = _parent.FirstChildElement("History");
    if ( not child )
      return;
    Print::xmg << Print::Xmg::comment << "Track History" << Print::endl;
    history = new History< t_Individual>;
    eostates.storeFunctor(history);
  }
  
  // creates history object if required
  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: Load_Method(const TiXmlElement &_parent)
  {
    // checks if there are more than one taboo list
    const TiXmlElement *child = _parent.FirstChildElement("Objective");
    if ( not child ) child = _parent.FirstChildElement("Method");
    objective = t_ObjectiveType :: new_from_xml( *child );
    if( not objective )
      throw std::runtime_error(" Could not find Objective tag... ");

    const TiXmlElement *store_xml = _parent.FirstChildElement("Store");
    if ( store_xml )
        store = new typename t_Store::FromObjective( evaluator, *store_xml );
    else
      store = new typename t_Store::FromObjective( evaluator, objective, *child );
    if ( not store )
      store = new typename t_Store::Optima( evaluator, *child );
    if ( not store )
      throw std::runtime_error( "Memory Error: could not create Store object \n");
    Print::xmg << Print::Xmg::comment << "Store: " << store->what_is() << Print::endl;

    if( history )
      evaluation = new Evaluation::WithHistory<t_GATraits>
                                              ( evaluator, *objective, *store, history );
    if ( not evaluation )
      evaluation = new Evaluation::Base<t_GATraits>( evaluator, *objective, *store );

    ranking = Ranking::new_from_xml<t_GATraits>( _parent );
    if( ranking ) Print::xmg << Print::Xmg::comment << ranking->what_is() << Print::endl;
      
  }
  
  // create Taboos
  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: Load_Taboos( const TiXmlElement &_node )
  {
    // checks if there are more than one taboo list
    if ( not _node.FirstChildElement("Taboos") ) return;
    const TiXmlElement &parent = *_node.FirstChildElement("Taboos");
    const TiXmlElement *child = parent.FirstChildElement();

    // creates Taboo container if there are more than one taboo list
    taboos = new Taboos<t_Individual>;

    // creates pop taboo
    child = parent.FirstChildElement("Population");
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
    if (child)
    {
      Print::xmg << Print::Xmg::comment << "Offspring Taboo" << Print::endl;
      OffspringTaboo<t_GATraits> *offspringtaboo 
         = new OffspringTaboo<t_GATraits>( &offsprings );
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
    if ( not breeder_ops )
    {
      std::cerr << "Error while creating breeding operators" << std::endl;
      throw "Error while creating operators in  Darwin<T_EVALUATOR>  :: make_GenOp ";
    }
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
    GenCount &generation_counter = continuator->get_generation_counter();

    // The following synchronizes Results::nb_val 
    // all procs. Should always be there!!
#ifdef _MPI
        Synchronize<types::t_unsigned> *synchro = new Synchronize<types::t_unsigned>( evaluation->nb_eval );
        eostates.storeFunctor( synchro );
        continuator->add( *synchro );

        synchro = new Synchronize<types::t_unsigned>( evaluation->nb_grad );
        eostates.storeFunctor( synchro );
        continuator->add( *synchro );
#endif

    // Creates SaveEvery object if necessary
    if (      _parent.FirstChildElement("Save") 
         and  _parent.FirstChildElement("Save")->Attribute("every") )
    {
      types::t_int n;
      _parent.FirstChildElement("Save")->Attribute("every", &n);
      
      if ( n >= 0 and do_save )
      {
        Print::xmg << Print::Xmg::comment << "Will Save Every " << n << " Generations " << Print::endl;
        SaveEvery<t_Darwin> *save = new SaveEvery<t_Darwin>( *this, &Darwin::Save, std::abs(n) );
        eostates.storeFunctor( save );
        continuator->add( *save );
      }
    }
 
    // Creates Statistics -- only one type available...
    const TiXmlElement *child = _parent.FirstChildElement("Statistics");
    if( child )
    {
      continuator->add( eostates.storeFunctor( new TrueCensus< t_GATraits >() ) );
      Print::xmg << Print::Xmg::comment << "Statistics: True population size, discounting twins" << Print::endl;
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
      
      terminator = new Terminator< types::t_unsigned, std::less<types::t_unsigned>, t_GATraits >
                                 ( evaluation->nb_eval, (types::t_unsigned) abs(max),
                                   std::less<types::t_unsigned>(), "nb_eval < term" );
      eostates.storeFunctor( terminator );
      continuator->add( *terminator );
      Print::xmg << Print::Xmg::comment << "Terminating after " << max << " evaluations" << Print::endl;
      
      // end if max
    }
    
    // Load object specific continuator
    eoF<bool> *specific = evaluator.LoadContinue( _parent );
    if ( specific )
    {
      eostates.storeFunctor(specific); 
      continuator->add( eostates.storeFunctor(new Continuator<t_GATraits>(*specific)) );
    }

    // Creates Print object
    {
      typedef PrintGA< Store::Base<t_GATraits>, Evaluation::Base<t_GATraits> > t_PrintGA;
      t_PrintGA* printga = new t_PrintGA( *store, *evaluation, generation_counter, do_print_each_call);
      eostates.storeFunctor(printga);
      continuator->add( *printga );
    }

    // Print Offsprings
    child = _parent.FirstChildElement("PrintOffsprings");
    if ( child )
    {
      PrintFitness<t_GATraits> *printfitness = new PrintFitness<t_GATraits> ( generation_counter );
      continuator->add( *printfitness );
      eostates.storeFunctor( printfitness );
    }

  } // end of make_check_point
  
  template<class T_EVALUATOR>
  bool Darwin<T_EVALUATOR> :: Restart()
  {
    TiXmlDocument doc( restart_filename.c_str() );
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
                     << " individuals into island " << islands.size() << "" << Print::endl;
        }
      } // end of "Population" Restart
    }
    return true;
  }

  template<class T_EVALUATOR>
  bool Darwin<T_EVALUATOR> :: Save()
  {
    if ( not do_save )
      return true;

    TiXmlDocument doc;
    TiXmlElement *node;
    if ( save_filename == filename )
    {
      doc.LoadFile(save_filename.c_str());
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

    if ( not doc.SaveFile(save_filename.c_str() ) )
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
    eoTruncate<t_Individual>* truncate       = new  eoTruncate<t_Individual>;
    eoMerge<t_Individual>* merge             = new  eoPlus<t_Individual>;
    eoReduceMerge<t_Individual>* reducemerge = new eoReduceMerge<t_Individual>( *truncate, *merge );
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
    if (not sibling)
      return NULL;
 
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
        eoGenOp<t_Individual> *taboo_op = make_genetic_op( *sibling->FirstChildElement(), NULL);
        if ( not taboo_op )
          Print::xmg << Print::Xmg::removelast;
        else
        {
          Print::xmg << Print::Xmg::comment << "TabooOp End" << Print::endl;
          eoMonOp<t_Individual> *utterrandom;
          utterrandom = new mem_monop_t<t_GATraits>
                              ( evaluator, &t_Evaluator::initialize, std::string( "Initialize" ) );
          eostates.storeFunctor(utterrandom);
          this_op = new TabooOp<t_Individual> ( *taboo_op, *taboos, 
                                                (types::t_unsigned)(pop_size+1),
                                                *utterrandom );
          eostates.storeFunctor( static_cast< TabooOp<t_Individual> *>(this_op) );
          is_gen_op = true;
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
                         ( evaluator, &t_Evaluator::initialize, std::string( "UtterRandom" ) );
        eostates.storeFunctor( static_cast< TabooOp<t_Individual> *>(this_op) );
      }
      else if ( str.compare("Minimizer") == 0 )
      {
        Minimizer_Functional<t_GATraits> *func = 
          new Minimizer_Functional<t_GATraits>( *evaluation, *taboos ); 
        if ( not func )
          throw std::runtime_error( "Memory Allocation Error");
        MinimizerGenOp<t_GATraits> *mingenop
            = new MinimizerGenOp<t_GATraits>( *func );
        if ( not mingenop )
          { delete func; this_op = NULL; }
        else if ( not mingenop->Load( *sibling ) )
          { delete mingenop; delete func; this_op = NULL; }
        else
        {
          eostates.storeFunctor( mingenop );
          this_op = mingenop;
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
              Print::out << " Failure, or No Operators Found in \"And\" operator \n";
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
          Print::xmg << Print::Xmg::addtolast << " period = " << period << Print::endl;
          this_op = new PeriodicOp<t_GATraits>( *this_op, (types::t_unsigned) abs(period),
                                              continuator->get_generation_counter(), eostates );
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
          static_cast< SequentialOp<t_GATraits>* >(current_op)->add( *this_op,
                                                                   static_cast<double>(prob) );
        else if ( current_op->className().compare("GA::ProportionalOp") == 0 )
          static_cast< ProportionalOp<t_GATraits>* >(current_op)->add( *this_op, 
                                                                     static_cast<double>(prob) );
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

    select = new eoDetTournamentSelect<t_Individual>(tournament_size);
    breeder = new Breeder<t_GATraits>(*select, *breeder_ops, continuator->get_generation_counter() );
    breeder->set_howmany(replacement_rate);

    eostates.storeFunctor(breeder);
    eostates.storeFunctor(select);
  }
  
  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: populate ()
  {
    islands.resize( nb_islands );
    typename t_Islands :: iterator i_pop = islands.begin();
    typename t_Islands :: iterator i_end = islands.end();
    types::t_unsigned target = pop_size;
#ifdef _MPI
    types::t_int residual = target % (mpi::main.size());
    target = target / mpi::main.size();
    if ( mpi::main.rank() < residual )
      ++target;
#endif
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
  void Darwin<T_EVALUATOR> :: random_populate ( t_Population &_pop, types::t_unsigned _size)
  {
    while ( _pop.size() < _size )
    {
      t_Individual indiv;
      evaluator.initialize(indiv);
      if ( not ( taboos and (*taboos)(indiv) ) )
        _pop.push_back(indiv);
    } // while ( i_pop->size() < target )
    _pop.resize( _size );
  }
  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: partition_populate ( t_Population &_pop, types::t_unsigned _size)
  {
    while ( _pop.size() < _size )
    {
      // first random object
      t_Individual indiv;
      evaluator.initialize(indiv);
      types::t_unsigned start = 0;  
      types::t_unsigned end = indiv.Object().size();  
    
      // indiv: ----
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;
      
      indiv.Object().mask( start, end );         // indiv: ****
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;
      
      indiv.Object().mask( start, end/2 );       // indiv: --**
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end );         // indiv: **--
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end/4 );       // indiv: -*--
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end );         // indiv: *-**
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( end/2, end*3/4 );     // indiv: *--*
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end );         // indiv: -**-
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( end/4, end/2 );       // indiv: --*-
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end );         // indiv: **-*
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end/2 );       // indiv: ---*
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
      if ( _pop.size() >= _size ) break;

      indiv.Object().mask( start, end );         // indiv: *---
      if ( taboos and (*taboos)(indiv) ) continue;
      _pop.push_back(indiv);
    }
    _pop.resize( _size );
  }


  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: run()
  {
#ifdef _MPI 
    for( types::t_int i = 0; i < mpi::main.size(); ++i )
    {
      if ( mpi::main.rank() == i )
      {
        struct timeval tv;
        struct timezone tz;
        gettimeofday(&tv,&tz);
        rng.reseed(tv.tv_usec);
      }
      mpi::main.barrier();
    }
#endif
    Print::xmg << Print::flush;
    populate();
    offsprings.clear();
    typename t_Islands :: iterator i_island_begin = islands.begin();
    typename t_Islands :: iterator i_island_end = islands.end();
    typename t_Islands :: iterator i_island;
    for ( i_island = i_island_begin; i_island != i_island_end; ++i_island )
    {
      (*evaluation)(offsprings, *i_island); // A first eval of pop.
#ifdef _MPI // "All Gather" new population
      breeder->synchronize_offsprings( *i_island );
#endif 
    }
    types::t_unsigned n = 0;

    do
    {
      try
      {
        Print::out << "Iteration " << n << Print::endl;
        Print::xmg << Print::flush;
        ++n;
        i_island = i_island_begin;
        for (int i=0; i_island != i_island_end; ++i, ++i_island )
        {

          types::t_unsigned pSize = i_island->size();
          offsprings.clear(); // new offsprings
          
          (*breeder)(*i_island, offsprings);
          
          (*evaluation)(*i_island, offsprings); // eval of parents + offsprings if necessary

#ifdef _MPI // "All Gather" offsprings -- note that replacement scheme is deterministic!!
          breeder->synchronize_offsprings( offsprings );
          if(history) history->synchronize();
#endif 
         
          (*replacement)(*i_island, offsprings); // after replace, the new pop. is in population
          
          if (pSize > i_island->size())
              throw std::runtime_error("Population shrinking!");
          else if (pSize < i_island->size())
              throw std::runtime_error("Population growing!");
        }
      }
      catch (std::exception& e)
      {
            std::string s = e.what();
            s.append( " in eoEasyEA");
            throw std::runtime_error( s );
      }
    } while ( continuator->apply( i_island_begin, i_island_end ) );

    Print::xmg << Print::flush;
    Save();
  }


#ifdef _MPI
  template<class T_EVALUATOR>
  void Darwin<T_EVALUATOR> :: LoadAllInputFiles(std::string &_input, 
                                                std::string &_restart, 
                                                std::string &_evaluator ) 
  {
    if ( mpi::main.is_root_node() )
    { 
      OPENXMLINPUT
      TIXML_OSTREAM stream;
      const TiXmlElement *parent = doc.RootElement();
      stream << *parent;
      _input = stream.c_str();
      _restart = _input; 
      _evaluator = _input;

      if ( restart_filename == filename )
      {
        doc.LoadFile( restart_filename.c_str() );
        if  ( !doc.LoadFile() ) 
        { 
          std::cerr << doc.ErrorDesc() << std::endl; 
          std::cerr << "Could not load restart file " << std::endl;
          _restart = "";
          goto nextfilename;
        } 

        parent = doc.RootElement();
        TIXML_OSTREAM stream;
        _restart = stream.c_str();
      }
nextfilename:
      if ( evaluator_filename != filename )
      {
        doc.LoadFile( evaluator_filename.c_str() );
        if  ( !doc.LoadFile() ) 
        { 
          std::cerr << doc.ErrorDesc() << std::endl; 
          throw std::runtime_error("Could not load restart file"); 
        } 

        parent = doc.RootElement();
        TIXML_OSTREAM stream;
        _evaluator = stream.c_str();
      }
    }

    mpi::BroadCast bc( mpi::main );
    bc << _input << _restart << _evaluator
       << mpi::BroadCast::allocate
       << _input << _restart << _evaluator
       << mpi::BroadCast::broadcast
       << _input << _restart << _evaluator;
  }
#endif

  template<class T_EVALUATOR>
  bool Darwin<T_EVALUATOR> :: Load(const std::string &_filename) 
  {
    filename = _filename;
    evaluator_filename = _filename;
    restart_filename = _filename;
    save_filename = _filename;

    // first loads all inputfiles 
#ifdef _MPI
    if ( mpi::main.is_root_node() )
    {
#endif
      OPENXMLINPUT
      const TiXmlElement *parent = docHandle.FirstChild("Job")
                                            .FirstChild("GA").Element();
      if ( not parent )
        return false;
      const TiXmlElement *child = parent->FirstChildElement("Filenames");
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
          Print::out.init( child->Attribute("out") );
        if ( child->Attribute("xmgrace") )
        {
          std::string f = Print::reformat_home(child->Attribute("xmgrace"));
          Print::xmg.init( f );
        }
      }

      Print::xmg << Print::Xmg :: comment << "new GA run" << Print::endl;
#ifdef _MPI
    }
    Print::out.sync_filename();
#endif

    Print::out << "Starting new genetic algorithm runs\n\n"
               << "GA Input file is located at " << filename << "\n"
               << "Functional Input file is located at " << evaluator_filename << "\n"
               << "Restart Input file is located at " << restart_filename << "\n"
               << "Xmgrace output file is located at " << Print::xmg.get_filename() << "\n"
               << "Will Save to file located at " << save_filename << "\n\n";

#ifdef _MPI
    // broadcasts all input files and filenames
    std::string input_str, restart_str, evaluator_str;
    LoadAllInputFiles(input_str, restart_str, evaluator_str);

    TiXmlDocument doc;
    if ( evaluator_filename == filename ) // works for all nodes!!
      doc.Parse( evaluator_str.c_str() );
    TiXmlHandle docHandle(&doc);
#endif
      
    // Loads evaluator first 
    if ( evaluator_filename != filename )
    { 
      if ( not evaluator.Load(evaluator_filename) )
        return false;
    } // for following, we know docHandle.FirstChild("Job") exists
    else if ( not evaluator.Load(*docHandle.FirstChild("Job").Element() ) )
      return false;
          
    // finds <GA> ... </GA> block 
#ifdef _MPI
    doc.Parse( input_str.c_str() );
    const TiXmlElement *parent = docHandle.FirstChild("Job")
                                          .FirstChild("GA").Element();
#endif
    if ( not Load_Parameters( *parent ) )
      return false;



    make_History( *parent );
    Load_Taboos( *parent );
    Load_Method( *parent );
    if ( not  Load_Mating( *parent ) )
      return false;
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
#ifdef _MPI
      if ( not mpi::main.is_root_node() )
        doc.Parse( restart_str.c_str() );
#endif
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
#ifdef _MPI
    if ( not mpi::main.is_root_node() )
    {
      do_save = 0;
      return true;
    }
#endif 
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
out:  Print::xmg << Print::endl;

    {
      Print::xmg << Print::Xmg::comment << "Mating Tournament Size: " << tournament_size << Print::endl
                 << Print::Xmg::comment << "Offspring Replacement Rate: " << replacement_rate << Print::endl
                 << Print::Xmg::comment << "Population Size: " << pop_size << Print::endl;
      if ( max_generations )
        Print::xmg << Print::Xmg::comment << "Maximum Number of Generations: " << max_generations
                   << Print::endl;
      else
        Print::xmg << Print::Xmg::comment << "Unlimited Number of Generations" << Print::endl;
      Print::xmg << Print::Xmg::comment << "Number of Islands: " << nb_islands << Print::endl 
                 << Print::Xmg::comment << "Will save to " << save_filename << Print::endl;
    }
    return true;
  }


#ifdef _MPI
  template<class T_EVALUATOR>
  bool Darwin<T_EVALUATOR> :: broadcast_islands( mpi::BroadCast &_bc )
  {
    types::t_int n = islands.size();
    if ( not _bc.serialize( n ) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      islands.resize(n);
    typename t_Islands :: iterator i_pop = islands.begin();
    typename t_Islands :: iterator i_pop_end = islands.end();
    for( ; i_pop != i_pop_end; ++i_pop )
    {
      n = i_pop->size();
      if ( not _bc.serialize( n ) ) return false;
      if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
        i_pop->resize(n);
      typename t_Population :: iterator i_indiv = i_pop->begin();
      typename t_Population :: iterator i_indiv_end = i_pop->end();
      for( ; i_indiv != i_indiv_end; ++i_indiv )
        if( not i_indiv->broadcast( _bc ) ) return false;
    }
    return true;
  }
#endif

} // namespace GA


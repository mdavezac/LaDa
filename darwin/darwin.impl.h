#include<iostream>

#include <eo/eoDetTournamentSelect.h>
#include <eo/eoReduceMerge.h>
#include <tinyxml/tinystr.h>


#include "print_xmgrace.h"
#include "functors.h"
#include "statistics.h"

namespace darwin
{
#ifdef _MPI
# define OPENXMLINPUT \
    TiXmlDocument doc( filename.c_str() ); \
    TiXmlHandle docHandle( &doc ); \
    if  ( mpi::main.rank() == mpi::ROOT_NODE ) \
    if  ( !doc.LoadFile() ) \
    { \
      std::cout << doc.ErrorDesc() << std::endl; \
      throw "Could not load input file in  Darwin<T_INDIVIDUAL,T_EVALUATOR> "; \
    } 
#else
# define OPENXMLINPUT \
    TiXmlDocument doc( filename.c_str() ); \
    TiXmlHandle docHandle( &doc ); \
    if  ( !doc.LoadFile() ) \
    { \
      std::cout << doc.ErrorDesc() << std::endl; \
      throw "Could not load input file in  Darwin<T_INDIVIDUAL,T_EVALUATOR> "; \
    } 
#endif

  template<class T_INDIVIDUAL, class T_EVALUATOR>
  Darwin<T_INDIVIDUAL,T_EVALUATOR> :: ~Darwin()
  {
    if ( results ) delete results; results = NULL; 
    typename t_minimizer_list :: iterator i_it = minimizers.begin();
    typename t_minimizer_list :: iterator i_end = minimizers.end();
    for( ; i_it != i_end; ++i_it )
      delete *i_it;
  }
  // reads in different parameters
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_Parameters(const TiXmlElement &_parent)
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
        std::ostringstream sstr; 
        sstr << "Seed: " << d;
        std::string str = sstr.str();
        printxmg.add_comment( str );
        rng.reseed( std::abs(d) * (mpi::main.rank()+1) );
      }
      else if ( str.compare("populate")==0 ) // seed from given number
      {
        std::string string = att->Value();
        std::ostringstream sstr; 
        sstr << "Populate Style: ";
        if ( string.compare("partition") == 0 )
        {
          populate_style = PARTITION_POPULATE;
          sstr << "Partition Populate";
        }
        else
          sstr << "Random Populate";
        printxmg.add_comment( sstr.str() );
      }
      else
        evaluator.LoadAttribute( *att );
    }

    // some checking
    if ( std::floor( pop_size * replacement_rate ) == 0 )
    {
      printxmg.add_comment("Error: replacement_rate is too small.");
      printxmg.add_comment("Error: setting replacement_rate too 1.0 / pop_size .");
      replacement_rate = 1.00 / pop_size + 10*types::tolerance;
    }

    return true;
  }

  // creates history object if required
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_History(const TiXmlElement &_parent)
  {
    // checks if there are more than one taboo list
    const TiXmlElement *child = _parent.FirstChildElement("History");
    if ( not child )
      return;
    printxmg.add_comment( "Track History" );
    history = new History< t_Individual, std::list<t_Individual> >;
    eostates.storeFunctor(history);
  }
  
  // creates history object if required
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_Method(const TiXmlElement &_parent)
  {
    // checks if there are more than one taboo list
    const TiXmlElement *child = _parent.FirstChildElement("Method");
    if( child and child->Attribute("type") ) 
    {
      std::string str = child->Attribute("type");
      if (    str.compare("Best Of")==0 
           or str.compare("BO")==0 )
        results = new BestOf<t_Individual, t_Evaluator, t_Population>( *child, evaluator, taboos );
      if ( str.compare("Target")==0 )
        results = new Target<t_Individual, t_Evaluator, t_Population>( *child, evaluator, taboos );
      if (    str.compare("ConvexHull")==0 
           or str.compare("CH")==0 )
        results = new ConvexHull<t_Individual, t_Evaluator, t_Population>( *child, evaluator, taboos );
    }
    if( not results )
      results = new Optimum<t_Individual, t_Evaluator, t_Population>( *child, evaluator, taboos );
  }
  
  // create Taboos
  template<class T_INDIVIDUAL, class T_EVALUATOR >
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_Taboos( const TiXmlElement &_node )
  {
    // checks if there are more than one taboo list
    if ( not _node.FirstChildElement("Taboos") )
      return;
    const TiXmlElement &parent = *_node.FirstChildElement("Taboos");
    const TiXmlElement *child = parent.FirstChildElement();

    // creates Taboo container if there are more than one taboo list
    types::t_unsigned nb_taboos = 0;
    for( ; child and nb_taboos < 2; child = child->NextSiblingElement() )
    {
      std::string name = child->Value();
      if (    name.compare("PopTaboo") == 0
           or name.compare("AgeTaboo") == 0
           or name.compare("OffspringTaboo") == 0 
           or ( name.compare("HistoryTaboo") == 0 and history) )
        ++nb_taboos;
      else
      {
        eoMonOp<const t_Object> *_op = evaluator.LoadTaboo( *child );
        if ( not _op )
          continue;
        delete _op;
        ++nb_taboos;
      }
    }
    if ( nb_taboos < 1 )
      return; // no Taboo tags in Taboos tag
    if ( nb_taboos > 1 ) // creates a container
    {
      taboos = new Taboos<t_Individual>;
      if ( not taboos )
      {
        std::cerr << "Could not allocate memory for taboos!!" << std::endl;
        throw "";
      }
      eostates.storeFunctor(taboos);
    }

    IslandsTaboos<t_Individual> *poptaboo = NULL;
    Taboo<t_Individual> *offspringtaboo = NULL;
    agetaboo = NULL;

    // creates age taboo -- updater is created in make_checkpoint
    child = parent.FirstChildElement("AgeTaboo");
    if (child)
    {
      printxmg.add_comment("Age Taboo");
      agetaboo = new Taboo< t_Individual, std::list<t_Individual> >;
      eostates.storeFunctor(agetaboo);
      if ( not taboos )
      {
        taboos = agetaboo;
        return;
      }
      static_cast< Taboos<t_Individual>* >(taboos)->add( agetaboo );
    }
    
    // creates pop taboo
    child = parent.FirstChildElement("PopTaboo");
    if (child)
    {
      printxmg.add_comment("Population Taboo");
      poptaboo = new IslandsTaboos<t_Individual>( islands );
      eostates.storeFunctor(poptaboo);
      if ( not taboos ) 
      {
        taboos = poptaboo;
        return;
      }
      static_cast< Taboos<t_Individual>* >(taboos)->add( poptaboo );
    }
    
    // creates offspring taboo
    child = parent.FirstChildElement("OffspringTaboo");
    if (child)
    {
      printxmg.add_comment("Offspring Taboo");
      offspringtaboo = new OffspringTaboo<t_Individual>( &offsprings );
      eostates.storeFunctor(offspringtaboo);
      if ( not taboos )
      {
        taboos = offspringtaboo;
        return;
      }
      static_cast< Taboos<t_Individual>* >(taboos)->add( offspringtaboo );
    }

    // makes history a taboo list if it exists
    child = parent.FirstChildElement("History Taboo");
    if (child and history)
    {
      printxmg.add_comment("History Taboo");
      if ( not taboos )
      {
        taboos = history;
        return;
      }
      static_cast< Taboos<t_Individual>* >(taboos)->add( history );
    }
    else if (child)
      std::cerr << "HistoryTaboo found in Taboos tags, but not History tag found!!" << std::endl
                << "Include History tag if you want HistoryTaboo" << std::endl;


    child = parent.FirstChildElement();
    for( ; child ; child = child->NextSiblingElement() )
    {
      eoMonOp<const t_Object> *op = evaluator.LoadTaboo( *child );
      if ( not op )
        continue;
      eostates.storeFunctor( op );
      TabooFunction<t_Individual> *func = new TabooFunction<t_Individual>( *op );
      eostates.storeFunctor(func);
      if ( not taboos )
      {
        taboos = func;
        return;
      }
      static_cast< Taboos<t_Individual>* >(taboos)->add( func );
    }

      
  }

  // Loads mating operations
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_Mating (const TiXmlElement &_parent)
  {
    const TiXmlElement *child = _parent.FirstChildElement("Breeding");
    if ( not child )
    {
      std::cerr << "No Breeding operator found in input. " << std::endl;
      std::cerr << "Giving Up" << std::endl;
      return false;
    }
    printxmg.add_comment("Breeding Operator Begin");
    breeder_ops = new SequentialOp<t_Individual>;
    eostates.storeFunctor( breeder_ops );
    breeder_ops = make_genetic_op(*child->FirstChildElement(), breeder_ops);
    printxmg.add_comment("Breeding Operator End");
    if ( not breeder_ops )
    {
      std::cout << "Error while creating breeding operators" << std::endl;
      throw "Error while creating operators in  Darwin<T_INDIVIDUAL,T_EVALUATOR>  :: make_GenOp ";
    }
    return true;
  }
  
  // Loads CheckPoints
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load_CheckPoints (const TiXmlElement &_parent)
  {
    // contains all checkpoints
    std::string str = "stop"; 
    if (      _parent.FirstChildElement("Filenames") 
         and  _parent.FirstChildElement("Filenames")->Attribute("stop") )
      str = _parent.FirstChildElement("Filenames")->Attribute("stop");

    continuator = new IslandsContinuator<t_Individual>(max_generations, str );
    eostates.storeFunctor( continuator );
    GenCount &generation_counter = continuator->get_generation_counter();

    // Creates SaveEvery object if necessary
    if (      _parent.FirstChildElement("Save") 
         and  _parent.FirstChildElement("Save")->Attribute("every") )
    {
      types::t_int n;
      _parent.FirstChildElement("Save")->Attribute("every", &n);
      
      if ( n >= 0 and do_save )
      {
        std::ostringstream sstr;
        sstr << "Will Save Every " << n << " Generations ";
        printxmg.add_comment( sstr.str() );
        SaveEvery<t_Darwin> *save = new SaveEvery<t_Darwin>( *this, &Darwin::Save, std::abs(n) );
        eostates.storeFunctor( save );
        continuator->add( *save );
      }
    }
 
    // Creates Print object
    {
      typedef Print< Results<t_Individual, t_Evaluator, t_Population> > t_Print;
      t_Print* print = new t_Print( *results, generation_counter, do_print_each_call);
      eostates.storeFunctor(print);
      continuator->add( *print );
    }
    
    // Print Offsprings
    const TiXmlElement *child = _parent.FirstChildElement("PrintOffsprings");
    if ( child )
    {
      PrintFitness<t_Individual> *printfitness = new PrintFitness<t_Individual> ( generation_counter );
      continuator->add( *printfitness );
      eostates.storeFunctor( printfitness );
    }

    // Creates Statistics -- only one type available...
    child = _parent.FirstChildElement("Statistics");
    if( child )
    {
      continuator->add( eostates.storeFunctor( new TrueCensus< t_Individual, t_Population >() ) );
      printxmg.add_comment("Statistics: True population size, discounting twins");
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

      if ( ref.compare("evaluation") == 0 )
      {
        terminator = new Terminator< types::t_unsigned, std::less<types::t_unsigned>, t_Individual >
                                   ( results->nb_eval, (types::t_unsigned) abs(max),
                                     std::less<types::t_unsigned>(), "nb_eval < term" );
        eostates.storeFunctor( terminator );
        continuator->add( *terminator );
        std::ostringstream sstr;
        sstr << "Terminating after " << max << " evaluations";
        printxmg.add_comment(sstr.str());
      }
      else if ( ref.compare("gradient") == 0 )
      {
        terminator = new Terminator< types::t_unsigned, std::less<types::t_unsigned>, t_Individual >
                                   ( results->nb_grad, (types::t_unsigned) abs(max),
                                     std::less<types::t_unsigned>(), "nb_grad < term" );
        eostates.storeFunctor( terminator );
        continuator->add( *terminator );
        std::ostringstream sstr;
        sstr << "Terminating after " << max << " gradient calls";
        printxmg.add_comment(sstr.str());
      }
      // end if max
    }
    
    // Creates Age taboo Updater
    child = _parent.FirstChildElement("Taboos");
    if( child ) child = child->FirstChildElement("AgeTaboo");
    if ( child and agetaboo )
    {
      int d = 0;
      child->Attribute("lifespan", &d );
      types::t_unsigned length = ( d >=0 ) ? (types::t_int) abs(d) : UINT_MAX;
      bool print_out = false;
      if ( child->Attribute("printout") )
      {
        std::string str = child->Attribute("printout");
        if ( str.compare("true") == 0 )
          print_out = true;
      }
      UpdateAgeTaboo< t_Individual > *updateagetaboo
             = new UpdateAgeTaboo<t_Individual > ( *agetaboo, generation_counter,
                                                   length, print_out);
      std::ostringstream sstr;
      sstr << "# Age Taboo, lifespan=" << d << std::endl;
      std::string str = sstr.str();
      printxmg.add_comment(str);
      eostates.storeFunctor(updateagetaboo);
      continuator->add(*updateagetaboo);
    }


    // Nuclear Winter -- only for agetaboo
    child = _parent.FirstChildElement("Taboos");
    if ( child ) child = child->FirstChildElement("NuclearWinter");
    if ( child and agetaboo)
    {
      eoGenOp<t_Individual> *nuclear_op;
    
      // first creates the nuclear op from input
      printxmg.add_comment("Nuclear Operator Begin");
      nuclear_op = make_genetic_op( *child->FirstChildElement() );
      printxmg.add_comment("Nuclear Operator End");
      if ( not nuclear_op )
        throw "Error while creating operators in  Darwin<t_Individual>  :: make_GenOp ";
      
      // creates the NuclearWinter 
      nuclearwinter = new NuclearWinter<t_Individual, t_Population >
                                       ( *taboos, *breeder_ops, *nuclear_op,
                                         replacement_rate );
      printxmg.add_comment("NuclearWinter");
      eostates.storeFunctor( nuclearwinter );
      continuator->add(*nuclearwinter);
    }


    // Load object specific continuator
    eoF<bool> *specific = evaluator.LoadContinue( _parent );
    if ( not specific )
      return;
    eostates.storeFunctor(specific); 
    continuator->add( eostates.storeFunctor(new Continuator<t_Individual>(*specific)) );

  } // end of make_check_point
  
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Restart()
  {
    TiXmlDocument doc( restart_filename.c_str() );
    TiXmlHandle docHandle( &doc ); 
    if  ( !doc.LoadFile() ) 
    { 
      std::cout << doc.ErrorDesc() << std::endl; 
      return false;
    } 
    TiXmlElement *parent = docHandle.FirstChild("Restart").Element();
    if ( not parent ) return false;

    return Restart( *parent );
  }

  // Loads restarts operations
  // Should be in a Restart XML Section
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Restart (const TiXmlElement &_node)
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
        results->Restart(*child);
      else if (     name.compare("History") == 0
                and ( do_restart & SAVE_HISTORY ) and history )
      {
        LoadObject<t_Object, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
        history->clear();
        if ( LoadIndividuals( *child, loadop, *history) )
        {
          std::ostringstream sstr;
          sstr << "Reloaded " << history->size() 
               << " individuals into history ";
          printxmg.add_comment( sstr.str() );
        }
        
      } // end of "Population" Restart
      else if (     name.compare("Population") == 0
                and ( do_restart & SAVE_POPULATION ) )
      {
        LoadObject<t_Object, t_Evaluator> loadop( evaluator, &t_Evaluator::Load, LOADSAVE_SHORT);
        islands.clear();
        const TiXmlElement *island_xml = child->FirstChildElement("Island");
        for(; island_xml; island_xml = island_xml->NextSiblingElement("Island") )
        {
          if( islands.size() > nb_islands )
            break;
          t_Population pop;
          LoadIndividuals( *island_xml, loadop, pop, pop_size);
          std::ostringstream sstr;
          islands.push_back( pop );
          sstr << "Loaded " << pop.size() 
               << " individuals into island " << islands.size();
          printxmg.add_comment( sstr.str() );
        }
      } // end of "Population" Restart
    }
    return true;
  }

  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Save()
  {
    if ( not do_save )
      return true;

    TiXmlDocument doc;
    TiXmlElement *node;
    if ( restart_filename == filename )
    {
      doc.LoadFile(restart_filename.c_str());
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

    if ( do_save & SAVE_RESULTS )
      results->Save( *node );
    if ( do_save & SAVE_HISTORY and history)
    {
      TiXmlElement *xmlhistory = new TiXmlElement("History");
      SaveObject<t_Object, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
      SaveIndividuals( *xmlhistory, saveop, history->begin(), history->end());
      node->LinkEndChild( xmlhistory );
    }
    if ( do_save & SAVE_POPULATION )
    {
      SaveObject<t_Object, t_Evaluator> saveop( evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
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

    doc.SaveFile(restart_filename.c_str());
    return true;
  }
  
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  eoReplacement<T_INDIVIDUAL>* Darwin<T_INDIVIDUAL,T_EVALUATOR> :: make_replacement()
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
  template<class T_INDIVIDUAL, class T_EVALUATOR >
  eoGenOp<T_INDIVIDUAL>* Darwin<T_INDIVIDUAL,T_EVALUATOR> :: make_genetic_op( const TiXmlElement &el,
                                                                  eoGenOp<t_Individual> *current_op )
  {
    eoOp<t_Individual>* this_op;
    const TiXmlElement *sibling = &el;
    if (not sibling)
      return NULL;
 
    printxmg.indent();

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
        printxmg.add_comment("TabooOp Begin");
        eoGenOp<t_Individual> *taboo_op = make_genetic_op( *sibling->FirstChildElement(), NULL);
        if ( not taboo_op )
          printxmg.remove_last();
        else
        {
          printxmg.add_comment("TabooOp End");
          eoMonOp<t_Individual> *utterrandom;
          utterrandom = new mem_monop_indiv_t<t_Evaluator, t_Individual>
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
        printxmg.deindent();
        this_op = make_genetic_op( *sibling->FirstChildElement(), NULL);
        printxmg.indent();
      }
      else if ( str.compare("UtterRandom") == 0 )
      {
        printxmg.add_comment("UtterRandom");
        this_op = new mem_monop_indiv_t<t_Evaluator, t_Individual>
                         ( evaluator, &t_Evaluator::initialize, std::string( "UtterRandom" ) );
        eostates.storeFunctor( static_cast< TabooOp<t_Individual> *>(this_op) );
      }
      else if ( str.compare("Operators") == 0 )
      {
        if ( sibling->Attribute("type") )
        {
          std::string sstr = sibling->Attribute("type");
          if ( sstr.compare("and") == 0 ) 
          {
            printxmg.add_comment("And Begin");
            SequentialOp<t_Individual> *new_branch = new SequentialOp<t_Individual>;
            this_op = make_genetic_op( *sibling->FirstChildElement(), new_branch);
            if ( not this_op )
            {
              printxmg.remove_last();
              delete new_branch;
            }
            else
            {
              printxmg.add_comment("And End");
              eostates.storeFunctor( new_branch );
            }
          }
        }
        if ( not this_op )
        {
          printxmg.add_comment("Or Begin");
          ProportionalOp<t_Individual> *new_branch = new ProportionalOp<t_Individual>;
          this_op = make_genetic_op( *sibling->FirstChildElement(), new_branch);
          if ( not this_op )
          {
            printxmg.remove_last();
            delete new_branch;
          }
          else
          {
            printxmg.add_comment("Or End");
            eostates.storeFunctor( new_branch );
          }
        }
        is_gen_op = true;
      }
      else // operators specific to t_Individual 
      {
        eoOp<typename t_Individual::t_Object> *eoop = evaluator.LoadGaOp( *sibling );
        // eoop is stored in Wrap_ObjectOp_To_GenOp
        // since eoOp does not derive from eoFunctor_Base, 
        // whereas eoMonOp, eoBinOp, ... do
        // eoFunctor_Base inheritance is necessary to store 
        // as in eostates.storeFunctor( functor )
        if ( eoop )
          this_op = Wrap_ObjectOp_To_GenOp<t_Individual>( eoop, eostates);
        else
        {
          minimizer::Base< typename t_Evaluator::t_Functional > 
            *minimizer = ( minimizer::Base< typename t_Evaluator::t_Functional >* )
                         evaluator.LoadMinimizer( *sibling );
          // MinimizerGenOp owns object pointed to by minimizer
          if ( minimizer )
          {
            minimizers.push_back( minimizer );
            this_op = new MinimizerGenOp<t_Individual, t_Evaluator, t_Population>( *minimizer, *results); 
            eostates.storeFunctor(static_cast< eoGenOp<t_Individual>* >(this_op) );
          }
        }
        is_gen_op = true;
      }
      if ( this_op and sibling->Attribute("period", &period) )
      {
        if (     (    ( not max_generations )
                   or (types::t_unsigned)std::abs(period) < max_generations ) 
             and period > 0 )
        {
          std::ostringstream sstr; 
          sstr << " period = " << period;
          std::string str = sstr.str();
          printxmg.add_to_last(str);
          this_op = new PeriodicOp<t_Individual>( *this_op, (types::t_unsigned) abs(period),
                                              continuator->get_generation_counter(), eostates );
          eostates.storeFunctor( static_cast< PeriodicOp<t_Individual> *>(this_op) );
          is_gen_op = true;
        }
      }
      if ( this_op and current_op )
      {
        if (not sibling->Attribute("prob", &prob) )
          prob = 1.0;
        std::ostringstream sstr;
        sstr << " prob = " << prob;
        std::string str = sstr.str();
        printxmg.add_to_last(str);
        if ( current_op->className().compare("darwin::SequentialOp") == 0 )
          static_cast< SequentialOp<t_Individual>* >(current_op)->add( *this_op,
                                                                   static_cast<double>(prob) );
        else if ( current_op->className().compare("darwin::ProportionalOp") == 0 )
          static_cast< ProportionalOp<t_Individual>* >(current_op)->add( *this_op, 
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

    printxmg.deindent();
    
    return current_op;
  }

  template<class T_INDIVIDUAL, class T_EVALUATOR>
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: make_breeder()
  {
    eoSelectOne<t_Individual> *select;

    select = new eoDetTournamentSelect<t_Individual>(tournament_size);
    breeder = new Breeder<t_Individual>(*select, *breeder_ops, continuator->get_generation_counter() );
    if ( nuclearwinter )
    {
      nuclearwinter->set_op_address( breeder->get_op_address() );
      nuclearwinter->set_howmany( breeder->get_howmany_address() ) ;
    }
    else
      breeder->set_howmany(replacement_rate);

    eostates.storeFunctor(breeder);
    eostates.storeFunctor(select);
  }
  
  template<class T_INDIVIDUAL, class T_EVALUATOR >
  void Darwin<T_INDIVIDUAL, T_EVALUATOR> :: populate ()
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

  template<class T_INDIVIDUAL, class T_EVALUATOR >
  void Darwin<T_INDIVIDUAL, T_EVALUATOR> :: random_populate ( t_Population &_pop, types::t_unsigned _size)
  {
    while ( _pop.size() < _size )
    {
      t_Individual indiv;
      evaluator.initialize(indiv);
      if ( taboos and ( not (*taboos)(indiv) ) )
        _pop.push_back(indiv);
    } // while ( i_pop->size() < target )
    _pop.resize( _size );
  }
  template<class T_INDIVIDUAL, class T_EVALUATOR >
  void Darwin<T_INDIVIDUAL, T_EVALUATOR> :: partition_populate ( t_Population &_pop, types::t_unsigned _size)
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


  template<class T_INDIVIDUAL, class T_EVALUATOR >
  void Darwin<T_INDIVIDUAL,T_EVALUATOR> :: run()
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
    printxmg.flush();
    populate();
    offsprings.clear();
    typename t_Islands :: iterator i_island_begin = islands.begin();
    typename t_Islands :: iterator i_island_end = islands.end();
    typename t_Islands :: iterator i_island;
    for ( i_island = i_island_begin; i_island != i_island_end; ++i_island )
    {
      results->evaluate(offsprings, *i_island); // A first eval of pop.
#ifdef _MPI // "All Gather" new population
      breeder->synchronize_offsprings( *i_island );
#endif 
    }

    do
    {
      try
      {
         i_island = i_island_begin;
         for (int i=0; i_island != i_island_end; ++i, ++i_island )
         {
           types::t_unsigned pSize = i_island->size();
           offsprings.clear(); // new offsprings
           
           (*breeder)(*i_island, offsprings);
           
           results->evaluate(*i_island, offsprings); // eval of parents + offsprings if necessary

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

    printxmg.flush();
    Save();
  }

  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: Load(const std::string &_filename) 
  {
    printxmg.add_comment("new GA run");

    filename = _filename;
    evaluator_filename = _filename;
    restart_filename = _filename;

    // first load filenames
    std::string ga_string;
#ifdef _MPI
    if ( mpi::main.rank() == mpi::ROOT_NODE )
    {
#endif
      OPENXMLINPUT
      const TiXmlElement *parent = docHandle.FirstChild("Job")
                                            .FirstChild("GA").Element();
      if ( not parent )
        return false;
      TIXML_OSTREAM stream;
      stream << *parent;
      ga_string = stream.c_str();
      const TiXmlElement *child = parent->FirstChildElement("Filenames");
      for( ; child; child = child->NextSiblingElement("Filenames") )
      {
        if (     child->Attribute("evaluator") 
             and evaluator_filename != filename )
          evaluator_filename = child->Attribute("evaluator");
        else if (     child->Attribute("restart") 
                  and restart_filename != filename )
          restart_filename = child->Attribute("restart");
      }
      
      // Loads evaluator first 
      if ( evaluator_filename != filename )
      { 
        if ( not evaluator.Load(evaluator_filename) )
          return false;
      } // for following, we know docHandle.FirstChild("Job") exists
      else if ( not evaluator.Load(*docHandle.FirstChild("Job").Element() ) )
        return false;
            
      // finds <GA> ... </GA> block 
      if ( not Load_Parameters( *parent ) )
        return false;

#ifdef _MPI
    }

    // broadcast what we have up to now
    mpi::BroadCast bc( mpi::main );
    broadcast(bc);
    bc.serialize(ga_string);
    bc.allocate_buffers();
    broadcast(bc);
    bc.serialize(ga_string);
    bc();
    broadcast(bc);
    bc.serialize(ga_string);
    bc.reset();


    for( types::t_int i=0; i < mpi::main.size(); ++i )
    {
      if ( mpi::main.rank() == i )
      {
        TiXmlDocument doc;
        doc.Parse( ga_string.c_str() );
        TiXmlHandle docHandle(&doc);
        const TiXmlElement *parent = docHandle.FirstChild("GA").Element();
        if ( not parent )
          return false;
#endif
        Load_History( *parent );
        Load_Taboos( *parent );
        Load_Method( *parent );
        if ( not  Load_Mating( *parent ) )
          return false;
        Load_CheckPoints( *parent );
        
        make_breeder();
        replacement = make_replacement();
#ifdef _MPI
      }
      mpi::main.barrier();
    }
#endif

    results->set_history(history);

#ifdef _MPI
    if ( mpi::main.rank() == mpi::ROOT_NODE )
    {
      OPENXMLINPUT
      const TiXmlElement *parent = docHandle.FirstChild("Job")
                                            .FirstChild("GA").Element();
      if ( not parent )
        return false;
#endif
      
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
        if ( restart_filename != filename )
        { 
          if ( not Restart() )
          {
            std::cerr << "Could not load restart from file" << restart_filename;
            printxmg.add_comment("Starting from scratch");
          }
        } // for following, we know docHandle.FirstChild("Job") exists
        else if ( not Restart( *docHandle.FirstChild("Job").Element() ) )
        {
          std::cerr << "Could not load restart from file" << restart_filename;
          printxmg.add_comment("Starting from scratch");
        }
      }
      restart_xml = parent->FirstChildElement("Save");
      std::ostringstream sstr("Will Save Results");
      if ( restart_xml and restart_xml->Attribute("what") )
      {
        std::string str = restart_xml->Attribute("what");
        if ( str.find("all") != std::string::npos )
        {
          do_save |= SAVE_POPULATION | SAVE_HISTORY;
          sstr << ", Population";
          if ( history )
            sstr << ", and History";
          goto out;
        }
        if ( str.find("pop") != std::string::npos )
        {
          do_save |= SAVE_POPULATION;
          sstr << ", Population";
        }
        if ( str.find("history") != std::string::npos and history)
        {
          do_save |= SAVE_HISTORY;
          sstr << ", History";
        }
      }
out:  printxmg.add_comment( sstr.str() );
#ifdef _MPI
    }

    if ( mpi::main.rank() != mpi::ROOT_NODE )
      do_save = 0;

    if ( history ) history->broadcast(bc);
    results->broadcast(bc);
    broadcast_islands(bc);
    bc.allocate_buffers();
    if ( history ) history->broadcast(bc);
    results->broadcast(bc);
    broadcast_islands(bc);
    bc();
    if( history ) history->broadcast(bc);
    results->broadcast(bc);
    broadcast_islands(bc);
#endif 

    {
      std::ostringstream sstr; 
      sstr << "Mating Tournament Size: " << tournament_size;
      printxmg.add_comment( sstr.str() );
      sstr.str(""); sstr << "Offspring Replacement Rate: " << replacement_rate;
      printxmg.add_comment( sstr.str() );
      sstr.str(""); sstr << "Population Size: " << pop_size;
      printxmg.add_comment( sstr.str() );
      sstr.str("");
      if ( max_generations )
        sstr << "Maximum Number of Generations: " << max_generations;
      else
        sstr << "Unlimited Number of Generations";
      printxmg.add_comment( sstr.str() );
      sstr.str(""); sstr << "Number of Islands: " << nb_islands;
      printxmg.add_comment( sstr.str() );
    }
    return true;
  }


  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: broadcast( mpi::BroadCast &_bc )
  {
#ifdef _MPI
    if( not _bc.serialize( filename ) ) return false;
    if( not _bc.serialize( evaluator_filename ) ) return false;
    if( not _bc.serialize( restart_filename ) ) return false;
    if( not _bc.serialize( tournament_size ) ) return false;
    if( not _bc.serialize( pop_size ) ) return false;
    if( not _bc.serialize( max_generations ) ) return false;
    if( not _bc.serialize( nb_islands ) ) return false;
    if( not _bc.serialize( restarts ) ) return false;
    if( not _bc.serialize( replacement_rate ) ) return false;
    if( not _bc.serialize( do_print_each_call ) ) return false;
    types::t_unsigned n = populate_style;
    if ( not _bc.serialize( n ) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      switch ( n ) 
      {
        case RANDOM_POPULATE: 
          populate_style = RANDOM_POPULATE; break;
        case PARTITION_POPULATE: 
          populate_style = PARTITION_POPULATE; break;
      }
    return _bc.serialize( evaluator );
#else
    return true;
#endif
  }
  template<class T_INDIVIDUAL, class T_EVALUATOR>
  bool Darwin<T_INDIVIDUAL,T_EVALUATOR> :: broadcast_islands( mpi::BroadCast &_bc )
  {
#ifdef _MPI
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
#endif
    return true;
  }

} // namespace darwin


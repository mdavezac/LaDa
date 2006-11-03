#include<eo/eoOpContainer.h>
#include<eo/eoDetTournamentSelect.h>
#include<eo/eoGeneralBreeder.h>
#include<eo/eoReduceMerge.h>


#include "gencount.h"
#include "generator.h"
#include "taboo_minimizers.h"
#include "statistics.h"
#include <lamarck/convex_hull.h>
#include <opt/opt_minimize.h>

#include <math.h>
#include <cmath>

using opt::NO_MINIMIZER;
using opt::LINEAR_MINIMIZER;
using opt::SA_MINIMIZER;

namespace LaDa 
{
# define OPENXMLINPUT \
    TiXmlDocument doc( filename.c_str() ); \
    TiXmlHandle docHandle( &doc ); \
    TiXmlElement *child; \
    if  ( !doc.LoadFile() ) \
    { \
      std::cout << doc.ErrorDesc() << std::endl; \
      throw "Could not load input file in  Darwin<t_Object, t_Lamarck> "; \
    } 

#  define OPENXMGRACE \
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 

#  define CLOSEXMGRACE \
     xmgrace_file.flush(); \
     xmgrace_file.close();

  const t_unsigned svn_revision = 152;
  template<class t_Object, class t_Lamarck> 
    const t_unsigned Darwin<t_Object, t_Lamarck> :: DARWIN  = 0;
  template<class t_Object, class t_Lamarck> 
    const t_unsigned Darwin<t_Object, t_Lamarck> :: LAMARCK = 1;
  template<class t_Object, class t_Lamarck> 
    const t_unsigned Darwin<t_Object, t_Lamarck> :: DEBUG   = 2;
    

  template< class t_Object, class t_Lamarck >
  Darwin<t_Object, t_Lamarck> :: Darwin ( t_Lamarck *_lam )
  {
    lamarck = _lam;

    // default values
    crossover_value = 0.5;
    mutation_value = 0.05;
    tournament_size = 2;
    replacement_rate = 0.1;
    max_generations = 200;
    pop_size = 100;
    method = DARWIN;
    evolve_from_start = false; 
    multistart = false; 
    minimizer = NO_MINIMIZER;
    max_calls = UINT_MAX;
    is_one_point_hull = false;
    minimize_best = 0;
    minimize_best_every = 5;

    extra_popalgo = NULL;
    taboos = NULL;
    nuclearwinter = NULL;
    agetaboo = NULL;
    colonize = NULL;

    print_strings.reserve(10);
    t_Object :: is_using_phenotype = false;
    nb_islands = 1;
  }

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: run()
  {
    make_algo();
    populate();
    offsprings.clear();
    typename t_Islands :: iterator i_island_begin = islands.begin();
    typename t_Islands :: iterator i_island_end = islands.end();
    typename t_Islands :: iterator i_island;
    for ( i_island = i_island_begin; i_island != i_island_end; ++i_island )
      (*popEval)(offsprings, *i_island); // A first eval of pop.

    do
    {
      try
      {
         i_island = i_island_begin;
         for (int i=0; i_island != i_island_end; ++i, ++i_island )
         {
           t_unsigned pSize = i_island->size();
           offsprings.clear(); // new offsprings
           
           (*breeder)(*i_island, offsprings);
           
           (*popEval)(*i_island, offsprings); // eval of parents + offsprings if necessary
          
           (*replace)(*i_island, offsprings); // after replace, the new pop. is in population
           
           if ( extra_popalgo )
             (*extra_popalgo)(*i_island); // minimizes best for instance
           
           if (pSize > i_island->size())
               throw std::runtime_error("Population shrinking!");
           else if (pSize < i_island->size())
               throw std::runtime_error("Population growing!");
         }

         // creates colonies if requested
         if ( colonize )
           (*colonize)( islands );

      }
      catch (std::exception& e)
      {
            std::string s = e.what();
            s.append( " in eoEasyEA");
            throw std::runtime_error( s );
      }
    } while ( continuator->apply( i_island_begin, i_island_end ) );
  }

  template<class t_Object, class t_Lamarck >
  bool Darwin<t_Object, t_Lamarck> :: Load(const std::string &_filename) 
  {
    filename = _filename;

    OPENXMLINPUT

    TiXmlElement *element = docHandle.FirstChild("LaDa")
                                     .FirstChild("GA").Element();
    TiXmlElement *parent;
    if (not element) 
      return false;


    { // first, opens output file
      xmgrace_filename = "convex_hull.agr";
      child = docHandle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
      if ( child and child->Attribute("xmgrace") )
        xmgrace_filename = child->Attribute("xmgrace");
    }
    
    // finds if in <GA> ... </GA> block 
    {
      std::string str = element->Value();
      parent = element;
      if ( str.compare("GA" ) != 0 )
        parent = element->FirstChildElement("GA");
    }
    

    // tournament size when selecting parents
    child = parent->FirstChildElement( "Selection" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("value", &d) and d > 1 )
        tournament_size = abs(d);
    }

    // Offsprings
    child = parent->FirstChildElement( "Offsprings" );
    if ( child )
    {
      // rate
      double d = 0;
      if ( child->Attribute("rate", &d) )
        if ( d <= 1.0 and d > 0.0 )
          replacement_rate = types::t_real(d);
    }

    // population size
    child = parent->FirstChildElement( "Population" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("size", &d) )
        if ( d > 0 )
          pop_size = (eotypes::t_unsigned) abs(d);
    }

    // method and nb steps
    {
      int d = 0;
      if ( parent->Attribute("maxgen", &d) )
        if ( d > 0 )
          max_generations = (types::t_unsigned) abs(d);

      if ( parent->Attribute("method") )
      {
        std::string str = parent->Attribute("method");
        if ( str.find("multistart") != std::string::npos )
        {
          multistart = true;
          eoHowMany nb(replacement_rate);
          max_generations = (types::t_unsigned) pop_size
                            + (types::t_unsigned) nb(pop_size) * max_generations;
          pop_size = 1;
          replacement_rate = 1.0;
          evolve_from_start = true;
        }
      } // if attribute "method" exists
      // number of Islands
      if ( parent->Attribute("islands", &d ) )
        nb_islands = ( d > 0 ) ? (types::t_int) abs(d) : 1;
    }   

      
    // Phenotype vs Genotype
    child = parent->FirstChildElement( "Phenotype" );
    if ( child )
      t_Object :: is_using_phenotype = true;

    write_xmgrace_header();

    // reseed
    { 
      t_int seed;
      child = parent->FirstChildElement("Seed");
      if (child and child->Attribute("nb", &seed) )
      {
        OPENXMGRACE
        xmgrace_file << "# Seed: " << seed << std::endl;
        rng.reseed( abs(seed) );
        CLOSEXMGRACE
      }
    }

    return true;
  }

  template<class t_Object, class t_Lamarck > 
  MinimizationOp< t_Object, Darwin<t_Object, t_Lamarck> >* 
    Darwin<t_Object, t_Lamarck> :: Load_Minimizer( const TiXmlElement* el, std::ofstream &_f )
  {
    t_unsigned type = LINEAR_MINIMIZER;
    if ( el->Attribute("type") )
    {
      std::string str = el->Attribute("type");
      if ( str.compare("linear" ) == 0 ) // Wang
      {
        _f << "Linear ";
        type = LINEAR_MINIMIZER;
      }
      else if ( str.compare("SA" ) == 0 ) // Wang
      {
        _f << "SA ";
        type = SA_MINIMIZER;
      }
    }
    t_unsigned n = UINT_MAX;
    int i = 0;
    if ( el->Attribute("maxeval", &i) )
      n = ( i <= 0 ) ? UINT_MAX : (types::t_unsigned) abs(i);
    if ( type == SA_MINIMIZER or type == LINEAR_MINIMIZER )
    {
      if ( n == UINT_MAX ) 
        _f << "with unlimited evaluations ";
      else
        _f << "with at most " << n << " evaluations ";
    }


    MinimizationOp< t_Object, t_Darwin >*  minop 
      = new MinimizationOp< t_Object, t_Darwin >( lamarck->add_minimizer( type, n), *this );
    eostates.storeFunctor(minop);
    return minop;
  }


  // makes genetic operators
  template< class t_Object, class t_Lamarck >
  eoGenOp<t_Object>* Darwin<t_Object, t_Lamarck>
      :: make_genetic_op( const TiXmlElement &el, std::ofstream &_f,
                          std::string &_special, std::string &_base,
                          eoGenOp<t_Object> *current_op = NULL)
  {
    eoOp<t_Object>* this_op;
    const TiXmlElement *sibling = &el;
    if (not sibling)
      throw "Error while creating Operator ";

    for ( ; sibling; sibling = sibling->NextSiblingElement() )
    {
      std::string str = sibling->Value();
      double prob = 0.0;
      int period = 0;
      this_op = NULL;
      bool is_gen_op = false;
      

      // then creates sibling
      if ( str.compare("Crossover" ) == 0 )
      {
        double d; 
        sibling->Attribute("value", &d);
        if ( d <= 0 and d > 1 )
          d = double(crossover_value);
        this_op = new Crossover<t_Object>( types::t_real(d) );
        eostates.storeFunctor( static_cast< Crossover<t_Object> *>(this_op) );
        _f << "# " << _special << _base << "Crossover: value=" << d;
      }
      else if ( str.compare("Krossover" ) == 0 )
      {
        this_op = new Krossover<t_Lamarck, t_Object>( *lamarck );
        eostates.storeFunctor( static_cast< Krossover<t_Lamarck, t_Object> *>(this_op) );
        _f << "# " << _special << _base << "Krossover ";
      }
      else if ( str.compare("Mutation" ) == 0 )
      {
        double d; 
        sibling->Attribute("value", &d);
        if ( d <= 0 and d > 1 )
          d = 1.0 / (double) lamarck->get_pb_size();
        this_op = new Mutation<t_Object>(  types::t_real(d) );
        eostates.storeFunctor( static_cast< Mutation<t_Object> *>(this_op) );
        _f << "# " << _special << _base << "Mutation: value=" << d;
      }
      else if ( str.compare("Minimizer") == 0 )
      {
        _f << "# " << _special << _base << "Minimizer: ";
        this_op = Load_Minimizer( sibling, _f );
      }
      else if ( str.compare("TabooMinimizer") == 0 and taboos )
      {
        std :: string type = "SA";
        if ( sibling->Attribute("type") )
          type = sibling->Attribute("type");
        if ( type.compare("GradientSA") != 0 )
          type = "SA";
        int i=0; t_unsigned maxeval = UINT_MAX;
        if ( sibling->Attribute("maxeval", &i ) )
          maxeval = ( i > 0 ) ? (types::t_unsigned) abs(i) : UINT_MAX;
        _f << "# " << _special << _base << "TabooMinimizer: " 
           << type << " maxeval " << maxeval;
        if ( type.compare("SA") == 0 )
          this_op = new SA_TabooOp<t_Darwin>( *this, *taboos, maxeval );
        else
          this_op = new GradientSA_TabooOp<t_Darwin>( *this, *taboos, maxeval );

        eostates.storeFunctor( static_cast< eoMonOp<t_Object> *>(this_op) );
      }
      else if ( str.compare("UtterRandom") == 0 )
      {
        this_op = new UtterRandom<t_Object>;
        eostates.storeFunctor( static_cast< UtterRandom<t_Object> *>(this_op) );
        _f << "# " << _special << _base << "UtterRandom ";
      }
      else if ( str.compare("TabooOp") == 0  and taboos )
      {
        _f << "# " << _special << _base << "TabooOp begin " << std::endl;
        std :: string special = _special + _base;
        eoGenOp<t_Object> *taboo_op;
        taboo_op = make_genetic_op( *sibling->FirstChildElement(), _f,  special, _base, NULL);
        _f << "# " << _special << _base << "TabooOp end";
        this_op = new TabooOp<t_Object> ( *taboo_op, *taboos, 
                                          (types::t_unsigned)(pop_size+1),
                                          eostates );
        eostates.storeFunctor( static_cast< TabooOp<t_Object> *>(this_op) );
        is_gen_op = true;
      }
      else if ( str.compare("TabooOp") == 0 )
        this_op = make_genetic_op( *sibling->FirstChildElement(), _f,  _special, _base, NULL);
      else if ( str.compare("Operators") == 0 )
      {
        if (     sibling->Attribute("type") )
        {
          std::string sstr = sibling->Attribute("type");
          if ( sstr.compare("and") == 0 ) 
          {
            _f << "# " << _special << _base << "And begin " << std::endl;
            std :: string special = _special + _base;
            eoSequentialOp<t_Object> *new_branch = new eoSequentialOp<t_Object>;
            eostates.storeFunctor( new_branch );
            this_op = make_genetic_op( *sibling->FirstChildElement(), _f,  special, _base, new_branch);
            _f << "# " << _special << _base << "And end";
          }
        }
        if ( not this_op )
        {
          _f << "# " << _special << _base << "Or begin " << std::endl;
          std :: string special = _special + _base;
          eoProportionalOp<t_Object> *new_branch = new eoProportionalOp<t_Object>;
          eostates.storeFunctor( new_branch );
          this_op = make_genetic_op( *sibling->FirstChildElement(), _f,  special, _base, new_branch);
          _f << "# " << _special << _base << "Or end";
        }
        is_gen_op = true;
      }
      if ( this_op and sibling->Attribute("period", &period) )
      {
        if (period > 0 and abs(period) < max_generations )
        {
          _f << " period= " << prob;
          this_op = new PeriodicOp<t_Object>( *this_op, (types::t_unsigned) abs(period),
                                              continuator->get_generation_counter(), eostates );
          eostates.storeFunctor( static_cast< PeriodicOp<t_Object> *>(this_op) );
          is_gen_op = true;
        }
      }
      if ( this_op and current_op != NULL )
      {
        if (not sibling->Attribute("prob", &prob) )
          prob = 1.0;
        _f << " prob= " << prob << std::endl;
        if ( current_op->className().compare("SequentialOp") == 0 )
          static_cast< eoSequentialOp<t_Object>* >(current_op)->add( *this_op, types::t_real(prob) );
        else if ( current_op->className().compare("ProportionalOp") == 0 )
          static_cast< eoProportionalOp<t_Object>* >(current_op)->add( *this_op, types::t_real(prob) );
      }
      else if ( this_op )
      {
        if ( is_gen_op )
          current_op = static_cast<eoGenOp<t_Object>*> (this_op);
        else 
          current_op = &wrap_op<t_Object>(*this_op, eostates);
        _f << std::endl;
      }
    }
    
    if ( not current_op  )
    {
      std::cerr << " Error while creating genetic operators " << std::endl;
      throw;
    }
    return current_op;
  }

  template<class t_Object, class t_Lamarck>
  eoBreed<t_Object>* Darwin<t_Object, t_Lamarck> :: make_breeder()
  {
    eoBreed<t_Object> *breed;
    eoSelectOne<t_Object> *select;

    select = new eoDetTournamentSelect<t_Object>(tournament_size);
    if ( nuclearwinter )
    {
      breed = new Breeder<t_Object>(*select, *breeder_ops, continuator->get_generation_counter() );
      nuclearwinter->set_op_address( static_cast<Breeder<t_Object>*>(breed)->get_op_address() );
      nuclearwinter->set_howmany( static_cast<Breeder<t_Object>*>(breed)->get_howmany_address() ) ;
    }
    else
      breed = new eoGeneralBreeder<t_Object>(*select, *breeder_ops, replacement_rate);

    eostates.storeFunctor(breed);
    eostates.storeFunctor(select);

    return breed;
  }
  
  template<class t_Object, class t_Lamarck>
  eoReplacement<t_Object>* Darwin<t_Object, t_Lamarck> :: make_replacement()
  {
    eoTruncate<t_Object>* truncate = new  eoTruncate<t_Object>;
    eoMerge<t_Object>* merge = new  eoPlus<t_Object>;
    eoReduceMerge<t_Object>* reducemerge = new eoReduceMerge<t_Object>( *truncate, *merge );
    eostates.storeFunctor(truncate);
    eostates.storeFunctor(merge);
    eostates.storeFunctor(reducemerge);

    return reducemerge;
  }

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: make_algo()
  {
    evaluation = new Evaluation<t_Darwin>(*this);
    popEval = new EvaluatePop<t_Object>(*evaluation);
    
    make_taboos();      // order counts
    make_breeder_ops(); // order counts
    make_checkpoint();  // order counts

    breeder = make_breeder();
    replace = make_replacement();
    make_extra_algo();
    eostates.storeFunctor(evaluation);
    eostates.storeFunctor(popEval);
  }

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: make_checkpoint()
  {
    OPENXMLINPUT;
    OPENXMGRACE;

    // contains all checkpoints
    continuator = new IslandsContinuator<t_Object>(max_generations);
    eostates.storeFunctor( continuator );
    GenCount &generation_counter = continuator->get_generation_counter();
 
    // our very own updater wrapper to print stuff
    PrintXmgrace< t_Darwin > *printxmgrace = new PrintXmgrace< t_Darwin >(this);
    eostates.storeFunctor(printxmgrace);
    continuator->add(*printxmgrace);
    
    // some statistics
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Statistics").Element();
    for( ; child; child = child->NextSiblingElement("Statistics") )
    {
      std::string str = "accumulated";
      eoStatBase<t_Object> *average = NULL;

      if ( child->Attribute("type" ) )
        str = child->Attribute( "type" );

      if ( str.compare("accumulated") == 0 )
      {
        average = new AccAverage< t_Darwin >( *this, generation_counter,
                                              lamarck->get_pb_size() );
        xmgrace_file << "# Statistics: accumulated average in real space" << std::endl;
      }
      else if ( str.compare("population") == 0 )
      {
        xmgrace_file << "# Statistics: population average";
        std::string att = "false";
        if ( child->Attribute("kspace") )
          att = child->Attribute("kspace");
        if ( att.compare("true") == 0 )
        {
          average = new kPopAverage< t_Darwin >( *this );
          xmgrace_file << "in k-space" << std::endl;
        }
        else 
          average = new PopAverage< t_Darwin >( *this, lamarck->get_pb_size() );
        xmgrace_file << std::endl;
      }
      else if ( str.compare("diversity") == 0 )
      {
        std::string att = "non-identical";
        if ( child->Attribute("count") )
          att = child->Attribute("identical");
        xmgrace_file << "# Statistics: diversity " << std::endl;
        if ( att.compare("identical") == 0 )
        {
          average = new Diversity< t_Darwin >( *this, false );
          xmgrace_file << " including count over identical object " << std::endl;
        }
        else
          average = new Diversity< t_Darwin >( *this, true );
        xmgrace_file << std::endl;
      }
      else if ( str.compare("true census") == 0 )
      {
        average = new TrueCensus< t_Darwin >( *this );
        xmgrace_file << "# Statistics: True population size, discounting twins" << std::endl;
      }

      if ( not average  )
      {
        std::cerr << "Error while reading Statistics tags " << std::endl;
        throw;
      }
      eostates.storeFunctor( average );
      continuator->add( *average );
    }

    // Creates Age taboo
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Taboos")
                     .FirstChild("AgeTaboo").Element();
    if ( child and agetaboo )
    {
      int d = 0;
      child->Attribute("lifespan", &d );
      t_unsigned length = ( d >=0 ) ? (types::t_int) abs(d) : UINT_MAX;
      bool print_out = false;
      if ( child->Attribute("printout") )
      {
        std::string str = child->Attribute("printout");
        if ( str.compare("true") == 0 )
          print_out = true;
      }
      UpdateAgeTaboo< t_Darwin > *updateagetaboo
             = new UpdateAgeTaboo<t_Darwin > ( *agetaboo, generation_counter,
                                               *this, length, print_out);
      xmgrace_file << "# Age Taboo, lifespan=" << d << std::endl;
      eostates.storeFunctor(updateagetaboo);
      continuator->add(*updateagetaboo);
    }

    // Nuclear Winter -- only for agetaboo
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Taboos")
                     .FirstChild("NuclearWinter").Element();
    if ( child and agetaboo)
    {
      eoGenOp<t_Object> *nuclear_op;
    
      // first creates the nuclear op from input
      std::string str = "  ";
      xmgrace_file << "# Nuclear Operator begin " << std::endl;
      nuclear_op = make_genetic_op( *child->FirstChildElement(), xmgrace_file, str,str );
      xmgrace_file << "# Nuclear Operator end " << std::endl;
      if ( not nuclear_op )
        throw "Error while creating operators in  Darwin<t_Object, t_Lamarck>  :: make_GenOp ";
      
      // creates the NuclearWinter 
      nuclearwinter = new NuclearWinter<t_Darwin >
                                       ( *taboos, *breeder_ops, *nuclear_op, *this,
                                         replacement_rate );
      xmgrace_file << "# NuclearWinter " << std::endl;
      eostates.storeFunctor( nuclearwinter );
      continuator->add(*nuclearwinter);
    }

    CLOSEXMGRACE
  } // end of make_check_point

    // create Taboos

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: make_taboos()
  {
    OPENXMLINPUT
    OPENXMGRACE 

    IslandsTaboos<t_Object> *poptaboo = NULL;
    Taboo<t_Object> *offspringtaboo = NULL;
    agetaboo = NULL;

    // creates age taboo -- updater is created in make_checkpoint
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Taboos")
                     .FirstChild("AgeTaboo").Element();
    if (child)
    {
      agetaboo = new Taboo< t_Object, std::list<t_Object> >;
      eostates.storeFunctor(agetaboo);
    }
    
    // creates pop taboo
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Taboos")
                     .FirstChild("PopTaboo").Element();
    if (child)
    {
      xmgrace_file << "# Pop Taboo " << std::endl; 
      poptaboo = new IslandsTaboos<t_Object>( islands );
      eostates.storeFunctor(poptaboo);
      if ( agetaboo ) 
      {
        taboos = new Taboos<t_Object>;
        eostates.storeFunctor(taboos);
        static_cast< Taboos<t_Object>* >(taboos)->add( agetaboo );
        static_cast< Taboos<t_Object>* >(taboos)->add( poptaboo );
      }
    }
    
    // creates offspring taboo
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Taboos")
                     .FirstChild("OffspringTaboo").Element();
    if (child)
    {
      xmgrace_file << "# Offspring Taboo " << std::endl; 
      offspringtaboo = new Taboo<t_Object>( &offsprings );
      eostates.storeFunctor(offspringtaboo);
      if ( (not taboos) and (agetaboo or poptaboo) ) 
      {
        taboos = new Taboos<t_Object>;
        eostates.storeFunctor(taboos);
        if ( agetaboo )
          static_cast< Taboos<t_Object>* >(taboos)->add( agetaboo );
        if ( poptaboo )
          static_cast< Taboos<t_Object>* >(taboos)->add( poptaboo );
      }
    }
    if ( taboos )
      static_cast< Taboos<t_Object>* >(taboos)->add( offspringtaboo );
    else if (agetaboo) 
      taboos = agetaboo;
    else if (poptaboo) 
      taboos = poptaboo;
    else if (offspringtaboo) 
      taboos = offspringtaboo;

    CLOSEXMGRACE
  }

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: make_colonize()
  {
    OPENXMLINPUT
    OPENXMGRACE 

    // creates age taboo -- updater is created in make_checkpoint
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("Colonize").Element();

    t_int every = 0;
    if (not child and not child->Attribute("every", &every) )
      return;
    if ( every <= 0 and every >= max_generations )
      return;
    xmgrace_file << "# Colonize every " << every;

    bool is_pop_stable = false;
    if ( child->Attribute("stablepop") )
    {
      xmgrace_file << " keeping population size stable ";
      is_pop_stable = true;
    }
    xmgrace_file << std::endl;

    colonize = new Colonize<t_Object>( evaluation, breeder,
                                       abs(every), is_pop_stable );
    eostates.storeFunctor(colonize);
    
    CLOSEXMGRACE
  }

  template<class t_Object, class t_Lamarck>
  void Darwin<t_Object, t_Lamarck> :: make_breeder_ops ()
  {
    OPENXMLINPUT
    OPENXMGRACE

    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA").Element();
    xmgrace_file << "# Breeding Operator begin " << std::endl;
    std::string str = "  ";
    breeder_ops = make_genetic_op(*child->FirstChildElement(), xmgrace_file, str, str);
    xmgrace_file << "# Breeding Operator end " << std::endl;
    if ( not breeder_ops )
      throw "Error while creating operators in  Darwin<t_Object, t_Lamarck>  :: make_GenOp ";

    CLOSEXMGRACE
  }

  template<class t_Object, class t_Lamarck>
  void Darwin<t_Object, t_Lamarck> :: populate ()
  {
    Generator generator;

    t_Object indiv;
    indiv.resize( lamarck->get_pb_size() );
    indiv.set_age( continuator->age() ); 

    eoPop<t_Object> population;
    population.reserve(pop_size);

    for( eotypes::t_unsigned n = 0; n < nb_islands; ++n )
    {
      for( eotypes::t_unsigned i = 0; i < pop_size; ++i )
      {
        typename t_Object :: iterator i_var = indiv.begin();
        typename t_Object :: iterator i_end = indiv.end();
        for ( ; i_var != i_end; ++i_var)
          *i_var = generator();
        if ( t_Object :: is_using_phenotype )
          indiv.set_phenotype_to_genotype();
        population.push_back(indiv);
      }
      islands.push_back( population );
      population.clear();
    }
  }

  template< class t_Object, class t_Lamarck >
  void Darwin <t_Object, t_Lamarck> :: write_xmgrace_header()
  {
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    xmgrace_file << "# LaDa svn revision: " << svn_revision << std::endl;
    xmgrace_file << "# Number of Islands: " << nb_islands << std::endl;
    xmgrace_file << "# population size: " << pop_size << std::endl;
    xmgrace_file << "# replacement rate: " << replacement_rate << std::endl;
    xmgrace_file << "# max generations: " << max_generations << std::endl;
    xmgrace_file << "# Using Phenotype: ";
    if ( t_Object :: is_using_phenotype )
      xmgrace_file << "true" << std::endl;
    else
      xmgrace_file << "false" << std::endl;
    if ( minimize_best > 0 and minimize_best <= 1 )
    {
      xmgrace_file << "# minimize best: rate " << minimize_best
                   << " every " << minimize_best_every << std::endl;
    }
    lamarck->write_xmgrace_header( xmgrace_file );
    xmgrace_file.flush();
    xmgrace_file.close();
  }

    
  template< class t_Object, class t_Lamarck >
  void Darwin <t_Object, t_Lamarck> :: print_xmgrace( bool is_last_call = false )
  {
    OPENXMGRACE 

    bool print_ch = not t_Object :: is_baseline_valid();
    std::string special = " ";
    if ( not print_ch )
      special = " ? ";
    if ( is_last_call )
      xmgrace_file << " # last call for alcohol" << std::endl; 
    std::vector< std::string > :: const_iterator i_str = print_strings.begin();
    std::vector< std::string > :: const_iterator i_end = print_strings.end();
    if ( i_str != i_end )
    {
      xmgrace_file << " #" << special <<  "iteration:" << continuator->age() << std::endl; 
      for ( ; i_str != i_end; ++i_str )
        xmgrace_file << " #" << special << (*i_str) << std::endl;
      print_strings.clear();
    }
    else if ( print_ch or is_last_call )
      xmgrace_file << " #" << special <<  "iteration:" << continuator->age() << std::endl; 

    if ( print_ch or is_last_call )
    {
      xmgrace_file << " # " << special << "TabooGradientSA: " 
                   << GradientSA_TabooOp<t_Darwin>::nb_evals << " "
                   << GradientSA_TabooOp<t_Darwin>::nb_grad_evals << std::endl;
      xmgrace_file << " # " << special << "TabooSA: " 
                   << SA_TabooOp<t_Darwin>::nb_evals << std::endl;
      lamarck->print_xmgrace( xmgrace_file,  print_ch );
      t_Object :: validate_baseline();
    }

    CLOSEXMGRACE
  }


  template< class t_Object, class t_Lamarck >
  void Darwin <t_Object, t_Lamarck> :: make_extra_algo()
  {
    OPENXMLINPUT
    OPENXMGRACE

    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("PopAlgo").Element();
    if ( not child )
      return;
    
    // get parameters
    {
      minimize_best = 0.1;
      double d=0;
      if ( child->Attribute( "rate", &d ) )
        if ( d > 0 and d <= 1 )
          minimize_best = types::t_real(d);
      int u = 0;
      if ( child->Attribute( "every", &u ) )
        minimize_best_every = ( u > 0 and (types::t_unsigned) abs(u) <= max_generations ) ?
                              (types::t_unsigned) abs(u)  : 0 ;
      std::string str;
    }

    xmgrace_file << "# PopAlgorithm: rate="<< minimize_best 
                 << " every=" << minimize_best_every << std::endl;
    std::string str = "    ";
    std::string base = "  ";
    eoGenOp<t_Object> *op = make_genetic_op(*child->FirstChildElement(), xmgrace_file, str, base);

    if ( not op )
    {
      std::cout << "Error while creating Genetic Operator for ExtraPopAlgo " << std::endl;
      throw; 
    }
    extra_popalgo = new Extra_PopAlgo< t_Darwin > 
                                     ( *op, *this, minimize_best, minimize_best_every );
    eostates.storeFunctor( extra_popalgo );

    CLOSEXMGRACE
  }

} // namespace LaDa


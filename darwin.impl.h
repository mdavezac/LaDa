#include<eo/eoOpContainer.h>
#include<eo/eoDetTournamentSelect.h>
#include<eo/eoGeneralBreeder.h>
#include<eo/eoReduceMerge.h>


#include "checkpoint.h"
#include "generator.h"
#include <lamarck/convex_hull.h>
#include <opt/opt_minimize.h>
using opt::NO_MINIMIZER;
using opt::WANG_MINIMIZER;
using opt::PHYSICAL_MINIMIZER;
using opt::LINEAR_MINIMIZER;
using opt::SA_MINIMIZER;

namespace LaDa 
{
  template<class t_Object, class t_Lamarck> 
    const unsigned Darwin<t_Object, t_Lamarck> :: DARWIN  = 0;
  template<class t_Object, class t_Lamarck> 
    const unsigned Darwin<t_Object, t_Lamarck> :: LAMARCK = 1;
  template<class t_Object, class t_Lamarck> 
    const unsigned Darwin<t_Object, t_Lamarck> :: DEBUG   = 2;

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
  }

  template< class t_Object, class t_Lamarck >
  void Darwin<t_Object, t_Lamarck> :: run()
  {
    eoPop<t_Object> offsprings;

    make_algo();
    populate();
    popEval->operator()(offsprings, population); // A first eval of pop.

    do
    {
      try
      {
         unsigned pSize = population.size();
         offsprings.clear(); // new offsprings

         breed->operator()(population, offsprings);

         popEval->operator()(population, offsprings); // eval of parents + offsprings if necessary

         replace->operator()(population, offsprings); // after replace, the new pop. is in population

         if ( extra_popalgo )
           extra_popalgo->operator()(population); // minimizes best for instance

         if (pSize > population.size())
             throw std::runtime_error("Population shrinking!");
         else if (pSize < population.size())
             throw std::runtime_error("Population growing!");

      }
      catch (std::exception& e)
      {
            std::string s = e.what();
            s.append( " in eoEasyEA");
            throw std::runtime_error( s );
      }
    } while ( continuator->operator()( population ) );
  }

  template<class t_Object, class t_Lamarck >
  bool Darwin<t_Object, t_Lamarck> :: Load(const std::string &_filename) 
  {
    filename = _filename;
    TiXmlDocument doc( filename.c_str() );
    
    if  ( !doc.LoadFile() )
    {
      std::cerr << "error while opening input file " << filename << std::endl
                << doc.ErrorDesc() << std::endl; 
      return false;
    }

    TiXmlHandle docHandle( &doc );
    TiXmlElement *element = docHandle.FirstChild("LaDa")
                                     .FirstChild("GA").Element();
    TiXmlElement *child, *parent;
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
        tournament_size = d;
    }

    // Offsprings
    child = parent->FirstChildElement( "Offsprings" );
    if ( child )
    {
      // rate
      double d = 0;
      if ( child->Attribute("rate", &d) )
        if ( d <= 1.0 and d > 0.0 )
          replacement_rate = d;
    }

    // population size
    child = parent->FirstChildElement( "Population" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("size", &d) )
        if ( d > 0 )
          pop_size = d;
    }

    // method and nb steps
    {
      int d = 0;
      if ( parent->Attribute("maxgen", &d) )
        if ( d > 0 )
          max_generations = d;

      if ( parent->Attribute("method") )
      {
        std::string str = parent->Attribute("method");
        if ( str.find("multistart") != std::string::npos )
        {
          multistart = true;
          eoHowMany nb(replacement_rate);
          max_generations = pop_size + nb(pop_size) * max_generations;
          pop_size = 1;
          replacement_rate = 1.0;
          evolve_from_start = true;
        }
      } // if attribute "method" exists
    }   

      

    write_xmgrace_header();

    return true;
  }

  template<class t_Object, class t_Lamarck > 
  MinimizationOp< t_Object, Darwin<t_Object, t_Lamarck> >* 
    Darwin<t_Object, t_Lamarck> :: Load_Minimizer( const TiXmlElement* el, std::ofstream &_f )
  {
    unsigned type = LINEAR_MINIMIZER;
    if ( el->Attribute("type") )
    {
      std::string str = el->Attribute("type");
      if ( str.compare("wang" ) == 0 ) // Wang
      {
        _f << "Wang ";
        type = WANG_MINIMIZER;
      }
      else if ( str.compare("physical" ) == 0 ) // Wang
      {
        _f << "Physical ";
        type = PHYSICAL_MINIMIZER;
      }
      else if ( str.compare("linear" ) == 0 ) // Wang
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
    unsigned n = UINT_MAX;
    int i = 0;
    if ( el->Attribute("maxeval", &i) )
      n = ( i <= 0 ) ? UINT_MAX : abs(i);
    if ( type == SA_MINIMIZER or type == LINEAR_MINIMIZER )
    {
      if ( n == UINT_MAX ) 
        _f << "with unlimited evaluations ";
      else
        _f << "with at most " << n << " evaluations ";
    }


    MinimizationOp< t_Object, Darwin<t_Object, t_Lamarck> >*  minop 
      = new MinimizationOp< t_Object, Darwin<t_Object, t_Lamarck> >( lamarck->add_minimizer( type, n), *this );
    eostates.storeFunctor(minop);
    return minop;
  }


  // makes genetic operators
  template< class t_Object, class t_Lamarck >
  eoGenOp<t_Object>* Darwin<t_Object, t_Lamarck> :: make_GenOp( const TiXmlElement *_el,
                                                                std::ofstream  &_f)
  {
    // creates operator from recurrent input
    std::string str = "  ";
    eoGenOp<t_Object>* this_op = make_recurrent_op(*_el, _f, str, str);
    if ( not this_op )
      throw "Error while creating operators in  Darwin<t_Object, t_Lamarck>  :: make_GenOp ";

    return  this_op;
  }

  template< class t_Object, class t_Lamarck >
  eoGenOp<t_Object>* Darwin<t_Object, t_Lamarck> :: make_recurrent_op(const TiXmlElement &el,
                                                                      std::ofstream &_f,
                                                                      std::string &_special,
                                                                      std::string &_base)
  {
    eoGenOp<t_Object>* this_op;
    bool is_sequential;
    if ( el.Attribute( "type" ) )
    {
       std::string str =  el.Attribute( "type" );
       is_sequential = true;
       if ( str.compare("and" ) == 0 ) // operators applied sequentially
       {
         double d;
         if ( el.Attribute("prob", &d ) )
           _f << "# " << _special << "Sequential: prob " << d << std::endl;
         else
           _f << "# " << _special << "Sequential" << std::endl;
         this_op = new  eoSequentialOp<t_Object>;
       }
       else
       {
         is_sequential = false;
         double d;
         if ( el.Attribute("prob", &d ) )
           _f << "# " << _special << "Proportional: prob " << d << std::endl;
         else
           _f << "# " << _special << "Proportional" << std::endl;
         this_op = new  eoProportionalOp<t_Object>;
       }
       eostates.storeFunctor(this_op);
    }
    else
      throw "Error while creating Operator ";

    const TiXmlElement *child = el.FirstChildElement();
    if (not child)
      throw "Error while creating Operator ";

    for ( ; child; child = child->NextSiblingElement() )
    {
      std::string str = child->Value();
      double prob = 0.0;
      
      // gets probability for applying child 
      if ( not child->Attribute("prob", &prob) )
        prob = 1.0;

      // then creates child
      if ( str.compare("Crossover" ) == 0 )
      {
        double d; 
        child->Attribute("value", &d);
        if ( d <= 0 and d > 1 )
          d = crossover_value;
        Crossover<t_Object>* crossover = new Crossover<t_Object>( d );
        eostates.storeFunctor(crossover);
        _f << "# " << _special << _base << "Crossover: value=" << d 
           << " prob="<< prob << std::endl;
        if ( is_sequential )
          static_cast< eoSequentialOp<t_Object>* >(this_op)->add( *crossover, prob );
        else
          static_cast< eoProportionalOp<t_Object>* >(this_op)->add( *crossover, prob );
      }
      else if ( str.compare("Mutation" ) == 0 )
      {
        double d; 
        child->Attribute("value", &d);
        if ( d <= 0 and d > 1 )
          d = 1.0 / (double) lamarck->get_pb_size();
        Mutation<t_Object> *mutation = new Mutation<t_Object>( d );
        eostates.storeFunctor(mutation);
        _f << "# " << _special << _base << "Mutation: value=" << d 
           << " prob="<< prob << std::endl;
        if ( is_sequential )
          static_cast< eoSequentialOp<t_Object>* >(this_op)->add( *mutation, prob );
        else
          static_cast< eoProportionalOp<t_Object>* >(this_op)->add( *mutation, prob );
      }
      else if ( str.compare("Minimizer") == 0 )
      {
        _f << "# " << _special << _base << "Minimizer: ";
        eoMonOp<t_Object>* minop = Load_Minimizer( child, _f );
        _f << ", prob=" << prob << std::endl;
        if ( is_sequential )
          static_cast< eoSequentialOp<t_Object>* >(this_op)->add( *minop, prob );
        else
          static_cast< eoProportionalOp<t_Object>* >(this_op)->add( *minop, prob );
      }
      else if ( str.compare("UtterRandom") == 0 )
      {
        UtterRandom<t_Object>* utterrandom = new UtterRandom<t_Object>;
        eostates.storeFunctor(utterrandom);
        _f << "# " << _special << _base << "UtterRandom "
           << " prob "<< prob << std::endl;
        if ( is_sequential )
          static_cast< eoSequentialOp<t_Object>* >(this_op)->add( *utterrandom, prob );
        else
          static_cast< eoProportionalOp<t_Object>* >(this_op)->add( *utterrandom, prob );
      }
      else if ( str.compare("Operators") == 0 )
      {
        std :: string special = _special + _base;
        eoGenOp<t_Object> *add_genop = make_recurrent_op( *child, _f,  special, _base);
        if ( is_sequential )
          static_cast< eoSequentialOp<t_Object>* >(this_op)->add( *add_genop, prob );
        else
          static_cast< eoProportionalOp<t_Object>* >(this_op)->add( *add_genop, prob );
      }
      else
        throw "Unknown operator in  Darwin<t_Object, t_Lamarck>  :: make_recurrent_op(...) ";
    }
    
    return this_op;
  }

  template< class t_Object, class t_Lamarck >
  eoMonOp<t_Object>* Darwin<t_Object, t_Lamarck> :: make_MonOp(const TiXmlElement &el,
                                                               std::ofstream &_f,
                                                               std::string &_special,
                                                               std::string &_base)
  {
    eoMonOp<t_Object>* this_op;
    bool is_sequential;
    if ( el.Attribute( "type" ) )
    {
       std::string str =  el.Attribute( "type" );
       is_sequential = true;
       if ( str.compare("and" ) == 0 ) // operators applied sequentially
       {
         double d;
         if ( el.Attribute("prob", &d ) )
           _f << "# " << _special << "Sequential: prob " << d << std::endl;
         else
           _f << "# " << _special << "Sequential" << std::endl;
         this_op = new  Sequential<t_Object>;
       }
       else
       {
         is_sequential = false;
         double d;
         if ( el.Attribute("prob", &d ) )
           _f << "# " << _special << "Proportional: prob " << d << std::endl;
         else
           _f << "# " << _special << "Proportional" << std::endl;
         this_op = new  Proportional<t_Object>;
       }
       eostates.storeFunctor(this_op);
    }
    else
      throw "Error while creating Operator ";

    const TiXmlElement *child = el.FirstChildElement();
    if (not child)
      throw "Error while creating Operator ";

    for ( ; child; child = child->NextSiblingElement() )
    {
      std::string str = child->Value();
      double prob = 0.0;
      
      // gets probability for applying child 
      if ( not child->Attribute("prob", &prob) )
        prob = 1.0;

      // then creates child
      if ( str.compare("Mutation" ) == 0 )
      {
        double d; 
        child->Attribute("value", &d);
        if ( d <= 0 and d > 1 )
          d = 1.0 / (double) lamarck->get_pb_size();
        Mutation<t_Object> *mutation = new Mutation<t_Object>( d );
        eostates.storeFunctor(mutation);
        _f << "# " << _special << _base << "Mutation: value=" << d 
           << " prob="<< prob << std::endl;
        if ( is_sequential )
          static_cast< Sequential<t_Object>* >(this_op)->add( mutation, prob );
        else
          static_cast< Proportional<t_Object>* >(this_op)->add( mutation, prob );
      }
      else if ( str.compare("Minimizer") == 0 )
      {
        _f << "# " << _special << _base << "Minimizer: ";
        eoMonOp<t_Object>* minop = Load_Minimizer( child, _f );
        _f << ", prob=" << prob << std::endl;
        if ( is_sequential )
          static_cast< Sequential<t_Object>* >(this_op)->add( minop, prob );
        else
          static_cast< Proportional<t_Object>* >(this_op)->add( minop, prob );
      }
      else if ( str.compare("UtterRandom") == 0 )
      {
        UtterRandom<t_Object>* utterrandom = new UtterRandom<t_Object>;
        eostates.storeFunctor(utterrandom);
        _f << "# " << _special << _base << "UtterRandom "
           << " prob "<< prob << std::endl;
        if ( is_sequential )
          static_cast< Sequential<t_Object>* >(this_op)->add( utterrandom, prob );
        else
          static_cast< Proportional<t_Object>* >(this_op)->add( utterrandom, prob );
      }
      else if ( str.compare("Operators") == 0 )
      {
        std :: string special = _special + _base;
        eoMonOp<t_Object> *add_genop = make_MonOp( *child, _f,  special, _base);
        if ( is_sequential )
          static_cast< Sequential<t_Object>* >(this_op)->add( add_genop, prob );
        else
          static_cast< Proportional<t_Object>* >(this_op)->add( add_genop, prob );
      }
      else
        throw "Unknown operator in  Darwin<t_Object, t_Lamarck>  :: make_MonOp(...) ";
    }
    
    return this_op;
  }

  template<class t_Object, class t_Lamarck>
  eoBreed<t_Object>* Darwin<t_Object, t_Lamarck> :: make_breeder()
  {
    eoGeneralBreeder<t_Object> *breed;
    eoGenOp<t_Object> *op;
    eoSelectOne<t_Object> *select;

    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app );

    select = new eoDetTournamentSelect<t_Object>(tournament_size);

    TiXmlDocument doc( filename.c_str() );
    TiXmlHandle docHandle( &doc );
    
    if  ( !doc.LoadFile() )
    {
      std::cout << doc.ErrorDesc() << std::endl; 
      throw "Could not load input file in  Darwin<t_Object, t_Lamarck>  :: make_GenOp";
    }
    TiXmlElement *child = docHandle.FirstChild("LaDa")
                                   .FirstChild("GA")
                                   .FirstChild("Operators").Element();
    if ( not child )
      throw "Could not find Operators in input file ";

    op = make_GenOp( child, xmgrace_file );
    xmgrace_file << "# Breeding Operator begin " << std::endl;
    breed = new eoGeneralBreeder<t_Object>(*select,*op , replacement_rate);
    xmgrace_file << "# Breeding Operator end " << std::endl;
    eostates.storeFunctor(breed);

    xmgrace_file.flush();
    xmgrace_file.close();
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
    Evaluation<t_Object, Darwin< t_Object, t_Lamarck> > *evaluation
       = new Evaluation<t_Object, Darwin< t_Object, t_Lamarck> >(*this);
    popEval = new EvaluatePop<t_Object>(*evaluation);
    
    continuator = make_checkpoint();
    breed = make_breeder();
    replace = make_replacement();
    make_extra_algo();
    eostates.storeFunctor(evaluation);
    eostates.storeFunctor(popEval);
  }

  template< class t_Object, class t_Lamarck >
  eoCheckPoint<t_Object>* Darwin<t_Object, t_Lamarck> :: make_checkpoint()
  {
    // continuator
    eoGenContinue<t_Object> *gen_continue = new eoGenContinue<t_Object>(max_generations);
    eostates.storeFunctor( gen_continue );

   
    // gen_continue
    eoCheckPoint<t_Object> *check_point = new eoCheckPoint<t_Object>(*gen_continue);
    eostates.storeFunctor( check_point );

    // our very own updater wrapper to print stuff
    Monitor< Darwin<t_Object, t_Lamarck> > *updater = new Monitor< Darwin<t_Object, t_Lamarck> >(this);
    eostates.storeFunctor(updater);
    check_point->add(*updater);
    
    // gen_continue -- should be last updater to be added
    nb_generations = new eoIncrementorParam<unsigned>("Gen.");
    eostates.storeFunctor(nb_generations);
    check_point->add(*nb_generations);
    
    return check_point;
  }

  template<class t_Object, class t_Lamarck>
  void Darwin<t_Object, t_Lamarck> :: populate ()
  {
    Generator generator;
    t_Object indiv;
    indiv.resize( lamarck->get_pb_size() );
    population.clear();
    population.reserve(pop_size);
    for( unsigned i = 0; i < pop_size; ++i )
    {
      typename t_Object :: iterator i_var = indiv.begin();
      typename t_Object :: iterator i_end = indiv.end();
      for ( ; i_var != i_end; ++i_var)
        *i_var = generator();
      population.push_back(indiv);
    }
  }

  template< class t_Object, class t_Lamarck >
  void Darwin <t_Object, t_Lamarck> :: write_xmgrace_header()
  {
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    xmgrace_file << "# population size: " << pop_size << std::endl;
    xmgrace_file << "# replacement rate: " << replacement_rate << std::endl;
    xmgrace_file << "# max generations: " << max_generations << std::endl;
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
  void Darwin <t_Object, t_Lamarck> :: print_xmgrace()
  {
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 
    bool print_ch = not population.begin()->is_baseline_valid();
    if ( not print_ch )
      xmgrace_file<< " # ? iteration:" << nb_generations->value(); 
    else 
      xmgrace_file << " # iteration:" << nb_generations->value(); 
    lamarck->print_xmgrace( xmgrace_file,  print_ch );
    if ( print_ch )
      population.begin()->validate_baseline();
    xmgrace_file.flush();
    xmgrace_file.close();
  }

  template< class t_Object, class t_Lamarck >
  void Darwin <t_Object, t_Lamarck> :: make_extra_algo()
  {
    TiXmlDocument doc( filename.c_str() );
    TiXmlHandle docHandle( &doc );
    
    if  ( !doc.LoadFile() )
    {
      std::cout << doc.ErrorDesc() << std::endl; 
      throw "Could not load input file in  Darwin<t_Object, t_Lamarck>  :: make_GenOp";
    }
    TiXmlElement *child = docHandle.FirstChild("LaDa")
                                   .FirstChild("GA")
                                   .FirstChild("PopAlgo")
                                   .FirstChild("Operators").Element();
    if ( not child )
      throw "Could not find PopAlgo Operators in input file ";

    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("PopAlgo").Element();
    if ( not child )
      throw "Could not find PopAlgo in input file ";
    
    if (  not child )
      extra_popalgo = NULL;
    
    // get parameters
    {
      minimize_best = 0.1;
      double d=0;
      if ( child->Attribute( "rate", &d ) )
        if ( d > 0 and d <= 1 )
          minimize_best = d;
      int u = 0;
      if ( child->Attribute( "every", &u ) )
        minimize_best_every = ( u > 0 and abs(u) <= max_generations ) ? (unsigned) u  : 0 ;
      std::string str;
    }

    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app );
    xmgrace_file << "# PopAlgorithm: rate="<< minimize_best 
                 << " every=" << minimize_best_every << std::endl;
    child = docHandle.FirstChild("LaDa")
                     .FirstChild("GA")
                     .FirstChild("PopAlgo")
                     .FirstChild("Operators").Element();
    std::string str = "    ";
    std::string base = "  ";
    eoMonOp<t_Object> *op = make_MonOp(*child, xmgrace_file, str, base);

    extra_popalgo = new Extra_PopAlgo< t_Object, Darwin<t_Object, t_Lamarck> > 
                                     ( *op, *this, minimize_best, minimize_best_every );
    eostates.storeFunctor( extra_popalgo );
    xmgrace_file.flush();
    xmgrace_file.close();
  }
} // namespace LaDa


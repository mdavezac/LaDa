
#include "motu.h"
#include "evaluation.h"
#include "operators.h"
#include "generator.h"
#include "checkpoint.h"

namespace LaDa 
{

  const double MotU :: ZERO_TOLERANCE = 1e-6;
  const int MotU :: DARWINISTIC          = 0;
  const int MotU :: LAMARCKIAN           = 1;
  const int MotU :: MULTISTART           = 2;
  const int MotU :: LAMARCKIAN_EVOLUTION = 3;

  MotU :: GA_Params :: GA_Params()
  {
    crossover_vs_mutation = 0.85;
    crossover_probability = 0.5;
    mutation_probability = 0.0;
    sequential_op = false;
    tournament_size = 2;
    replacement_rate = 0.1;
    max_generations = 200;
    pop_size = 100;
    method = DARWINISTIC;
    utter_random = false;
  }
  
  MotU :: MotU(const std::string &_filename)
  {
    MotU::MotU();
    TiXmlDocument doc( filename.c_str() );
    
    if  ( !doc.LoadFile() )
    {
      std::cout << doc.ErrorDesc() << std::endl; 
      return;
    }

    TiXmlHandle docHandle( &doc );
    if ( Load ( docHandle ) )
      read_CH();
    else
      std::cerr << " Error while loading Lamarck parameters from "
                << filename
                << std::endl;
  }

  bool MotU :: read_CH()
  {
    TiXmlDocument doc( xml_filename.c_str() );
    
    if  ( !doc.LoadFile() ) // no convex hull to read...
      return false;
  
    TiXmlHandle handle( &doc );
    TiXmlElement *child = handle.FirstChild( "LaDa" ).FirstChild( "ConvexHull" ).Element();
    if ( not child )
      return false;
    
    return convex_hull.Load(child, *axes);
  }


  // loads input from filename file
  bool MotU :: Load( TiXmlHandle &handle )
  {
    TiXmlElement *child;
    rVector3d vec;

    // clusters, lattice, harmonics ....
    Functional_Builder :: Load (handle );
    
  
    xml_filename = "convex_hull.xml";
    child = handle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
    if ( child and child->Attribute("xml") )
      xml_filename = child->Attribute("xml");
      
    xmgrace_filename = "convex_hull.agr";
    child = handle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
    if ( child and child->Attribute("xmgrace") )
      xmgrace_filename = child->Attribute("xmgrace");

    // reads structure from input
    child = handle.FirstChild( "LaDa" ).FirstChild( "Structure" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Structure in " << filename << std::endl;
      return false;
    }
    structure.Load(child, *axes);
    ga_params.mutation_probability = 1.0 / ((double) structure.atoms.size());

    // finds GA parameters
    child = handle.FirstChild( "LaDa" ).FirstChild( "GA" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find input file " << filename << std::endl;
      return false;
    }
    ga_params.Load(child);

    // restart content of xmgrace file 
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    #if defined(WANG_CONSTRAINTS)
      xmgrace_file << "# Wang Constraints " << std::endl; 
    #elif defined(LINEAR_SOLVE)
      xmgrace_file << "# Linear Solver " << std::endl; 
    #else
      xmgrace_file << "# S^2 - 1 = 0 Constraints " << std::endl; 
    #endif
    #ifdef ONE_POINT
      xmgrace_file << "# One point convex hull " << std::endl; 
    #else
      xmgrace_file << "# N-point convex hull " << std::endl; 
    #endif
    switch( ga_params.method )
    {
      case LAMARCKIAN : xmgrace_file << "# Lamarckian GA" << std::endl; break;
      case MULTISTART : xmgrace_file << "# Multistart" << std::endl; break;
      default:
      case DARWINISTIC : xmgrace_file << "# Darwinistic GA" << std::endl; break;
    }
    xmgrace_file.flush();
    xmgrace_file.close();

    return true;
  }


  // GA specific paremeters
  bool MotU :: GA_Params :: Load( TiXmlElement *element)
  {
    TiXmlElement *child;
    
    // lamarckian method 
    child = element->FirstChildElement( "method" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("type", &d) )
        switch ( d )
        {
          case LAMARCKIAN : method = LAMARCKIAN; break;
          case LAMARCKIAN_EVOLUTION : method = LAMARCKIAN_EVOLUTION; break;
          default:
          case DARWINISTIC : method = DARWINISTIC; break;
        }
    }

    // crossover_vs_mutation
    child = element->FirstChildElement( "crossover_vs_mutation" );
    if ( child )
    {
      double d = 0;
      if ( child->Attribute("value", &d) )
        crossover_vs_mutation = d;
    }

    // crossover_probability
    child = element->FirstChildElement( "crossover_probability" );
    if ( child )
    {
      double d = 0;
      if ( child->Attribute("value", &d) )
        crossover_probability = d;
    }

    // mutation_probability
    child = element->FirstChildElement( "mutation_probability" );
    if ( child )
    {
      double d = 0;
      if ( child->Attribute("value", &d) )
        mutation_probability = d;
    }

    // crossover and mutation applied sequentialy ?
    if ( element->FirstChildElement( "sequential" ) )
      sequential_op = true;
    
    // crossover and mutation applied sequentialy ?
    if ( element->FirstChildElement( "utter random" ) )
      utter_random = true;

    // tournament size when selecting parents
    child = element->FirstChildElement( "tournament size" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("value", &d) )
        if ( d > 2 )
          tournament_size = d;
    }

    // tournament size when selecting parents
    child = element->FirstChildElement( "replacement rate" );
    if ( child )
    {
      double d = 0;
      if ( child->Attribute("value", &d) )
        if ( d <= 1.0 and d > 0.0 )
          replacement_rate = d;
    }

    // population size
    child = element->FirstChildElement( "population size" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("value", &d) )
        if ( d > 0 )
          pop_size = d;
    }
    
    // max generations
    child = element->FirstChildElement( "max generations" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("value", &d) )
        if ( d > 0 )
          max_generations = d;
    }

    // population size
    child = element->FirstChildElement( "population size" );
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("value", &d) )
        if ( d > 0 )
          pop_size = d;
    }


    return true;
  }

  void MotU :: populate ()
  {
    Generator generator;
    eoInitFixedLength<t_individual> init( structure.atoms.size(), generator);
    population.append( ga_params.pop_size, init);
  }

  eoCheckPoint<MotU::t_individual>* MotU :: make_checkpoint()
  {
    // continuator
    eoGenContinue<t_individual> *gen_continue = new eoGenContinue<t_individual>(ga_params.max_generations);
    eostates.storeFunctor( gen_continue );

    // checkpoints
      // gen_continue
    eoCheckPoint<t_individual> *check_point = new eoCheckPoint<t_individual>(*gen_continue);
    eostates.storeFunctor( check_point );

      // gen_continue
    eoIncrementorParam<unsigned> *nb_generations = new eoIncrementorParam<unsigned>("Gen.");
    eostates.storeFunctor(nb_generations);
    check_point->add(*nb_generations);

      // our very own updater wrapper to print stuff
    Monitor<MotU> *updater = new Monitor<MotU>(this);
    eostates.storeFunctor(updater);
    check_point->add(*updater);
    
    return check_point;
  }
  
  void MotU :: print_xmgrace()
  {
    if( not population.begin()->is_baseline_valid() )
    {
       std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 
       xmgrace_file << " # iteration nb " << nb_generations->value() 
                    << " " << EvalCounter << std::endl;
       xmgrace_file << " # evaluation calls " 
                    << VA_CE::Polynome::nb_eval << " "
                    << VA_CE::Polynome::nb_eval_grad << " "
                    << VA_CE::Polynome::nb_eval_with_grad << std::endl;
       convex_hull.print_out(xmgrace_file, CONVEX_HULL::PRINT_XMGRACE);
       xmgrace_file.flush();
       xmgrace_file.close();
       population.begin()->validate_baseline();
    }
  }

  void MotU :: print_xml()
  {
    TiXmlDocument doc;
    TiXmlElement *LaDa_node = new TiXmlElement("LaDa");
    TiXmlElement *convex_hull_node = new TiXmlElement("ConvexHull");

    doc.SetTabSize(1);
    doc.LinkEndChild( new TiXmlDeclaration("1.0", "", "") );
    
    doc.LinkEndChild( LaDa_node );
    LaDa_node->LinkEndChild( convex_hull_node ); 
    convex_hull.print_xml( convex_hull_node, *axes );

    doc.SaveFile(xml_filename.c_str());
  }

  // returns polynomial value only
  double MotU :: evaluate( t_individual &individual )
  {
    ++EvalCounter;

    double result;
    functional.set_variables( individual.get_variables() );

    switch( ga_params.method )
    {
      case LAMARCKIAN_EVOLUTION:
        if ( nb_generations->value() == 0 )
          break;
      case LAMARCKIAN:
        minimizer.minimize();
        functional.set_to_closest_constraints();
      default:  break;
    }

    result = functional.evaluate();

    structure.set_atom_types( *individual.get_variables() );
    if ( convex_hull.add_structure(result, structure) )
      population.begin()->invalidate_baseline();
  }

  eoGenOp<MotU::t_individual>* MotU :: make_GenOp()
  {
    // completely randomizes offsprings
    if ( ga_params.utter_random )
    {
      UtterRandom<t_individual> *random = new UtterRandom<t_individual>;
      eoMonGenOp<t_individual> *op = new  eoMonGenOp<t_individual>( *random );
      eostates.storeFunctor( random );
      eostates.storeFunctor( op );
      return op;
    }

    // creates crossover and mutation operators
    Crossover<t_individual> *crossover = new Crossover<t_individual>( ga_params.crossover_probability );
    Mutation<t_individual> *mutation = new Mutation<t_individual>( ga_params.mutation_probability );
    eostates.storeFunctor(crossover);
    eostates.storeFunctor(mutation);

    // which are applied either sequentially
    if ( ga_params.sequential_op )
    {
      eoSequentialOp<t_individual> *op = new  eoSequentialOp<t_individual>;
      eostates.storeFunctor(op);
      op->add(*crossover, 1.0);
      op->add(*mutation, 1.0);
      return op;
    }

    // or propotionally
    eoProportionalOp<t_individual> *op = new  eoProportionalOp<t_individual>;
    eostates.storeFunctor(op);
    op->add(*crossover, ga_params.crossover_vs_mutation);
    op->add(*mutation, 1.0 - ga_params.crossover_vs_mutation );

    return op;
  }

  eoBreed<MotU::t_individual>* MotU :: make_breeder()
  {
    eoSelectOne<t_individual> *select;
    select = new eoDetTournamentSelect<t_individual>(ga_params.tournament_size);

    eoGeneralBreeder<t_individual> *breed;
    breed = new eoGeneralBreeder<t_individual>(*select, *make_GenOp(), ga_params.replacement_rate);
    eostates.storeFunctor(breed);
  }
  
  eoReplacement<MotU::t_individual>* MotU :: make_replacement()
  {
    eoTruncate<t_individual>* truncate = new  eoTruncate<t_individual>;
    eoMerge<t_individual>* merge = new  eoPlus<t_individual>;
    eoReduceMerge<t_individual>* reducemerge = new eoReduceMerge<t_individual>( *truncate, *merge );
    eostates.storeFunctor(truncate);
    eostates.storeFunctor(merge);
    eostates.storeFunctor(reducemerge);
  }

  void MotU :: make_algo()
  {
    Evaluation<t_individual, MotU> *evaluation = new Evaluation<t_individual, MotU>(this);
    EvaluatePop<t_individual> *evalpop = new EvaluatePop<t_individual>(*evaluation);
    algorithm = new eoEasyEA<t_individual>( *make_checkpoint(), 
                                            *evalpop, 
                                            *make_breeder(), 
                                            *make_replacement() );
    eostates.storeFunctor(evaluation);
    eostates.storeFunctor(evalpop);
    eostates.storeFunctor(algorithm);
  }

  void MotU :: run()
  {
    // generate functional
    add_equivalent_clusters();
    generate_functional( structure, functional );

    // fitness pointers  -- note that functiona->variables will be set
    // to Individual::variables at each evaluation
    fitness.set_baseline( &convex_hull );
    fitness.set_quantity( &functional ); 
    functional.get_Obj1()->destroy_variables(); 

    // minimizer
    minimizer.set_object( fitness );

    try 
    {
      // create algorithm
      make_algo();
      
      // create starting population
      population.clear();
      populate();
      
      // finally, runs algorithm
      algorithm->operator()( population );
    }
    catch (std::exception &e)
    {
       std::cerr << "Error while running LaDa::Motu::run() " << std::endl
                 << e.what() << std::endl; 
    }

  }
} // namespace LaDa

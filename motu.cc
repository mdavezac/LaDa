#include "motu.h"
#include "evaluation.h"
#include "operators.h"
#include "generator.h"
#include "checkpoint.h"

#undef min // idiots
#undef max
#include <eo/utils/eoHowMany.h>
#undef min // idiots
#undef max
#include <limits.h>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>

namespace LaDa 
{

  const unsigned MotU :: DARWIN  = 0;
  const unsigned MotU :: LAMARCK = 1;
  const unsigned MotU :: DEBUG   = 2;
  const unsigned MotU :: GA :: NO_MINIMIZER       = 0;
  const unsigned MotU :: GA :: WANG_MINIMIZER     = 1;
  const unsigned MotU :: GA :: PHYSICAL_MINIMIZER = 2;
  const unsigned MotU :: GA :: LINEAR_MINIMIZER   = 3;
  const unsigned MotU :: GA :: SA_MINIMIZER       = 4;

  MotU :: GA :: GA()
  {
    crossover_vs_mutation = 0.85;
    crossover_probability = 0.5;
    mutation_probability = 0.0;
    sequential_op = false;
    tournament_size = 2;
    replacement_rate = 0.1;
    max_generations = 200;
    pop_size = 100;
    method = DARWIN;
    utter_random = false;
    evolve_from_start = false; 
    multistart = false; 
    minimizer = NO_MINIMIZER;
    max_eval_calls = 0;
    max_grad_calls = 0;
  }
  
  bool MotU :: Load(const std::string &_filename) 
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
    if ( not Load ( docHandle ) )
    {
      std::cerr << "error while  loading Lamarck parameters from " << filename << std::endl
                << " tinyxml: " << doc.ErrorDesc() << std::endl;
      return false;
    }
    read_CH();
    return true;
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
    if ( not Functional_Builder :: Load (handle ) )
      return false;
     
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

    // sets hull type: one or many points
    convex_hull.is_flat = false;
    child = handle.FirstChild( "LaDa" ).FirstChild( "OnePointHull" ).Element();
    if ( child )
      convex_hull.is_flat = true;
    
    // finds GA parameters
    child = handle.FirstChild( "LaDa" ).FirstChild( "GA" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find GA in input file " << filename << std::endl;
      return false;
    }
    if ( not  ga_params.Load(child) )
      return false;

    // restart content of xmgrace file 
    std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::trunc ); 
    write_xmgrace_header( xmgrace_file );
    xmgrace_file.flush();
    xmgrace_file.close();

    return true;
  }

  void MotU :: write_xmgrace_header(std::ofstream  &_f)
  {
    switch( ga_params.method )
    {
      case DEBUG: _f << "# job: debug" << std::endl; break;
      case LAMARCK: _f << "# job: Lamarckian GA" << std::endl; break;
      default:
      case DARWIN: _f << "# job: Darwinistic GA" << std::endl; break;
    }
    switch( ga_params.minimizer )
    {
      case GA::WANG_MINIMIZER: 
        _f << "# Wang Constraints " << std::endl; 
        break;
      case GA::PHYSICAL_MINIMIZER: 
        _f << "# Physical constraints " << std::endl; 
        break;
      case GA::LINEAR_MINIMIZER: 
        _f << "# Linear minimizer " << std::endl; 
        break;
      case GA::SA_MINIMIZER: 
        _f << "# T=0 Simulated Annealing " << std::endl; 
        break;
      default:
        _f << "# No minimizer " << std::endl; 
    }
    if ( ga_params.multistart )
        _f << "# multistart " << std::endl; 
    if ( ga_params.evolve_from_start )
        _f << "# evolve starting population " << std::endl; 
    if (convex_hull.is_flat)
      _f << "# One point convex hull " << std::endl; 
    else
      _f << "# N-point convex hull " << std::endl; 

    _f << "# population size: " << ga_params.pop_size << std::endl;
    _f << "# replacement rate: " << ga_params.replacement_rate << std::endl;
    _f << "# max generations: " << ga_params.max_generations << std::endl;
    if ( ga_params.utter_random )
      _f << "# utter random operator: " << std::endl;
    else
    {
      _f << "# sequential operators: " << ga_params.sequential_op << std::endl;
      _f << "# crossover vs mutation: " << ga_params.crossover_vs_mutation << std::endl;
      _f << "# crossover prob: " << ga_params.crossover_probability << std::endl;
      _f << "# mutation prob: " << ga_params.mutation_probability << std::endl;
    }
  }


  // GA specific paremeters
  bool MotU :: GA :: Load( TiXmlElement *element)
  {
    TiXmlElement *child, *parent;
    
    // finds if in <GA> ... </GA> block 
    {
      std::string str = element->Value();
      parent = element;
      if ( str.compare("GA" ) != 0 )
        parent = element->FirstChildElement("GA");
    }
    

    // Checks Operators
    child = parent->FirstChildElement( "Operators" );
    if ( child )
    {
      if ( child->Attribute( "type" ) )
      {
        std::string str =  child->Attribute( "type" );
        if ( str.compare("and" ) == 0 ) // operators applied sequentially
          sequential_op = true;
      }
      
      // Crossover vs Mutation
      double d = 0;
      if ( child->Attribute("prob", &d) and d > 0 and d < 1)
        crossover_vs_mutation = d;
      
      // if tag present, then applies utter random
      if ( child->FirstChildElement( "utter random" ) )
        utter_random = true;

      // gets Mutation prob
      TiXmlElement *grandchild = child->FirstChildElement( "Mutation" ); d=0;
      if (     grandchild 
           and grandchild->Attribute("prob", &d) and d > 0 and d < 1)
        mutation_probability = d;
      
      // gets Crossover prob
      grandchild = child->FirstChildElement( "Crossover" ); d=0;
      if (     grandchild 
           and grandchild->Attribute("prob", &d) and d > 0 and d < 1)
        crossover_probability = d;
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
      //
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
        if ( str.find("lamarck") != std::string::npos )
        {
          method = LAMARCK;
          if ( str.find("from start") != std::string::npos )
            evolve_from_start = true;
          if ( minimizer == NO_MINIMIZER )
            minimizer = LINEAR_MINIMIZER;
        }
        if ( str.find("debug") != std::string::npos )
          method = DEBUG;
        if ( str.find("multistart") != std::string::npos )
        {
          multistart = true;
          eoHowMany nb(replacement_rate);
          max_generations = pop_size + nb(pop_size) * max_generations;
          pop_size = 1;
          utter_random = true;
          replacement_rate = 1.0;
          evolve_from_start = true;
        }
      } // if attribute "method" exists
    }   

      
    if ( method != DARWIN )  // no optimizer - keep default
    {
      child = parent->FirstChildElement( "Minimizer" );
      if ( child->Attribute("type") )
      {
        std::string str = child->Attribute("type");
        if ( str.compare("wang" ) == 0 ) // Wang
          minimizer = WANG_MINIMIZER;
        else if ( str.compare("physical" ) == 0 ) // Wang
          minimizer = PHYSICAL_MINIMIZER;
        else if ( str.compare("linear" ) == 0 ) // Wang
          minimizer = LINEAR_MINIMIZER;
        else if ( str.compare("SA" ) == 0 ) // Wang
          minimizer = SA_MINIMIZER;
      }
    }
    // checks for minimizer setup
    if (minimizer == LINEAR_MINIMIZER or minimizer == SA_MINIMIZER )
    {
      int d=0;
      if( child->Attribute("maxeval", &d) )
        max_eval_calls = ( d <= 0 ) ? UINT_MAX : (unsigned) d;
      if( child->Attribute("maxgrad", &d) )
        max_grad_calls = ( d <= 0 ) ? max_eval_calls : (unsigned) d;
    }

    return true;
  }

  void MotU :: populate ()
  {
    Generator generator;
    t_individual indiv;
    indiv.resize(structure.atoms.size());
    population.clear();
    population.reserve(ga_params.pop_size);
    for( unsigned i = 0; i < ga_params.pop_size; ++i )
    {
      Individual<> :: CONTAINER_ITERATOR i_var = indiv.begin();
      Individual<> :: CONTAINER_ITERATOR i_end = indiv.end();
      for ( ; i_var != i_end; ++i_var)
        *i_var = generator();
      population.push_back(indiv);
    }
  }

  eoCheckPoint<MotU::t_individual>* MotU :: make_checkpoint()
  {
    // continuator
    eoGenContinue<t_individual> *gen_continue = new eoGenContinue<t_individual>(ga_params.max_generations);
    ga_params.eostates.storeFunctor( gen_continue );

    // checkpoints
      // gen_continue
    eoCheckPoint<t_individual> *check_point = new eoCheckPoint<t_individual>(*gen_continue);
    ga_params.eostates.storeFunctor( check_point );

      // our very own updater wrapper to print stuff
    Monitor<MotU> *updater = new Monitor<MotU>(this);
    ga_params.eostates.storeFunctor(updater);
    check_point->add(*updater);
    
      // gen_continue -- should be last updater to be added
    ga_params.nb_generations = new eoIncrementorParam<unsigned>("Gen.");
    ga_params.eostates.storeFunctor(ga_params.nb_generations);
    check_point->add(*ga_params.nb_generations);
    
    return check_point;
  }
  
  void MotU :: print_xmgrace()
  {
    if( not population.begin()->is_baseline_valid() )
    {
       std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 
       xmgrace_file << " # iteration:" << ga_params.nb_generations->value() 
                    << "   evaluation calls: " << EvalCounter << std::endl;
       xmgrace_file << " # polynomial calls: " 
                    << VA_CE::Polynome::nb_eval << " "
                    << VA_CE::Polynome::nb_eval_grad << " "
                    << VA_CE::Polynome::nb_eval_with_grad << std::endl;
       convex_hull.print_out(xmgrace_file, VA_CE::Convex_Hull::PRINT_XMGRACE);
       xmgrace_file.flush();
       xmgrace_file.close();
       population.begin()->validate_baseline();
    }
//   else
//   {
//      std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 
//      xmgrace_file << " # iteration:" << ga_params.nb_generations->value() 
//                   << " no change " << std::endl;
//      xmgrace_file.flush();
//      xmgrace_file.close();
//   }
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
    fitness.set_variables( individual.get_variables() );
 
    if ( ga_params.evolve_from_start
         or ga_params.nb_generations->value() )
      minimizer->minimize();
 
    result = functional.evaluate();
 
    structure.set_atom_types( *individual.get_variables() );
    if ( convex_hull.add_structure(result, structure) )
      population.begin()->invalidate_baseline();
    
    return result;
  }

  eoGenOp<MotU::t_individual>* MotU :: make_GenOp()
  {
    // completely randomizes offsprings
    if ( ga_params.utter_random )
    {
      UtterRandom<t_individual> *random = new UtterRandom<t_individual>;
      eoMonGenOp<t_individual> *op = new  eoMonGenOp<t_individual>( *random );
      ga_params.eostates.storeFunctor( random );
      ga_params.eostates.storeFunctor( op );
      return op;
    }

    // creates crossover and mutation operators
    Crossover<t_individual> *crossover = new Crossover<t_individual>( ga_params.crossover_probability );
    Mutation<t_individual> *mutation = new Mutation<t_individual>( ga_params.mutation_probability );
    ga_params.eostates.storeFunctor(crossover);
    ga_params.eostates.storeFunctor(mutation);

    // which are applied either sequentially
    if ( ga_params.sequential_op )
    {
      eoSequentialOp<t_individual> *op = new  eoSequentialOp<t_individual>;
      ga_params.eostates.storeFunctor(op);
      op->add(*crossover, 1.0);
      op->add(*mutation, 1.0);
      return op;
    }

    // or propotionally
    eoProportionalOp<t_individual> *op = new  eoProportionalOp<t_individual>;
    ga_params.eostates.storeFunctor(op);
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
    ga_params.eostates.storeFunctor(breed);

    return breed;
  }
  
  eoReplacement<MotU::t_individual>* MotU :: make_replacement()
  {
    eoTruncate<t_individual>* truncate = new  eoTruncate<t_individual>;
    eoMerge<t_individual>* merge = new  eoPlus<t_individual>;
    eoReduceMerge<t_individual>* reducemerge = new eoReduceMerge<t_individual>( *truncate, *merge );
    ga_params.eostates.storeFunctor(truncate);
    ga_params.eostates.storeFunctor(merge);
    ga_params.eostates.storeFunctor(reducemerge);

    return reducemerge;
  }

  void MotU :: make_algo()
  {
    Evaluation<t_individual, MotU> *evaluation = new Evaluation<t_individual, MotU>(this);
    EvaluatePop<t_individual> *evalpop = new EvaluatePop<t_individual>(*evaluation);
    ga_params.algorithm = new eoEasyEA<t_individual>( *make_checkpoint(), 
                                            *evalpop, 
                                            *make_breeder(), 
                                            *make_replacement() );
    ga_params.eostates.storeFunctor(evaluation);
    ga_params.eostates.storeFunctor(evalpop);
    ga_params.eostates.storeFunctor(ga_params.algorithm);
  }

  // creates and sets functionals and minimizers,
  // then branches off to different jobs
  void MotU :: run()
  {
     // completes cluster list with equivalent clusters
     add_equivalent_clusters();
    
     // adds endpoints to convex hull
     init_convex_hull();
     
    // debug case: iterates over structures
    if( ga_params.method == DEBUG)
    {
      run_debug();
      return;
    }

    // generate functional
    generate_functional( structure, &functional );

    // fitness pointers  -- note that functiona->variables will be set
    // to Individual::variables at each evaluation
    fitness.set_baseline( &convex_hull );
    fitness.set_quantity( &functional ); 
    functional.destroy_variables(); 

    try 
    {
      // minimizer
      switch( ga_params.minimizer )
      {
        case GA::WANG_MINIMIZER: 
          minimizer = new opt::Minimize_Wang<FITNESS>( &fitness);
          break;
        case GA::PHYSICAL_MINIMIZER: 
          minimizer = new opt::Minimize_Ssquared<FITNESS>( &fitness );
          break;
        case GA::SA_MINIMIZER: 
          minimizer = new opt::Minimize_Linear<FITNESS>( &fitness );
          static_cast< opt::Minimize_Linear<FITNESS>* >(minimizer)->simulated_annealing = true;
          static_cast< opt::Minimize_Linear<FITNESS>* >(minimizer)->max_eval_calls = ga_params.max_eval_calls;
          break;
        case GA::LINEAR_MINIMIZER: 
          minimizer = new opt::Minimize_Linear<FITNESS>( &fitness );
          static_cast< opt::Minimize_Linear<FITNESS>* >(minimizer)->simulated_annealing = false;
          static_cast< opt::Minimize_Linear<FITNESS>* >(minimizer)->max_eval_calls = ga_params.max_eval_calls;
          static_cast< opt::Minimize_Linear<FITNESS>* >(minimizer)->max_grad_calls = ga_params.max_grad_calls;
          break;
        default:
          minimizer = new opt::Minimize_Base<FITNESS>( &fitness ); // just a dummy, doesn't minimize
          break;
      }


      // create algorithm
      make_algo();
      
      // create starting population
      population.clear();
      populate();
      
      // finally, runs algorithm
      ga_params.algorithm->operator()( population );
    }
    catch (std::exception &e)
    {
       std::cerr << "Error while running LaDa::Motu::run() " << std::endl
                 << e.what() << std::endl; 
    }

    delete minimizer;
    minimizer = NULL;
  }

  // Debug job: _ finds Structures tag in input.xml
  //            _ evaluate each structure cum atoms in turn
  void MotU :: run_debug()
  {
    // first creates doc handle
    TiXmlDocument doc( filename.c_str() );
    
    if  ( !doc.LoadFile() )
    {
      std::cout << doc.ErrorDesc() << std::endl; 
      return;
    }

    TiXmlHandle docHandle( &doc );
    
    // then finds <Structures>
    TiXmlElement *child = docHandle.FirstChild("LaDa")
                                   .FirstChild("Structures")
                                   .FirstChild("Structure").Element();
                              
    if ( not child )
    {
      std::cerr << "Could not find structures in " << filename << std::endl;
      return;
    }

    // now iterates over structures
    for(; child; child = child->NextSiblingElement("Structure"))
    {
      structure.Load(child, *axes);
      generate_functional( structure, &functional );
      {
        std::vector< Ising_CE::Atom > :: iterator i_atom = structure.atoms.begin();
        std::vector< Ising_CE::Atom > :: iterator i_last  = structure.atoms.end();
        VA_CE::Polynome :: CONTAINER_ITERATOR i_var  = functional.begin();
        
        for(; i_atom != i_last; ++i_atom, ++i_var )
          *i_var = i_atom->type; 
      }
      std::cout << functional.get_concentration() << " "
                << functional.get_Obj1()->evaluate() << "   "
                << functional.get_Obj2()->evaluate() << "   "
                << structure.energy << " - "
                << functional.evaluate() << " = "
                << structure.energy - functional.evaluate() << std::endl;
      functional.destroy_variables(); // variables must be destroyed explicitely!!
      delete functional.get_Obj1(); 
      delete functional.get_Obj2(); 
    }

  }

   // adds endpoints to convex hull -- cluster list should be complete
   // at this point!!
   // single site CEs only
   void MotU :: init_convex_hull()
   {
     Ising_CE::Structure struc; 
     FUNCTIONAL func;

     // creates structure...
     struc.cell = lattice->cell;
     struc.atoms.push_back( Ising_CE::Atom( lattice->atom_pos(0), -1.0 ) );

     // and functional - note that there is no CS 
     // (in fact CS may segfault here )
     generate_functional(struc, &func);
     

     // adds first endpoint
     {
       *(func.begin()) = -1.0;
       VA_CE::Breaking_Point bp(func.get_Obj1()->evaluate(), struc);
       convex_hull.force_add( bp );
     }

     // adds second endpoint
     {
       *(func.begin()) = 1.0;
       struc.atoms.clear();
       struc.atoms.push_back( Ising_CE::Atom( lattice->atom_pos(0), 1.0 ) );
       VA_CE::Breaking_Point bp(func.get_Obj1()->evaluate(), struc);
       convex_hull.force_add( bp );
     }

     // clean-up
     func.destroy_variables();
     delete func.get_Obj1();
     delete func.get_Obj2();
   }
} // namespace LaDa

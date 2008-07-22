#ifdef _IAGA_LAMARCK_
  #include <fstream>
  #include <iostream>
  #include <sstream>
  #include <iomanip>
  #include <algorithm>

  #include "lamarck_interface.h"

  #include <opt/ndim_iterator.h>
  #include <analysis/analyze_code.h>
  #include <analysis/call_log.h>

  #include <lamarck/structure.h>

  #undef min
  #undef max
  #undef MIN
  #undef MAX
  #include <atat/parse.h>
  #include <atat/stringo.h>
  #include <atat/mrefine.h>
  #undef min
  #undef max
  #undef MIN
  #undef MAX

  
  // *************************
  // Lamarck class subroutines 
  // *************************

  extern "C" { LatticeSpec *GetGlobalLattice(); }

  const double Lamarck :: ZERO_TOLERANCE = 1e-6;
  const int Lamarck :: XMGRACE_FORMAT       = 0;
  const int Lamarck :: XML_FORMAT           = 1;
  const int Lamarck :: DARWINISTIC          = 0;
  const int Lamarck :: LAMARCKIAN           = 1;
  const int Lamarck :: MULTISTART           = 2;
  const int Lamarck :: LAMARCKIAN_EVOLUTION = 3;

  Lamarck :: ~Lamarck()
  {
    functional.get_Obj1()->destroy_variables();
  }

  Lamarck :: Lamarck(const std::string &_filename)
  {
    convex_hull_has_changed = false;
    nb_iterations = 0;
    GA_style = DARWINISTIC;
    CalculatingOldPop = false;


    filename = _filename;
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

  double Lamarck :: optimize(AtomicConfig *pAtoms)
  {
    VA_CE::Polynome :: CONTAINER_ITERATOR i_begin  = functional.begin();
    VA_CE::Polynome :: CONTAINER_ITERATOR i_spin  = i_begin;
    VA_CE::Polynome :: CONTAINER_ITERATOR i_last  = functional.end();
    int* i_var = pAtoms->AtomList;

    // copies atom type from IAGA to Lamarck
    if ( GA_style == MULTISTART )
      for( ; i_spin != i_last; ++i_spin, ++i_var ) 
        *i_spin = ( rand() -  ( (double) (RAND_MAX>>1) ) > 0.0 ) ? -1.0 : 1.0; 
    else
      for( ; i_spin != i_last; ++i_spin, ++i_var ) 
        *i_spin = ( *i_var == 31 ) ? -1.0 : 1.0;

    // minimizes
    minimizer.minimize();


    // sets  to physical and updates convex hull
    #ifndef WANG_CONSTRAINTS // strict S_i^2 - 1 = 0 constraints
      double result;
      functional.get_Obj1()->set_to_closest_constraints();

      result = functional.evaluate();

      std::vector<Ising_CE::Atom> :: iterator i_atom = structure.atoms.begin();
      i_spin = i_begin; i_var = pAtoms->AtomList;
      for( ; i_spin != i_last; ++i_spin, ++i_var, ++i_atom ) 
      {
        *i_var = ( abs(*i_spin -1.0 ) < ZERO_TOLERANCE ) ? 49 : 31;
        i_atom->type = *i_spin;
      }
      if ( convex_hull.add_structure(result, structure) )
        convex_hull_has_changed = true;
      return result;
    #else
      return find_best_vertex( pAtoms );
    #endif


  }

  double Lamarck :: no_optimize(AtomicConfig* pAtoms)
  {
    VA_CE::Polynome :: CONTAINER_ITERATOR i_spin  = functional.begin();
    VA_CE::Polynome :: CONTAINER_ITERATOR i_last  = functional.end();
    std::vector<Ising_CE::Atom> :: iterator i_atom = structure.atoms.begin();
    int* i_var = pAtoms->AtomList;
    double result;
    
    // copies atom type from IAGA to Lamarck
    for( ; i_spin != i_last; ++i_spin, ++i_var, ++i_atom ) 
    {
      *i_spin = ( *i_var == 31 ) ? -1.0 : 1.0;
      i_atom->type = *i_spin;
    }

    if ( CalculatingOldPop )
      return functional.evaluate();

    result = functional.evaluate();
    if ( convex_hull.add_structure(result, structure) )
      convex_hull_has_changed = true;

    return result;
  }


  inline double Lamarck :: evaluate(AtomicConfig* pAtoms)
  {
    if ( CalculatingOldPop )
      return no_optimize(pAtoms);

    switch ( GA_style )
    {
      case MULTISTART:
      case LAMARCKIAN: return optimize(pAtoms); 
      case LAMARCKIAN_EVOLUTION: 
        if ( nb_iterations == 0 )
          return no_optimize(pAtoms);
        else
          return optimize(pAtoms);
      default:             
      case DARWINISTIC: return no_optimize(pAtoms);
    }
  }

  #ifdef WANG_CONSTRAINTS
    // looks around current position, and finds the best vertex.
    // For each spin, if it is closer to 0 than +/-1, then both +/-1 vertices are
    // tested, otherwise, only the closest vertex is used
    // Note that it first checks all these structures for convex hull
    // breaking points, and then finds the best structure among them
    // according to their fitness (energy - convex_hull)
    // returns the energy of the best structure
    double Lamarck :: find_best_vertex(const AtomicConfig *pAtoms)
    {
      LOG_PROCEDURE("Lamarck :: find_best_vertex")
      opt :: Ndim_Iterator< int, std::less_equal<int> > global_iterator;
      FUNCTIONAL :: CONTAINER_ITERATOR i_begin = functional.begin();
      FUNCTIONAL :: CONTAINER_ITERATOR i_var = i_begin;
      FUNCTIONAL :: CONTAINER_ITERATOR i_end = functional.end();
      std::vector<Ising_CE::Atom> :: iterator i_begin_atom = structure.atoms.begin();
      std::vector<Ising_CE::Atom> :: iterator i_atom = i_begin_atom;
      double result, concentration;
      std::vector<double> energies, xcons;

      // constructs global iterator
      for ( ; i_var != i_end ; ++i_var )
      {
        double d_to_one = fabs( *i_var - 1.0 ) * weight;
        double d_to_mone = fabs( *i_var + 1.0 ) * weight;
        double d_to_zero = fabs( *i_var );
        if ( d_to_one < d_to_zero and d_to_one < d_to_mone )
           global_iterator.add( 1, 1);
        else if ( d_to_mone < d_to_zero and d_to_mone < d_to_one )
           global_iterator.add( 0, 0);
        else 
           global_iterator.add( 0, 1);
      }
  
      // first adds all structures to convex hull
      // and adds
      do 
      {
        LOG_PROCEDURE("Lamarck :: find_best_vertex -- evals")

        global_iterator.init_loop();
        concentration = 0.0;
        for( i_var = i_begin, i_atom = i_begin_atom; i_var != i_end; 
             ++i_var, ++i_atom, global_iterator.next_iterator() )
        {
          *i_var = ( global_iterator.get_current_iterator() == 1 ) ? 1.0 : -1.0;
          i_atom->type = *i_var;
          concentration += (*i_var);
        }

        result = functional.evaluate();
        energies.push_back(result); 
        xcons.push_back(concentration); 
        if ( convex_hull.add_structure(result, structure) )
          convex_hull_has_changed = true; 

      }
      while( ++global_iterator );


      // finds best vertex
      {
        // computes fitness energies
        std::vector<double> :: iterator i_e_begin = energies.begin();
        std::vector<double> :: iterator i_e = i_e_begin;
        std::vector<double> :: iterator i_e_end = energies.end();
        std::vector<double> :: iterator i_x = xcons.begin();
        double best_e = *i_e - convex_hull.evaluate(*i_x);
        result = *i_e;
        int best_vertex = 0, vertex = 0;
        for ( ++i_e, ++i_x; i_e != i_e_end; ++i_e, ++i_x, ++vertex)
          if ( *i_e - convex_hull.evaluate(*i_x) <  best_e )
          {
            result = *i_e;
            best_e = result - convex_hull.evaluate(*i_x);
            best_vertex = vertex;
          }
          
        // finally, goes through the global iterator loop again ... 
        global_iterator.reset();
        for( vertex = 0; vertex != best_vertex; ++vertex )
          ++global_iterator;

        // ... and copies best vertex to iaga
        int *iaga_atom = pAtoms->AtomList;
        global_iterator.init_loop();
        do 
        {
          *iaga_atom = ( global_iterator.get_current_iterator() == 1 ) ? 49 : 31;
          ++iaga_atom;
        }
        while ( global_iterator.next_iterator() );

      }
  
      return result;
    }
  #endif // WANG_CONSTRAINTS

  bool Lamarck :: Load( TiXmlHandle &handle )
  {
    TiXmlElement *child;
    rVector3d vec;

    // clusters, lattice, harmonics ....
    Functional_Builder :: Load (handle );
    
    // gets parameters
    random_range = 1;
    child = handle.FirstChild( "LaDa" ).FirstChild( "Random" ).Element();
    if ( child )
    {
      child->Attribute("range", &random_range);
      centered = child->Attribute("centered") ? true : false;
    }

    #ifdef WANG_CONSTRAINTS
      child = handle.FirstChild( "LaDa" ).FirstChild( "method" ).Element();
      if ( child )
        child->Attribute("value", &weight);
    #endif
  
    xml_filename = "convex_hull.xml";
    child = handle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
    if ( child and child->Attribute("xml") )
      xml_filename = child->Attribute("xml");
      
    xmgrace_filename = "convex_hull.agr";
    child = handle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
    if ( child and child->Attribute("xmgrace") )
      xmgrace_filename = child->Attribute("xmgrace");

    // finds the type of GA done here -- Darwinistic by default
    child = handle.FirstChild( "LaDa" ).FirstChild( "GA" ).Element();
    if ( child )
    {
      int d = 0;
      if ( child->Attribute("type", &d) )
        switch ( d )
        {
          case LAMARCKIAN : GA_style = LAMARCKIAN; break;
          case LAMARCKIAN_EVOLUTION : GA_style = LAMARCKIAN_EVOLUTION; break;
          case MULTISTART : GA_style = MULTISTART; break;
          default:
          case DARWINISTIC : GA_style = DARWINISTIC; break;
        }
    }
    // deletes content of xmgrace file
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
    switch( GA_style )
    {
      case LAMARCKIAN : xmgrace_file << "# Lamarckian GA" << std::endl; break;
      case LAMARCKIAN_EVOLUTION : xmgrace_file << "# Lamarckian Evolution GA" << std::endl; break;
      case MULTISTART : xmgrace_file << "# Multistart" << std::endl; break;
      default:
      case DARWINISTIC : xmgrace_file << "# Darwinistic GA" << std::endl; break;
    }
    xmgrace_file.flush();
    xmgrace_file.close();

    // reads structure from input
    child = handle.FirstChild( "LaDa" ).FirstChild( "Structure" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Structure in " << filename << std::endl;
      return false;
    }
    structure.Load(child, *axes);

    return true;
  }

  bool Lamarck :: read_CH()
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

  // sets up functional fitness and minimizer
  void Lamarck :: iaga_setup(SolverContext* solver)
  {
    // generate functional
    add_equivalent_clusters();
    generate_functional( structure, functional );

    // fitness pointers 
    fitness.set_baseline( &convex_hull );
    fitness.set_quantity( &functional ); 

    // minimizer
    minimizer.set_object( fitness );
  }

  void Lamarck :: print_xml()
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

  void Lamarck :: print_xmgrace()
  {
     std::ofstream xmgrace_file( xmgrace_filename.c_str(), std::ios_base::out|std::ios_base::app ); 
     xmgrace_file << " # iteration nb " << nb_iterations << std::endl;
     convex_hull.print_out(xmgrace_file, CONVEX_HULL::PRINT_XMGRACE);
     xmgrace_file.flush();
     xmgrace_file.close();
  }

  void Lamarck :: print_out( int format = XMGRACE_FORMAT )
  {
    switch (format) 
    {
      case XML_FORMAT: print_xml(); break;
      default :
      case XMGRACE_FORMAT: print_xmgrace(); break;
    }
  }




  // *****************************
  // lamarck_interface subroutines
  // *****************************

  int VASolverSetup(SolverContext *solver);
  void VASolverCleanup(SolverContext *solver);
  void FillSolverCtxVA(SolverContext *ctx);
  void va_functional ( SolverContext *solver, AtomicConfig *pAtoms, 
                       double *functional_p, MPI_Comm *comm_funct, 
                       int *error, char *tag);

  extern "C" void FillSolverCtxVA(SolverContext *ctx)	
  {
      ctx->fnSolverSetup = (FnSolverSetup)VASolverSetup;
      ctx->fnSolverCleanup = (FnSolverCleanup)VASolverCleanup;
      ctx->fnSolverEvaluateAtoms = (FnSolverEvaluateAtoms)va_functional;
      ctx->fnSolverEvaluate = NULL;
  }


  Lamarck *lamarck = NULL;
  int VASolverSetup(SolverContext *solver)
  {
    // destroys previous solver if any
    if ( lamarck )
      delete lamarck; 
    
    // creates new solver, and reads input
    lamarck = new Lamarck("input.xml");
    if ( not lamarck )
      return 1;
    
    lamarck->iaga_setup(solver);

    return 0;
  }

  void VASolverCleanup(SolverContext *solver)
  {
    if ( lamarck )
      delete lamarck;
    lamarck = NULL;
  }

  void va_functional ( SolverContext *solver, AtomicConfig *pAtoms, 
                       double *functional_p, MPI_Comm *comm_funct, 
                       int *error, char *tag)
  {
    int mynode_funct, mynode_world;
    MPI_Comm_rank(*comm_funct,&mynode_funct);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode_world);

    *functional_p = lamarck->evaluate(pAtoms);
  }

  extern "C" double get_config_concentration ( AtomicConfig *pAtoms )
  {
    int *i_atom = pAtoms->AtomList; 
    int *i_last = i_atom + (pAtoms->NumAtoms >> 1);
    int result=0;
    for ( ; i_atom != i_last; ++i_atom )
      ( *i_atom == 31 ) ? --result : ++result;
    return ( (double) result / ( (double) (pAtoms->NumAtoms >> 1) ) );
  }

  extern "C" double evaluate_convex_hull( double concentration )
  {
     return lamarck->evaluate_convex_hull( concentration ); 
  }

  extern "C" int convex_hull_has_changed()
  {
    if ( lamarck->convex_hull_has_changed )
    {
      lamarck->convex_hull_has_changed = false;
      return true;
    }
    return false;
  }

  extern "C" void print_convex_hull( int format = Lamarck :: XMGRACE_FORMAT )
  {
    lamarck->print_out(format);
  }

  extern "C" void increment_iterations()
  { ++lamarck->nb_iterations; }
     
  extern "C" void set_recalculating(bool _recalc = true)
  { lamarck->CalculatingOldPop = _recalc; }
#endif // _IAGA_LAMARCK_


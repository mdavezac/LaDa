//
//  Version: $Id: pescan.cc 278 2007-09-12 23:04:56Z davezac $
//
#include <functional>
#include <algorithm>
#include <ext/algorithm>
#include <fstream>
#ifndef __PGI
  #include<ext/functional>
  using __gnu_cxx::compose1;
#else
  #include<functional>
  using std::compose1;
#endif

#include <print/stdout.h>
#include <print/manip.h>
#include <lamarck/atom.h>
#include <opt/va_minimizer.h>

namespace BandGap
{
  template< class T_INDIVIDUAL, class T_INDIVTRAITS >
  bool Evaluator<T_INDIVIDUAL, T_INDIVTRAITS> :: 
           Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    t_Object &object = _indiv.Object();
    _node.Attribute("CBM", &object.cbm);
    _node.Attribute("VBM", &object.vbm);
    _indiv.quantities() = object.cbm - object.vbm; 

    return t_Base::Load( _indiv, _node, _type );
  }
  template< class T_INDIVIDUAL, class T_INDIVTRAITS >
  bool Evaluator<T_INDIVIDUAL, T_INDIVTRAITS > ::
    Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  {
    const t_Object &object = _indiv.Object();
    _node.SetDoubleAttribute("CBM", (double) object.cbm);
    _node.SetDoubleAttribute("VBM", (double) object.vbm);

    return t_Base::Save( _indiv, _node, _type );
  }
  template< class T_INDIVIDUAL, class T_INDIVTRAITS >
  bool Evaluator<T_INDIVIDUAL, T_INDIVTRAITS> :: Load( const TiXmlElement &_node )
  {
    if ( not t_Base::Load( _node ) )
    {
      std::cerr << " Could not load TwoSites::Evaluator<Object> input!! " << std::endl; 
      return false;
    }
    if ( not vff.Load( _node ) )
    {
      std::cerr << " Could not load vff input!! " << std::endl; 
      return false;
    }
    if ( not vff.initialize_centers() )
    {
      std::cerr << " Could not initialize Atomic_Center list in vff!! " << std::endl
                << " Are you sure the lattice and the structure correspond? " << std::endl; 
      return false;
    }
    if ( not vff_minimizer.Load( _node ) )
    {
      std::cerr << " Could not load vff minimizer from input!! " << std::endl;
      return false;
    }
    // wrong order! just to check whether Refs are read from input
    pescan.set_references( -666.666, 666.666 );
    if ( not pescan.Load( _node ) )
    {
      std::cerr << " Could not load pescan interface from input!! " << std::endl; 
      return false;
    }
    types::t_real x, y;
    pescan.get_references(x,y); // band edges have not been read if below is true
    if (     std::abs(x + 666.666 ) < types::tolerance 
         and std::abs(y - 666.666 ) < types::tolerance )
      set_all_electron();

    if (     _node.FirstChildElement("Filenames") 
         and _node.FirstChildElement("Filenames")->Attribute("BandEdge") )
    {
      references_filename = _node.FirstChildElement("Filenames")->Attribute("BandEdge");
      references_filename = Print::reformat_home( references_filename );
    }

    return true;
  }

  template< class T_INDIVIDUAL, class T_INDIVTRAITS >
  void Evaluator<T_INDIVIDUAL, T_INDIVTRAITS> :: init( t_Individual &_indiv )
  {
    t_Base :: init( _indiv );
    // sets structure to this object 
    structure << *current_object;
  }

  template< class T_INDIVIDUAL, class T_INDIVTRAITS >
  void Evaluator<T_INDIVIDUAL, T_INDIVTRAITS> :: evaluate()
  {
    // Creates an mpi aware directory: one per proc
    std::ostringstream sstr; sstr << "escan" << nbeval; 
    ++nbeval;
#ifdef _MPI
    sstr << mpi::main.rank();
#endif
    pescan.set_dirname(sstr.str());

    // minimizes vff energy
    vff_minimizer.minimize();
    structure.energy = vff.energy();


    // creates an mpi aware file name for atomic configurations
    sstr.str("");
    sstr << "atom_config";
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
    // prints atomic configurations
    vff.print_escan_input(sstr.str());
    // tells pescan where to find atomic configurations
    pescan.set_atom_input( sstr.str() );

    // then evaluates band gap
    types::t_real result = pescan(structure);

    // copies band edges into object
    get_bands( current_object->vbm, current_object->cbm );
    current_individual->quantities() = current_object->cbm - current_object->vbm;

    // checks that non-zero band gap has been found
    if ( result < types::tolerance )
    {
      Print::out << " Found metallic or negative band gap!! " << result << "\n" 
                 << current_object << "\n"
                 << " Will Try and Recompute Band Gap \n";
      set_all_electron();
      evaluate();
      return;
    }
    if ( pescan.get_method() == Pescan::Interface::Escan::FOLDED_SPECTRUM ) return;
    pescan.set_method(); // resets to folded spectrum if necessary


    // writes referecnce
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return; // not root no read write
#endif 
    Print::out << " Writing band edges to file " << references_filename << "\n" 
               << " using band gap of object " <<  current_object << "\n";
    write_references();
  }

  template <class T_INDIVIDUAL, class T_INDIVTRAITS> 
  bool Evaluator<T_INDIVIDUAL,T_INDIVTRAITS> :: Continue()
  {
    // on first iteration, writes references... then read them on following iterations
    ++age;
    read_references();

#ifdef _MPI
    if ( not mpi::main.is_root_node() )
      return true;
#endif

    if ( check_ref_every != -1 ) // recomputes all electron references every so often, if required
    {
      if ( not ( age % check_ref_every ) )
        set_all_electron();
    }
    return true;
  }

  template <class T_INDIVIDUAL, class T_INDIVTRAITS> 
  void Evaluator<T_INDIVIDUAL,T_INDIVTRAITS> :: write_references()
  {
#ifdef _MPI 
    if ( not mpi::main.is_root_node() )
      return;
#endif
    std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() )
      return;
    types :: t_real a, b;
    pescan.get_references( a, b );
    file << a << "   "; if ( file.fail() ) return;
    file << b << std::endl; if ( file.fail() ) return;
    file.close();
    return;
  }
  template <class T_INDIVIDUAL, class T_INDIVTRAITS> 
  void Evaluator<T_INDIVIDUAL, T_INDIVTRAITS> :: read_references()
  {
    types :: t_real a, b;
#ifdef _MPI 
    mpi::BroadCast bc(mpi::main);
    if ( mpi::main.is_root_node() )
    {
#endif
      std::ifstream file( references_filename.c_str(), std::ios_base::in ); 
      if ( not file.is_open() ) goto failure;
      file >> a; if ( file.fail() ) goto failure;
      file >> b; if ( file.fail() ) goto failure;
      file.close();
#ifdef _MPI
    }
    bc << a << b << mpi::BroadCast::allocate 
       << a << b << mpi::BroadCast::broadcast
       << a << b;
#endif 

    if ( a >= b )  return;
    pescan.set_references( a, b ); 
    return;

failure:
#ifdef _MPI
    bc << a << a << mpi::BroadCast::allocate 
       << a << a << mpi::BroadCast::broadcast;
#endif
  }

} // namespace pescan




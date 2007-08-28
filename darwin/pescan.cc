//
//  Version: $Id$
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

#include "pescan.h"
#include "print_xmgrace.h"
#include "lamarck/atom.h"
#include "opt/va_minimizer.h"

namespace BandGap
{
  bool Evaluator :: Load( t_Individual &_indiv, const TiXmlElement &_node, bool _type )
  {
    t_Object &object = _indiv.Object();
    _node.Attribute("CBM", &object.CBM);
    _node.Attribute("VBM", &object.VBM);
    _indiv.quantities() = object.VBM - object.CBM; 

    return t_Base::Load( _indiv, _node, _type );
  }
  bool Evaluator :: Save( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const
  {
    const t_Object &object = _indiv.Object();
    _node.SetDoubleAttribute("CBM", (double) object.CBM);
    _node.SetDoubleAttribute("VBM", (double) object.VBM);

    return t_Base::Save( _indiv, _node, _type );
  }
  bool Evaluator :: Load( const TiXmlElement &_node )
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
    pescan.set_references( -666.666, 666.666 ); // wrong order! just to check whether Refs are read from input
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
      references_filename = reformat_home(_node.FirstChildElement("Filenames")->Attribute("BandEdge"));

    return true;
  }

  void Evaluator::init( t_Individual &_indiv )
  {
    t_Base :: init( _indiv );
    vff.set_variables( &current_object->bitstring );
  }
  void Evaluator::evaluate()
  {
    // sets structure to this object 
    structure << *current_object;

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
    types::t_real vbm, cbm;
    get_bands( vbm, cbm );
    current_object->VBM = vbm;
    current_object->CBM = cbm;
    current_individual->quantities() = cbm - vbm;

    // checks that non-zero band gap has been found
    if ( result < types::tolerance )
    {
      set_all_electron();
      evaluate();
      return;
    }
    if ( pescan.get_method() == Pescan::Interface::Escan::FOLDED_SPECTRUM ) return;
    pescan.set_method(); // resets to folded spectrum if necessary


    // writes referecnce
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return // not root no read write
#endif 
    write_references();
  }

  bool Evaluator :: Continue()
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

  void Evaluator::write_references()
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
  void Evaluator::read_references()
  {
    types :: t_real a, b;
#ifdef _MPI 
    if ( mpi::main.is_root_node() )
    {
#endif
      std::ifstream file( references_filename.c_str(), std::ios_base::in ); 
      if ( not file.is_open() )
        return;
      file >> a; if ( file.fail() ) return;
      file >> b; if ( file.fail() ) return;
      if ( a >= b )
       return;
      file.close();
#ifdef _MPI
    }
    mpi::BroadCast broadcast(mpi::main);
    broadcast.serialize(a);
    broadcast.serialize(b);
    broadcast.allocate_buffers();
    broadcast.serialize(a);
    broadcast.serialize(b);
    broadcast();
    broadcast.serialize(a);
    broadcast.serialize(b);
#endif 

    pescan.set_references( a, b ); 
    return;
  }

} // namespace pescan




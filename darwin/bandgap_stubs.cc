//
//  Version: $Id$
//
#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include "bandgap_stubs.h"

namespace BandGap
{
  bool Keeper :: Load ( const TiXmlElement &_node )
  {
    double d;

    if ( not _node.Attribute("cbm", &d ) ) goto errorout;
    cbm = types::t_real(d);
    if ( not _node.Attribute("vbm", &d ) ) goto errorout;
    vbm = types::t_real(d);

    return true;
errorout:
    std::cerr << "Could not Load BandGap::Keeper" << std::endl;
    return false;
  }
  bool Keeper :: Save( TiXmlElement &_node ) const
  {
    _node.SetDoubleAttribute("vbm", vbm );
    _node.SetDoubleAttribute("cbm", cbm );

    return true;
  }
  bool Darwin :: Load( const TiXmlElement &_node )
  {
    if ( not bandgap.Load( _node ) ) return false;
    
    if ( not _node.FirstChildElement("GA") ) return true;
    const TiXmlElement *child = _node.FirstChildElement("GA");
    if ( not _node.FirstChildElement("GA") ) return true;
    child = child->FirstChildElement("Filenames");
    for(; child; child = child->NextSiblingElement("Filenames") )
    {
      if( not child->Attribute("BandEdge") ) continue;
      references_filename = child->Attribute("BandEdge");
      references_filename = Print::reformat_home( references_filename );
      Print::out << "Will store Reference energies at: " << references_filename << "\n";
      Print::xmg << Print::Xmg::comment 
                 << "Will store Reference energies at: " << references_filename << Print::endl;
    }

    return true;
  }

  void Darwin::operator()()
  {
    // Creates an mpi aware directory: one per proc
    std::ostringstream sstr; sstr << "escan" << nbeval; 
    ++nbeval;
#ifdef _MPI
    sstr << mpi::main.rank();
#endif
    dirname =  sstr.str();
    bandgap.set_dirname( dirname );
    bandgap.set_atom_input( atomicconfig );

    // then evaluates band gap
    types::t_real result = bandgap(structure);

    // checks that non-zero band gap has been found (and non-negative!)
    if ( result < types::tolerance )
    {
      Print::out << " Found metallic or negative band gap!! " << result << "\n" 
                 << " Will Try and Recompute Band Gap \n";
      set_all_electron();
      Darwin::operator()();
      return;
    }
    if ( bandgap.get_method() == Pescan::Interface::FOLDED_SPECTRUM ) return;
    bandgap.set_method(); // resets to folded spectrum if necessary


    // writes referecnce
#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return; // not root no read write
#endif 
    Print::out << " Writing band edges to file " << references_filename << "\n";
    write_references();
  }

  bool Darwin :: Continue()
  {
    // on first iteration, writes references... then read them on following iterations
    ++age;
    read_references();

#ifdef _MPI
    if ( not mpi::main.is_root_node() ) return true;
#endif
    // recomputes all electron references every so often, if required
    if (     check_ref_every != -1 
         and age % check_ref_every == 0 )  set_all_electron();

    return true;
  }

  void Darwin::write_references()
  {
    types :: t_real a, b;

    bandgap.Eref = bandgap.bands;
    Print::out << "Reference Energies are: CBM=" << bandgap.Eref.cbm
               << ", VBM=" << bandgap.Eref.vbm << "\n";
    Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM=" << bandgap.Eref.cbm
               << ", VBM=" << bandgap.Eref.vbm << Print::endl;
#ifdef _MPI 
    if ( not mpi::main.is_root_node() ) return;
#endif
    std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() ) return;
    file << bandgap.Eref.cbm << "   "; if ( file.fail() ) return;
    file << bandgap.Eref.vbm << std::endl; if ( file.fail() ) return;
    file.close();
    return;
  }
  void Darwin::read_references()
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
       << a << b << mpi::BroadCast::clear;
#endif 

    if ( a >= b )  return;
    bandgap.Eref.vbm = a; 
    bandgap.Eref.cbm = b; 
    return;

failure:
#ifdef _MPI
    bc << a << a << mpi::BroadCast::allocate 
       << a << a << mpi::BroadCast::broadcast
       << a << a << mpi::BroadCast::clear;
#endif
    return;
  }


} // namespace Pescan




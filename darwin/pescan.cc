//
//  Version: $Id$
//
#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include "pescan.h"

namespace Pescan
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
    std::cerr << "Could not Load Pescan::Keeper" << std::endl;
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
    {
      set_all_electron();
      Print::out << "No reference energy on input, "
                 << "will first perform all electron calculation\n"; 
      Print::xmg << Print::Xmg::comment << "No reference energy on input, "
                 << "will first perform all electron calculation" << Print::endl; 
    }
    else
    {
      Print::out << "Reference Energies are: CBM=" << y
                 << ", VBM=" << x << "\n";
      Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM=" << y
                 << ", VBM=" << x << Print::endl;
    }

    if (     _node.FirstChildElement("Filenames") 
         and _node.FirstChildElement("Filenames")->Attribute("BandEdge") )
    {
      references_filename = _node.FirstChildElement("Filenames")->Attribute("BandEdge");
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
    pescan.set_dirname( dirname );
    pescan.set_atom_input( atomicconfig );

    // then evaluates band gap
    types::t_real result = pescan(structure);

    // checks that non-zero band gap has been found (and non-negative!)
    if ( result < types::tolerance )
    {
      Print::out << " Found metallic or negative band gap!! " << result << "\n" 
                 << " Will Try and Recompute Band Gap \n";
      set_all_electron();
      Darwin::operator()();
      return;
    }
    if ( pescan.get_method() == Pescan::Interface::Escan::FOLDED_SPECTRUM ) return;
    pescan.set_method(); // resets to folded spectrum if necessary


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
    pescan.get_references( a, b );

    Print::out << "Reference Energies are: CBM=" << b
               << ", VBM=" << a << "\n";
    Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM=" << b
               << ", VBM=" << a << Print::endl;
#ifdef _MPI 
    if ( not mpi::main.is_root_node() ) return;
#endif
    std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() ) return;
    pescan.get_references( a, b );
    file << a << "   "; if ( file.fail() ) return;
    file << b << std::endl; if ( file.fail() ) return;
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
    pescan.set_references( a, b ); 
    return;

failure:
#ifdef _MPI
    bc << a << a << mpi::BroadCast::allocate 
       << a << a << mpi::BroadCast::broadcast
       << a << a << mpi::BroadCast::clear;
#endif
    return;
  }

  void Darwin::operator<<( const Vff::Darwin &_vff )
  {
    // creates an mpi aware file name for atomic configurations
    std::ostringstream  sstr;
    sstr << "atom_config";
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
    // prints atomic configurations
    _vff.print_escan_input(sstr.str());
    // tells pescan where to find atomic configurations
    atomicconfig = sstr.str();
  }


} // namespace pescan




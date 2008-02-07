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

    // Sets to folded spectrum if Pescan::BandGap::Eref have been specified
    // correctly
    if( Fuzzy::gt( bandgap.Eref.cbm, bandgap.Eref.vbm ) )
    {
       bandgap.set_method(); 
       Print::xmg << Print::Xmg::comment 
                  << "Pescan References (from input): "
                  << bandgap.Eref.vbm << " " << bandgap.Eref.cbm 
                  << Print::endl;
       Print::out << "Pescan References (from input): "
                  << bandgap.Eref.vbm << " " << bandgap.Eref.cbm 
                  << "\n";
    }
    else
    {
       Print::xmg << Print::Xmg::comment 
                  << "No references given on input. " << Print::endl;
       Print::out << "No references given on input.\n"
                  << "Will start with an all-electron calculation of the band-gap\n";
    }

    
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
 
    // Writes band edges to file
    // In case the band edges are unknown (all-electron calculation)
    // just creates file.
    if ( bandgap.get_method() != Pescan::Interface::FOLDED_SPECTRUM )
    {
      __NOTMPIROOT( return true; )
      std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
      file.close();
    }
    else write_references();

    return true;
  }

  void Darwin::operator()()
  {
    // Creates an mpi aware directory: one per proc
    std::ostringstream sstr;
    sstr << "ESCAN" << nbeval __DOMPICODE( << "." << mpi::main.rank() ); 
    ++nbeval;
    dirname =  sstr.str();
    bandgap.set_dirname( dirname );
    bandgap.set_atom_input( atomicconfig );

    // then evaluates band gap
    types::t_real result = bandgap(structure);

    if ( bandgap.get_method() == Pescan::Interface::FOLDED_SPECTRUM ) return;
    bandgap.set_method(); // resets to folded spectrum if necessary
    bandgap.set_nbstates(1);
    bandgap.Eref = bandgap.bands;
    Print::out << "Reference Energies are: CBM=" << bandgap.Eref.cbm
               << ", VBM=" << bandgap.Eref.vbm << "\n";
    Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM=" << bandgap.Eref.cbm
               << ", VBM=" << bandgap.Eref.vbm << Print::endl;


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

    __NOTMPIROOT( return true; )

    // recomputes all electron references every so often, if required
    if (     check_ref_every != -1 
         and age % check_ref_every == 0 )  set_all_electron();

    return true;
  }

  void Darwin::write_references()
  {
    types :: t_real a, b;

    __NOTMPIROOT( return; )

    std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() ) return;
    file << bandgap.Eref.vbm << "   "; if ( file.fail() ) return;
    file << bandgap.Eref.cbm << std::endl; if ( file.fail() ) return;
    file.close();
    return;
  }

  void Darwin::read_references()
  {
    types :: t_real a, b;

    __NOTMPIROOT( goto broadcast; )

    { // only serial and root node
      std::ifstream file( references_filename.c_str(), std::ios_base::in ); 
      if ( not file.is_open() )     { a = b = 0; goto broadcast; }
      file >> a; if ( file.fail() ) { a = b = 0; goto broadcast; }
      file >> b; if ( file.fail() ) { a = b = 0; goto broadcast; }
      file.close();
    }

broadcast:
    __DOMPICODE( 
      mpi::BroadCast bc(mpi::main); 
      bc << a << b << mpi::BroadCast::allocate 
         << a << b << mpi::BroadCast::broadcast
         << a << b << mpi::BroadCast::clear;
    )

    if ( Fuzzy::geq(a, b) )  return;
    bandgap.Eref.vbm = a; 
    bandgap.Eref.cbm = b; 
  }

} // namespace Pescan




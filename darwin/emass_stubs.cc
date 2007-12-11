//
//  Version: $Id$
//
#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include "emass_stubs.h"



namespace eMassSL
{
  bool Keeper :: Load ( const TiXmlElement &_node )
  {
    std::istringstream tensor_txt;
    const TiXmlElement *child;
    double d;

    if ( not _node.Attribute("cbm", &d ) ) goto errorout;
    cbm = types::t_real(d);
    child = _node.FirstChildElement( "emass" );
    if ( not child  ) goto errorout;
    if ( not child->GetText() ) goto errorout;
    tensor_txt.str(child->GetText());

    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(0,0); 
    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(0,1); emass(1,0) = emass(0,1);
    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(0,2); emass(2,0) = emass(0,2);
    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(1,1); 
    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(1,2); emass(2,1) = emass(1,2);
    if( not tensor_txt.good() ) goto errorout;
    tensor_txt >> emass(1,2); 

    return true;
errorout:
    std::cerr << "Could not Load eMassSL::Keeper" << std::endl;
    return false;
  }
  bool Keeper :: Save( TiXmlElement &_node ) const
  {
    _node.SetDoubleAttribute("cbm", cbm );
    TiXmlElement *tensor_xml = new TiXmlElement("emass");
    if( not tensor_xml ) return false;
    std::ostringstream sstr;
    sstr << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,0) << " " 
         << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,1) << " " 
         << std::fixed << std::setw(8) << std::setprecision(4) << emass(0,2) << " "
         << std::fixed << std::setw(8) << std::setprecision(4) << emass(1,1) << " " 
         << std::fixed << std::setw(8) << std::setprecision(4) << emass(1,2) << " "
         << std::fixed << std::setw(8) << std::setprecision(4) << emass(2,2);
    TiXmlText *tensor_txt = new TiXmlText( sstr.str().c_str() );
    if( not tensor_txt )
    {
      delete tensor_xml;
      return false;
    }
    tensor_xml->LinkEndChild( tensor_txt );
    _node.LinkEndChild( tensor_xml );

    return true;
  }



  bool Darwin :: Load( const TiXmlElement &_node )
  {
    if ( not emass.Load( _node ) ) return false;
    
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
    emass.set_dirname( dirname );
    emass.set_atom_input( atomicconfig );

    // then evaluates band gap
    if ( not emass(structure) ) 
    {
      if( emass.get_method() == Pescan::Interface::ALL_ELECTRON )  return;
      Print::out << " Error while computing effective mass\n"
                 << " Will try all electron computation before giving up " 
                 << Print::endl;
      set_all_electron();
      Darwin::operator()();
    }
    if ( emass.get_method() == Pescan::Interface::FOLDED_SPECTRUM ) return;
    emass.set_method(); // resets to folded spectrum if necessary


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
    types :: t_real a;

    emass.set_reference( emass.eigenvalues.back() );
    Print::out << "Reference Energies are: CBM=" << emass.get_reference() << "\n";
    Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM=" 
               << emass.get_reference() << Print::endl;
#ifdef _MPI 
    if ( not mpi::main.is_root_node() ) return;
#endif
    std::ofstream file( references_filename.c_str(), std::ios_base::out | std::ios_base::trunc ); 
    if ( not file.is_open() ) return;
    file << emass.get_reference() << "   "; if ( file.fail() ) return;
    file.close();
    return;
  }
  void Darwin::read_references()
  {
    types :: t_real a;
#ifdef _MPI 
    mpi::BroadCast bc(mpi::main);
    if ( mpi::main.is_root_node() )
    {
#endif
      std::ifstream file( references_filename.c_str(), std::ios_base::in ); 
      if ( not file.is_open() ) goto failure;
      file >> a; if ( file.fail() ) goto failure;
      file.close();
#ifdef _MPI
    }
    bc << a << mpi::BroadCast::allocate 
       << a << mpi::BroadCast::broadcast
       << a << mpi::BroadCast::clear;
#endif 

    emass.set_reference(a); 
    return;

failure:
#ifdef _MPI
    bc << a << mpi::BroadCast::allocate 
       << a << mpi::BroadCast::broadcast
       << a << mpi::BroadCast::clear;
#endif
    return;
  }



} 





//
//  Version: $Id$
//
#ifdef _EMASS
  #include <vff/layered.h>
  #include <vff/va.h>
  #include "emass.h"
  typedef Pescan::eMassSL t_Pescan;
  typedef Vff::VABase<Vff::Layered> t_Vff;
#else
  #include "va.h"
  typedef Pescan::VirtualAtom t_Pescan;
  typedef Vff::VABase<Vff::Functional> t_Vff;
#endif
#include <tinyxml/tinyxml.h>
#include <print/manip.h>


bool evaluate( const TiXmlElement &_node,
               Ising_CE::Structure &_structure,
               t_Pescan &_pescan, t_Vff &_vff )
{
  if( not _structure.Load(_node) ) return false;
  _structure.set_site_indices();

  _vff.init(true);

#ifndef _EMASS
   Pescan::BandGap& bandgap = (Pescan::BandGap&) _pescan;
   bandgap.set_method( Pescan::Interface::ALL_ELECTRON );
  _structure.energy = _pescan.evaluate();
#else
  _vff.evaluate();
  _vff.print_escan_input();
  _pescan.set_method( Pescan::Interface::ALL_ELECTRON );
  _pescan( _structure );
#endif

  Ising_CE::Fourier( _structure.atoms.begin(), _structure.atoms.end(),
                     _structure.k_vecs.begin(), _structure.k_vecs.end() );

  return true;
}


int main(int argc, char *argv[]) 
{
  std::string filename("input.xml");
#ifdef _CHECK_RESULTS
  std::string checkfilename("");
#endif
  if( argc > 1 )
  {
    std::ostringstream sstr;
    for( types::t_int i = 1; i < argc; ++i )
      sstr << argv[i] << " "; 
    std::istringstream istr( sstr.str() );
    while ( istr.good() )
    {
      std::string is_op;
      istr >> is_op; is_op = Print::StripEdges( is_op );
      if( is_op.empty() ) continue;
      else if(     istr.good()
               and (is_op == "-i" or is_op == "--input") ) istr >> filename;
#ifndef _CHECK_RESULTS
      else if( is_op == "-h" or is_op == "--help" )
        std::cout << "Command-line options:\n\t -h, --help this message"
                  << "\n\t -i, --input XML input file (default: input.xml)\n\n";
#else
      else if (     istr.good()
                and (is_op == "-s" or is_op == "--save") ) istr >> checkfilename;
      else if( is_op == "-h" or is_op == "--help" )
        std::cout << "Command-line options:\n\t -h, --help this message"
                  << "\n\t -i, --input XML input file\n\n"
                  << "\n\t -s, --save XML file where GA results where saved"
                  << " (default: read from GA output tags )\n\n";
#endif

      else std::cerr << "I don't understand command-line option " << is_op << std::endl;
    }
    filename = Print::reformat_home( filename );
    if( filename != "input.xml" )
      std::cout << "Reading from input file " << filename << std::endl;
#ifdef _CHECK_RESULTS
    checkfilename = Print::reformat_home( checkfilename );
    if( checkfilename != filename )
      std::cout << "Reading from GA output file " << checkfilename << std::endl;
#endif
  }
  TiXmlDocument doc( filename.c_str() );
  
  if  ( !doc.LoadFile() )
  {
    std::cerr << "error while opening input file " << filename << std::endl
              << doc.ErrorDesc() << std::endl; 
    return false;
  }

  TiXmlHandle handle( &doc );

  Ising_CE::Lattice lattice;
  Ising_CE::Structure structure;
  Ising_CE::Structure::lattice = &lattice;

  // loads lattice
  TiXmlElement *child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
  if ( not child )
  {
    std::cerr << "Could not find Lattice in input" << std::endl;
    return false;
  }
  if ( not lattice.Load(*child) )
  {
    std::cerr << "Error while reading Lattice from input" << std::endl;
    return false;
  }

#ifdef _EMASS
  t_Pescan pescan;
  t_Vff vff( structure );
  child = handle.FirstChild( "Job" ).Element();
  if ( not vff.Load(*child) )
  {
    std::cerr << "Error while reading vff from input" << std::endl;
    return false;
  }

#else
  t_Pescan pescan( structure );
  t_Vff &vff = pescan.Vff();
#endif
  child = handle.FirstChild( "Job" ).Element();
  if ( not pescan.Load(*child) )
  {
    std::cerr << "Error while reading pescan from input" << std::endl;
    return false;
  }

#ifdef _CHECK_RESULTS
  if ( checkfilename == "" )
  { 
    child = handle.FirstChild( "Job" ).FirstChild( "GA" ).FirstChild( "Filenames" ).Element();
    for( ; child; child = child->NextSiblingElement("Filenames") )
    {
      if ( not child->Attribute("save")  ) continue;
      checkfilename = Print::reformat_home(child->Attribute("save"));
    }
  }
  if ( checkfilename != "" and checkfilename != filename )
    if( not doc.LoadFile(checkfilename.c_str()) )
    {
      std::cerr << "Could not loadfile " << checkfilename << std::endl;
      return false;
    }
#endif


#ifdef _CHECK_RESULTS
  TiXmlHandle docHandle( &doc ); 
  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "optimum" ).
                    FirstChild( "Indivividual" ).
                    FirstChild( "Structure" ).Element();
  if( evaluate( *child, structure, pescan, vff ) )
  {
    std::cout << "\n\n\nChecking Optimum\n";
    structure.print_out( std::cout ); 
    std::cout << "\n\n";
    structure.print_xcrysden( std::cout );
#ifndef _EMASS
    const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
    std::cout << "\nVBM: " << bandgap.bands.vbm
              << " -- CBM:" << bandgap.bands.cbm
              << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
#else
    std::cout << "Emass tensor:\n" << pescan.bands();
#endif
  }
  else  std::cout << "\n\n\nChecking Optimum: error ... \n";

  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "Individual").Element();
  for(; child; child = child->NextSiblingElement("Individual" ) )
  {
    if ( not child->FirstChildElement("Structure") ) continue;
#else
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  for(; child; child = child->NextSiblingElement("Structure" ) )
  {
#endif
    if( evaluate( *child, structure, pescan, vff ) )
    {
      std::cout << "\n\n\nNew Structure\n";
      structure.print_out( std::cout ); 
      std::cout << "\n\n";
      structure.print_xcrysden( std::cout );
#ifndef _EMASS
      const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
      std::cout << "\nVBM: " << bandgap.bands.vbm
                << " -- CBM:" << bandgap.bands.cbm
                << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
#else
      std::cout << "Emass tensor:\n" << pescan.tensor;
#endif
    }
    else  std::cout << "\n\n\nNew Structure: error ... \n";
  }

  return 0;
}

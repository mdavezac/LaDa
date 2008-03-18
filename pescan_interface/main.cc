//
//  Version: $Id$
//
#ifdef _EMASS
  #include <vff/layered.h>
  #include <vff/va.h>
  #include "emass.h"
  typedef Pescan::eMassSL t_Pescan;
  typedef Vff::VABase<Vff::Functional> t_Vff;
#define __PROGNAME__ "Empirical Pseudo-Potential Functional for 100 Effective Masses in GaInSb"
#else
  #include "va.h"
  typedef Pescan::VirtualAtom t_Pescan;
  typedef Vff::VABase<Vff::Functional> t_Vff;
#define __PROGNAME__ "Empirical Pseudo-Potential Functional"
#endif
#include <tinyxml/tinyxml.h>
#include <print/manip.h>
#include <revision.h>


void parse_cli( int argc, char *argv[],
                std::string &_filename, std::string &_checkfilename,
                bool &_docheckresult, bool &_doevaluate )
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
             and (is_op == "-i" or is_op == "--input") ) istr >> _filename;
    else if (     istr.good()
              and (is_op == "-s" or is_op == "--save") ) 
      { _docheckresult = true; istr >> _checkfilename; }
    else if (is_op == "-d" or is_op == "--donteval") _doevaluate = false;
    else if( is_op == "-h" or is_op == "--help" )
    {
      std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package.\n"
                << "Command-line options:\n\t -h, --help this message"
                << "\n\t -v, --version Subversion Revision and Package version"
                << "\n\t -i, --input XML input file"
                << "\n\t -s, --save XML file where GA results where saved"
                << " (default: read from GA output tags or structures from input)"
                << "\n\t -d, --donteval Won't run evaluations, just printouts structures\n\n";
      exit(1);
    }
    else if( is_op == "-v" or is_op == "--version" )
    {
      std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package\n"
                << "Subversion Revision: " << SVN::Revision << "\n\n"; 
      exit(1);
    }
  }
  _filename = Print::reformat_home( _filename );
  if( _filename != "input.xml" )
    std::cout << "Reading from input file " << _filename << std::endl;
  if( _docheckresult )
  { 
    _checkfilename = Print::reformat_home( _checkfilename );
    if( _checkfilename != _filename )
      std::cout << "Reading from GA output file " << _checkfilename << std::endl;
  }
}

#ifdef _EMASS
inline void operator<<( t_Pescan &_pescan, const t_Vff &_vff )
{
  // creates an mpi aware file name for atomic configurations
  std::ostringstream  sstr;
  sstr << "atom_config";
#ifdef _MPI
    sstr << "." << mpi::main.rank();
#endif
  // prints atomic configurations
  _vff.print_escan_input(sstr.str());
  // tells bandgap where to find atomic configurations
  _pescan.set_atom_input(sstr.str());
}
#endif

bool evaluate( const TiXmlElement &_node,
               Ising_CE::Structure &_structure,
               t_Pescan &_pescan, t_Vff &_vff,
               bool _doeval )
{
  if( not _structure.Load(_node) ) return false;
  std::cout << "Successfuly read structure input from file " << std::endl;
  _structure.set_site_indices();

  if( not _vff.init(true) )
  {
    std::cerr << "Error while initializing vff\n" 
              << "Will skip current structure" << std::endl;
    return false;
  }

#ifndef _EMASS
   Pescan::BandGap& bandgap = (Pescan::BandGap&) _pescan;
   bandgap.set_method( Pescan::Interface::ALL_ELECTRON );
   if( _doeval ) _structure.energy = _pescan.evaluate();
#else
   _vff.evaluate();
   _pescan << _vff;
   if( _doeval ) _structure.energy = _pescan( _structure );
#endif

  Ising_CE::Fourier( _structure.atoms.begin(), _structure.atoms.end(),
                     _structure.k_vecs.begin(), _structure.k_vecs.end() );

  return true;
}


int main(int argc, char *argv[]) 
{
  std::string filename("input.xml");
  std::string checkfilename("");
  bool do_check_results = false;
  bool do_evaluate = true;

  parse_cli( argc, argv, filename, checkfilename, do_check_results, do_evaluate );

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
  std::cout << "Successfuly read Lattice input from file " << filename << std::endl;

#ifdef _EMASS
  t_Pescan pescan;
  t_Vff vff( structure );
  child = handle.FirstChild( "Job" ).Element();
  if ( not vff.Load(*child) )
  {
    std::cerr << "Error while reading vff from input" << std::endl;
    return false;
  }
  std::cout << "Successfuly read vff input from file " << filename << std::endl;

#else
  t_Pescan pescan( structure );
  t_Vff &vff = pescan.Vff();
  __DIAGA( pescan.set_mpi( mpi::main ); )
#endif
  child = handle.FirstChild( "Job" ).Element();
  if ( not pescan.Load(*child) )
  {
    std::cerr << "Error while reading pescan from input" << std::endl;
    return false;
  }
  std::cout << "Successfuly read pescan input from file " << filename << std::endl;

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


  TiXmlHandle docHandle( &doc ); 
  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "optimum" ).
                    FirstChild( "Indivividual" ).
                    FirstChild( "Structure" ).Element();
  if( child and evaluate( *child, structure, pescan, vff, do_evaluate ) )
  {
    do_check_results = true;
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
    std::cout << "Emass tensor:\n" << pescan.tensor;
#endif
  }

  child = docHandle.FirstChild( "Restart" ).
                    FirstChild( "Results" ).
                    FirstChild( "Individual").Element();
  for(; child; child = child->NextSiblingElement("Individual" ) )
  {
    if ( not child->FirstChildElement("Structure") ) continue;
    if( not evaluate( *child, structure, pescan, vff, do_evaluate ) ) continue;
    do_check_results = true;
    
    std::cout << "\n\n\nNew Structure\n";
    structure.print_out( std::cout ); 
    std::cout << "\n\n";
    structure.print_xcrysden( std::cout );
    if ( not do_evaluate ) continue;
#ifndef _EMASS
    const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
    std::cout << "\nVBM: " << bandgap.bands.vbm
              << " -- CBM:" << bandgap.bands.cbm
              << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
#else
    std::cout << "Emass tensor:\n" << pescan.tensor;
#endif
  }
  if ( do_check_results ) return 0;
  child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
  for(; child; child = child->NextSiblingElement("Structure" ) )
  {
    if( not evaluate( *child, structure, pescan, vff, do_evaluate ) ) continue;
    
    std::cout << "\n\n\nNew Structure\n";
    structure.print_out( std::cout ); 
    std::cout << "\n\n";
    structure.print_xcrysden( std::cout );
    if ( not do_evaluate ) continue;
#ifndef _EMASS
    const Pescan::BandGap& bandgap = (const Pescan::BandGap&) pescan;
    std::cout << "\nVBM: " << bandgap.bands.vbm
              << " -- CBM:" << bandgap.bands.cbm
              << "    ---    Band Gap: " << bandgap.bands.gap() << std::endl;
#else
    std::cout << "Emass tensor:\n" << pescan.tensor;
#endif
  }

  return 0;
}

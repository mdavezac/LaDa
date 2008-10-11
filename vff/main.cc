//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <stdexcept>

#include <tinyxml/tinyxml.h>

#include <crystal/lattice.h>
#include <crystal/structure.h>
#include <revision.h>
#include <opt/debug.h>

#include <vff/va.h>
#ifdef _LAYERED
#include <vff/layered.h>
typedef Vff::VABase<Vff::Layered> t_Vff;
#define __PROGNAME__ "Valence Force Field Functional for Epitaxial Structures"
#else
#include <vff/functional.h>
typedef Vff::VABase<Vff::Functional> t_Vff;
#define __PROGNAME__ "Valence Force Field Functional"
#endif

#include <print/manip.h>

void parse_cli( int argc, char *argv[], std::string &_filename )
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
    else if( is_op == "-h" or is_op == "--help" )
    {
      std::cout << "\n" << __PROGNAME__ << " from the " << PACKAGE_STRING << " package.\n"
                << "Command-line options:\n\t -h, --help this message"
                << "\n\t -v, --version Subversion Revision and Package version"
                << "\n\t -i, --input XML input file (default: input.xml)\n\n";
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
}  
    
int main(int argc, char *argv[]) 
{  
  try
  {
    std::string filename("input.xml");
    parse_cli( argc, argv, filename );
    
    TiXmlElement *child;
    atat::rVector3d vec;
    Crystal::Structure structure;
    Crystal::Lattice lattice;
    
    Crystal :: Structure :: lattice = &lattice;
    
    TiXmlDocument doc( filename.c_str() );
    
    __DOASSERT( not doc.LoadFile(),
                "Error while opening input file.\n";)
    TiXmlHandle handle( &doc );
    
    // loads lattice
    child = handle.FirstChild( "Job" ).FirstChild( "Lattice" ).Element();
    __DOASSERT( not child, "Could not find Lattice tag in input file.\n")
    __DOASSERT( not lattice.Load(*child), 
                "Error while reading lattice description from input file.\n")
    
    // loads structure
    child = handle.FirstChild( "Job" ).FirstChild( "Structure" ).Element();
    __DOASSERT( not child, "No <Structure> tags were found in input file.\n")
    for (; child ; child = child->NextSiblingElement("Structure") )
    {
      __DOASSERT( not structure.Load(*child), 
                  "Error while reading Structure description from input.\n")
    
      // loads vff functional
      TiXmlElement *vff_xml = handle.FirstChild( "Job" ).Element();
      __DOASSERT( not vff_xml, "Could not find <Job> tag in input.\n")
      t_Vff vff(structure);
      __DOASSERT( not vff.Load(*vff_xml),
                  "Could not Load Valence Force Field Functional from input.\n")
      if( not vff.init(true) ) continue;
      
      structure.energy = vff.evaluate() / 16.0217733;
      const atat::rMatrix3d stress = vff.Vff().get_stress();
      std::cout << std::fixed << std::setprecision(5) 
                << "Energy [eV]: " << std::setw(12) << structure.energy << std::endl
                << "Stress Tensor: " << std::endl 
                << std::setw(12) << stress(0,0) << " " << std::setw(12)
                                 << stress(1,0) << " " << std::setw(12)
                                 << stress(2,0) << std::endl
                << std::setw(12) << stress(0,1) << " " << std::setw(12)
                                 << stress(1,1) << " " << std::setw(12)
                                 << stress(2,1) << std::endl
                << std::setw(12) << stress(0,2) << " " << std::setw(12)
                                 << stress(1,2) << " " << std::setw(12)
                                 << stress(2,2) 
                << std::endl << std::endl << std::endl;
      vff.print_escan_input( "atomic.config" );
    
    
      Crystal::Fourier( structure.atoms.begin(),  structure.atoms.end(),
                         structure.k_vecs.begin(), structure.k_vecs.end() );
      structure.print_out(std::cout );
    
    }
  }
  catch ( std::exception &_e )
  {
    std::cerr << "Caught error while running " << __PROGNAME__ 
              << "\n" << _e.what() << "\n";
  }

  return 0;
}

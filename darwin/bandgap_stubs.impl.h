//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <print/stdout.h>
#include <print/manip.h>
#include <print/xmg.h>

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/string.hpp>
#include <mpi/mpi_object.h>
#include <opt/tinyxml.h>

namespace LaDa
{
  namespace BandGap
  {
#   if defined( DARWINHEAD ) || defined(INDARWIN)
#     error "Macros with same names."
#   endif
#   define DARWINHEAD \
      Darwin<T_VFF> 
#    define INDARWIN( var ) \
       template< class T_VFF >  var DARWINHEAD

    INDARWIN( bool ) :: Load( const TiXmlElement &_node )
    {
      if ( not bandgap.Load( _node ) ) return false;

      // Sets to folded spectrum if Pescan::BandGap::Eref have been specified
      // correctly
      if( Fuzzy::gt( bandgap.BandGap().Eref.cbm, bandgap.BandGap().Eref.vbm ) )
      {
         bandgap.BandGap().set_method(); 
         Print::xmg << Print::Xmg::comment 
                    << "Pescan References (from input): "
                    << bandgap.BandGap().Eref.vbm << " "
                    << bandgap.BandGap().Eref.cbm 
                    << Print::endl;
         Print::out << "Pescan References (from input): "
                    << bandgap.BandGap().Eref.vbm << " "
                    << bandgap.BandGap().Eref.cbm 
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
        Print::out << "Will store Reference energies at: "
                   << references_filename << "\n";
        Print::xmg << Print::Xmg::comment 
                   << "Will store Reference energies at: "
                   << references_filename << Print::endl;
      }
   
      // Writes band edges to file
      // In case the band edges are unknown (all-electron calculation)
      // just creates file.
      if ( bandgap.BandGap().get_method() != Pescan::Interface::FOLDED_SPECTRUM )
      {
        __MPICODE( if( not MPI_COMM.rank() ) return true; )
        std::ofstream file( references_filename.c_str(),
                            std::ios_base::out | std::ios_base::trunc ); 
        file.close();
      }
      else write_references();

      return true;
    }


    INDARWIN( void ) :: operator()( LaDa::GA::Keepers::BandGap &_keeper )
    {
      Darwin::operator()();
      // copies band edges into object
      _keeper.vbm = bandgap.BandGap().bands.vbm; 
      _keeper.cbm = bandgap.BandGap().bands.cbm;
      _keeper.energy = structure.energy;
      bandgap.get_stress(_keeper.stress);
    }

    INDARWIN( void ) :: operator()()
    {
      // Creates an mpi aware directory: one per proc
      std::ostringstream sstr;
      sstr << "ESCAN" << nbeval __MPICODE( << suffix ); 
      ++nbeval;
      dirname =  sstr.str();
      bandgap.BandGap().set_dirname( dirname );

      // then evaluates band gap
      types::t_real result = bandgap.evaluate();

      if ( bandgap.BandGap().get_method() == Pescan::Interface::FOLDED_SPECTRUM )
        return;
      bandgap.BandGap().set_method(); // resets to folded spectrum if necessary
      bandgap.BandGap().set_nbstates(1);
      bandgap.BandGap().Eref = bandgap.BandGap().bands;
      Print::out << "Reference Energies are: CBM=" << bandgap.BandGap().Eref.cbm
                 << ", VBM=" << bandgap.BandGap().Eref.vbm << "\n";
      Print::xmg << Print::Xmg::comment << "Reference Energies are: CBM="
                 << bandgap.BandGap().Eref.cbm
                 << ", VBM=" << bandgap.BandGap().Eref.vbm << Print::endl;


      // writes referecnce
      __MPICODE( if( not MPI_COMM.rank() ) return; )

      Print::out << " Writing band edges to file " << references_filename << "\n";
      write_references();
    }

    INDARWIN( bool ) :: Continue()
    {
      // on first iteration, writes references... then read them on following iterations
      ++age;
      read_references();

      __MPICODE( if( not MPI_COMM.rank() ) return true; )

      // recomputes all electron references every so often, if required
      if (     check_ref_every != -1 
           and age % check_ref_every == 0 )  set_all_electron();

      return true;
    }

    INDARWIN( void ) :: write_references()
    {
      types :: t_real a, b;

      __MPICODE( if( not MPI_COMM.rank() ) return; )

      std::ofstream file( references_filename.c_str(),
                          std::ios_base::out | std::ios_base::trunc ); 
      if ( not file.is_open() ) return;
      file << bandgap.BandGap().Eref.vbm << "   "; if ( file.fail() ) return;
      file << bandgap.BandGap().Eref.cbm << std::endl; if ( file.fail() ) return;
      file.close();
      return;
    }

    INDARWIN( void ) :: read_references()
    {
      types :: t_real a, b;

      __MPICODE( if( not MPI_COMM.rank() ) goto broadcast; )

      { // only serial and root node
        std::ifstream file( references_filename.c_str(), std::ios_base::in ); 
        if ( not file.is_open() )     { a = b = 0; goto broadcast; }
        file >> a; if ( file.fail() ) { a = b = 0; goto broadcast; }
        file >> b; if ( file.fail() ) { a = b = 0; goto broadcast; }
        file.close();
      }

      broadcast:
#     ifdef _MPI
        mpi::Pair< types::t_real, types::t_real > pair(a,b);
        boost::mpi::broadcast( MPI_COMM, pair, 0 );
#     endif

      if ( Fuzzy::geq(a, b) )  return;
      bandgap.BandGap().Eref.vbm = a; 
      bandgap.BandGap().Eref.cbm = b; 
    }

#   undef DARWINHEAD
#   undef INDARWIN
  } // namespace Pescan
} // namespace LaDa




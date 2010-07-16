#ifndef _LADA_CRYSTAL_READ_STRUCTURE_H_
#define _LADA_CRYSTAL_READ_STRUCTURE_H_

#include "LaDaConfig.h"

#include <vector>
#include <ostream>
#include <fstream>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <opt/types.h>

#include "structure.h"
#include "fill_structure.h"
#include "smith.h"

namespace LaDa 
{
  namespace Crystal
  {
    //! Reads structure in NREL format.
    void read_structure( Structure &_struct, const boost::filesystem::path &_path );

    //! \brief reads lda energies and structures from NREL input files.
    //! \param[in] _path is the full or relative path to the "LDAs.dat" file.
    //!                  Structure files are expected to be found in the same
    //!                  directory as the LDAs.dat and is deduced from \a _path.
    //! \param[inout] _structures are added to this container.
    void read_ce_structures( const boost::filesystem::path &_dir,
                             std::vector<Crystal::Structure> &_structures );
    //! \brief Reads one structure from input file.
    //! \details Expects NREL pifile format. Returns true if succesfull, false if
    //!          reached eof.
    bool read_pifile_structure( std::istream &_sstr,
                                Crystal::Structure &_structure );
    //! Computes predictions for all structures in PI file.
    template< class T_FUNCTIONAL >
      void enumerate_pifile( const std::string &_file, T_FUNCTIONAL &_op );

    template< class T_FUNCTIONAL >
      void enumerate_pifile( const std::string &_file, T_FUNCTIONAL &_op )
      {
        __DEBUGTRYBEGIN
        Crystal :: Structure structure;
        std::ifstream file( _file.c_str(), std::ifstream::in );
        size_t i(0);
        do
        {
          if( not Crystal :: read_pifile_structure( file, structure ) ) continue;
          std::cout << "    @" << structure.name << "  " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
          foreach( Crystal::Structure::t_Atom &atom, structure.atoms )
            atom.type = math::gt( atom.type, 0e0 ) ? -1e0: 1e0;
          std::cout << "    @-" << structure.name << " " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
        }
        while( not file.eof() );
        __DEBUGTRYEND(, "Error while enumerating pifile.\n" )
      }

    template< class T_FUNCTIONAL >
      void enumerate_gusfile( const std::string &_file, T_FUNCTIONAL &_op )
      {
        __DEBUGTRYBEGIN
        Crystal :: Structure structure;
        __DOASSERT( structure.lattice == NULL, "Lattice not set in structure.\n" );
        std::ifstream file( _file.c_str(), std::ifstream::in );
        std::string line;
        do
        {
          getline(file, line);
          if( line.find("#tot") != std::string::npos ) break;
        } while( (not file.eof()) and file.good() );

        if( file.bad() ) return;
        do
        {
          getline(file, line);
          boost::algorithm::trim(line); 
          if( line[0] == '#' ) continue;
          if( line.size() == 0 ) continue;
          std::istringstream input( line );
          types::t_int dummy;
          input >> structure.name;
          for( size_t i(0); i < 6; ++i ) input >> dummy;
          structure.cell = Eigen::Matrix3d::Zero();
          input >> structure.cell(0,0) 
                >> structure.cell(1,0) >> structure.cell(1,1)
                >> structure.cell(2,0) >> structure.cell(2,1) >> structure.cell(2,2);
          for( size_t i(0); i < 9; ++i ) input >> dummy;
          input >> line;
          structure.cell = structure.lattice->cell * structure.cell;
          fill_structure( structure );
          __ASSERT( line.size() != structure.atoms.size(),
                    "labels and structure have different sizes.\n" )
          t_SmithTransform const transform = get_smith_transform( structure );
          Structure::t_Atoms::iterator i_atom = structure.atoms.begin();
          Structure::t_Atoms::iterator i_atom_end = structure.atoms.end();
          math::iVector3d const &smith( boost::tuples::get<1>(transform) );
          for(; i_atom != i_atom_end; ++i_atom )
          {
            math::iVector3d indices( get_smith_index(transform, i_atom->pos) );
            const size_t index
            ( 
              indices[2] + smith(2) * ( indices[1] + smith(1) * indices[0] )
            );
            i_atom->type = line[index] == '0' ? -1.0: 1.0; 
          }


          std::cout << "    @" << structure.name << "  " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
          
          for(i_atom = structure.atoms.begin(); i_atom != i_atom_end; ++i_atom )
            i_atom->type = i_atom->type > 0e0 ? -1e0: 1e0;
          std::cout << "    @-" << structure.name << " " 
                    << structure.get_concentration()
                    << " " << _op( structure ) << "\n";
        }
        while( (not file.eof()) and file.good() );
        __DEBUGTRYEND(, "Error while enumerating pifile.\n" )
      }


  } // namespace Crystal

} // namespace LaDa

#endif

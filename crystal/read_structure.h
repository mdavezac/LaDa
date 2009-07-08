//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_READ_STRUCTURE_H_
#define _LADA_CRYSTAL_READ_STRUCTURE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/lexical_cast.hpp>

#include <vector>
#include <ostream>
#include <fstream>
#include <string>


#include <opt/types.h>

#include "structure.h"
#include "fill_structure.h"

namespace LaDa 
{
  namespace Crystal {

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

  } // namespace Crystal

} // namespace LaDa

#endif

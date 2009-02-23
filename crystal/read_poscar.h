//
//  Version: $Id$
//

#ifndef _LADA_CRYSTAL_READ_POSCAR_H_
#define _LADA_CRYSTAL_READ_POSCAR_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <ostream>
#include <fstream>
#include <string>


#include <opt/types.h>

#include "structure.h"

namespace LaDa 
{
  namespace Crystal {

    //! Reads structure in NREL format.
    template<class T_TYPE> 
      void read_poscar( TStructure<T_TYPE> &_struct, 
                        const boost::filesystem::path &_path,
                        const std::vector< T_TYPE >& _types );

    template<class T_TYPE> 
      void read_poscar( TStructure<T_TYPE> &_struct, 
                        const boost::filesystem::path &_path,
                        const std::vector< T_TYPE >& _types )
      {
        __TRYBEGIN
        
        namespace fs = boost::filesystem;  
        __DOASSERT( not fs::exists( _path ), "Path " << _path << " does not exits.\n" )
        __DOASSERT( not( fs::is_regular( _path ) or fs::is_symlink( _path ) ),
                    _path << " is neither a regulare file nor a system link.\n" )
        std::ifstream file( _path.string().c_str(), std::ifstream::in );
        std::string line;
        
        // name
        std::getline( file, line );
        _structure.name();
      }

  } // namespace Crystal

} // namespace LaDa

#endif

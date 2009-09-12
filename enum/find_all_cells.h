//
//  Version: $Id$
//
#ifndef LADA_ENUM_FIND_ALL_CELLS_H_
#define LADA_ENUM_FIND_ALL_CELLS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <opt/debug.h>
#include <opt/types.h>
#include <atat/vectmac.h>


namespace LaDa
{
  //! \cond
  namespace Crystal
  {
    class Lattice;
  }
  //! \endcond

  namespace enumerate
  {
    //! \brief Finds all inequivalent supercells of size \a _nmax.
    void find_all_cells( Crystal::Lattice const &_lattice, size_t _nmax, 
                         std::vector<atat::rMatrix3d> &_out );

  }
}

#endif

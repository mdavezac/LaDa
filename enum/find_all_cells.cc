//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "find_all_cells.h"


namespace LaDa
{

  namespace enumerate
  {

    template<class T_CONTAINER>
      void find_all_cells_impl_( Crystal::Lattice const &_lattice, size_t _nmax, 
                                 T_CONTAINER &_out );

    void find_all_cells( Crystal::Lattice const &_lattice, size_t _nmax, 
                         std::vector<atat::rMatrix3d> &_out )
      { find_all_cells_impl_(_lattice, _nmax, _cells); }

  } // namespace enumerate

} // namespace LaDa

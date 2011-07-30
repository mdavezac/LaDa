#ifndef LADA_CRYSTAL_TRAITS_H
#define LADA_CRYSTAL_TRAITS_H

#include "LaDaConfig.h"

#include <vector>

namespace LaDa 
{
  namespace crystal
  {
    //! \cond
    template<class TYPE> class Atom;
    //! \endcond
    namespace traits
    {
      template<class T_TYPE> struct StructureData
      {
        //! \typedef Type of the species
        typedef T_TYPE t_Type;
        //! \typedef The type of the collection of atoms. 
        typedef Atom<t_Type> t_Atom;
        //! \typedef The type of the collection of atoms. 
        typedef std::vector<t_Atom> t_Atoms;
      };
    };
  }
}
#endif

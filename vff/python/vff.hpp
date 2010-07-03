#ifndef _LADA_PYTHON_VFF_HPP_
#define _LADA_PYTHON_VFF_HPP_

#include "LaDaConfig.h"

#include <utility>

#include "../functional.h"
#include "../layered.h"
#include "../va.h"

namespace LaDa
{
  namespace python
  {
    void expose_vff();
    void expose_layeredvff();

    //! A tuple class to hold the vff functional and a structure.
    template< class T >
      struct Pair
      {
        //! Type of the functional
        typedef T first_type;
        //! Type of the structure.
        typedef Crystal::Structure second_type;
        //! Structure.
        second_type second;
        //! Functional.
        first_type first;
        //! Initialization.
        Pair() : second(), first(second) {}
        //! Initialization.
        Pair(first_type const &_f, second_type &_s) : second(_s), first(second)
          { first.copy_parameters(_f); }
        //! Copy.
        Pair(Pair const &_c) : second(_c.second), first(second) {}
      };

    //! Vff "functional" for python.
    typedef Pair< Vff::VABase<Vff::Functional> >  t_Vff;
    //! Vff "functional" for python.
    typedef Pair< Vff::VABase<Vff::Layered> >  t_LayeredVff;

    template<class T>
      void print_escan_input( Pair<T> &_self, const std::string &_path,
                              Crystal::TStructure<std::string> const &_str) 
      {
        Crystal::convert_string_to_real_structure(_str, _self.second);
        _self.first.init(true, false);
        _self.first.print_escan_input( _path ); 
      }
  }
} // namespace LaDa
#endif

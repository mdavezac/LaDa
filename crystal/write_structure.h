#ifndef LADA_CRYSTAL_WRITE_STRUCTURE_H
#define LADA_CRYSTAL_WRITE_STRUCTURE_H

#include "LaDaConfig.h"

#include <sstream>

#include <boost/type_traits/is_fundamental.hpp>

#include "exceptions.h"
#include "structure.h"

namespace LaDa 
{
  namespace crystal
  {
    namespace details
    {
      template<T_TYPE> size_t size(std::string const & _in, boost::false_type) { return 1; }
      template<T_TYPE> size_t size(std::vector<T_TYPE> const & _in, boost::false_type) { return _in.size(); }
      template<T_TYPE> size_t size(std::list<T_TYPE> const & _in, boost::false_type) { return _in.size(); }
      template<T_TYPE> size_t size(T_TYPE const & _in, boost::true_type) { return 1; }

      template<T_TYPE> T_TYPE const& first(std::string const & _in, boost::false_type) { return _in; }
      template<T_TYPE> typename T_TYPE::reference size(std::vector<T_TYPE> const & _in, boost::false_type)
        { return _in.first(); }
      template<T_TYPE> typename T_TYPE::reference size(std::list<T_TYPE> const & _in, boost::false_type)
        { return _in.first(); }
      template<T_TYPE> T_TYPE const& first(T_TYPE const & _in, boost::true_type) { return _in; }
    }
    //! \brief Returns structure in xcrysden format.
    //! This object throws if there are more than one specie per atom (e.~.g~ a
    //! lattice).
    template<class T_TYPE>
      std::string xcrysden(TemplateStructure<T_TYPE> const &_in)
      {
        typedef typename TemplateStructure<T_TYPE>::const_iterator iterator;
        iterator i_first = _in.begin();
        iterator const i_end = _in.end();
        std::ostringstream result;

        result << "CRYSTAL\nPRIMVEC\n"
               << (_in.scale*_in.cell.transpose())
               << "\nPRIMCOORD\n" 
               << _in.size() << " 1 \n";  
        for(; i_first != i_end; ++i_first)
        {
          if(not details::size(_in, boost::is_fundamental<T_TYPE>()) == 1)
            BOOST_THROW_EXCEPTION(error::too_many_species());
          _int << " " << Physics::Atomic::Z(first(i_first->type))
               << " " << ( (*i_first)[0] * i_first->scale ) << " "
               << " " << ( (*i_first)[1] * i_first->scale ) << " "
               << " " << ( (*i_first)[2] * i_first->scale ) << "\n";
        }
        return sstr.string();
      }
  }
}

#endif

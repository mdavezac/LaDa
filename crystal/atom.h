#ifndef LADA_CRYSTAL_ATOM_H
#define LADA_CRYSTAL_ATOM_H

#include "LaDaConfig.h"

#include <boost/serialization/base_object.hpp>

#include "atom_base.h"
#include "atom_freeze.h"

namespace LaDa
{
  namespace Crystal
  {
    //! \brief Describes an atom where the type is a vector.
    //! \details An atom consists of a position, a type, and frozen status
    //!          variable. The position should always be in cartesian units. The
    //!          type can be anything, from a string with the symbol of the atom,
    //!          to an double wich codes for the atomic type somehow, to a vector
    //!          of strings which describe the possible occupations of the atomic
    //!          position. To this end, the type is a template type \a T_TYPE.
    //!          The frozen status variable indicate which, if any, coordinate
    //!          should not be touched deuring optimization. There are a number
    //!          of possibilities:
    //!            - frozen::NONE indicate that no coordinate is
    //!                                     frozen.
    //!            - frozen::X indicate that the cartesian x coordinate
    //!                                  is frozen.
    //!            - frozen::Y indicate that the cartesian y coordinate
    //!                                  is frozen.
    //!            - frozen::Z indicate that the cartesian z coordinate
    //!                                  is frozen.
    //!            - frozen::T indicate that the occupation is frozen.
    //!            - Any combination of the above.
    //!            .
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation).
    class Atom : public AtomBase< std::vector<std::string> >, public AtomFreezeMixin
    {
      friend class boost::serialization::access;
#     ifdef LADA_WITH_LNS
        //! To load and save to xml-like input.
        friend class load_n_save::access; 
#     endif
      public:
        //! Type of the atomic type.
        typedef std::vector<std::string> t_Type;
        
        //! Constructor
        Atom() : AtomBase<t_Type>(), AtomFreezeMixin() {};
        //! Constructor and Initializer
        explicit  Atom_Type( const math::rVector3d &_pos, t_Type _type) 
                    : AtomBase<t_Type>(_pos, _type), AtomFreezeMixin(_type) {};
        //! Copy Constructor
        Atom(const Atom &_c ) : AtomBase<t_Type>(_c), AtomFreezeMixin(_c) {};
    
      private:
        //! Serializes an atom.
        template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version);
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version);
#       endif
    };

    template< class ARCHIVE >
      void Atom :: serialize( ARCHIVE & _ar, const unsigned int _version)
      {
        namespace bs = boost::serialization;
        _ar & bs::base_object< AtomBase<t_Type> >(*this);
        _ar & bs::base_object<AtomFreezeMixin>(*this); 
      }
#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_ARCHIVE>
        bool Atom :: lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          namespace lns = LaDa :: load_n_save;
          typedef AtomBase<t_Type> t1;
          typedef AtomFreezeMixin  t2;
          return _ar & (lns::base<t1>(*this) << lns::base<t2>(*this)); 
        }
#   endif

    inline std::ostream& operator<< (std::ostream &stream, Atom const &_atom)
    {
      return stream << static_cast< AtomBase<std::string> >(_atom) 
                    << "  freeze: " << freeze;
    }

  } // namespace Crystal
} // namespace LaDa
  
#endif

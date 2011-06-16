#ifndef _ATOM_H_
#define _ATOM_H_

#include "LaDaConfig.h"

#include <boost/serialization/base_object.hpp>

# ifdef LADA_WITH_LNS
#  include <load_n_save/xpr/base.h>
# endif

#include "atom_base.h"

namespace LaDa
{
  namespace Crystal
  {

    //! \brief Describes an atom.
    //! \details An atom consists of a position, a type, and frozen status
    //!          variable. The position should always be in cartesian units. The
    //!          type can be anything, from a string with the symbol of the atom,
    //!          to an double wich codes for the atomic type somehow, to a vector
    //!          of strings which describe the possible occupations of the atomic
    //!          position. To this end, the type is a template type \a T_TYPE.
    //!          The frozen status variable indicate which, if any, coordinate
    //!          should not be touched deuring optimization. There are a number
    //!          of possibilities:
    //!            - Atom_Type::FREEZE_NONE indicate that no coordinate is
    //!                                     frozen.
    //!            - Atom_Type::FREEZE_X indicate that the cartesian x coordinate
    //!                                  is frozen.
    //!            - Atom_Type::FREEZE_Y indicate that the cartesian y coordinate
    //!                                  is frozen.
    //!            - Atom_Type::FREEZE_Z indicate that the cartesian z coordinate
    //!                                  is frozen.
    //!            - Atom_Type::FREEZE_T indicate that the occupation is frozen.
    //!            - Any combination of the above.
    //!            .
    //! \warning The default equality comparison operator compares positions only (not
    //!          occupation).
    class Atom : public AtomBase<std::string> 
    {
      friend class boost::serialization::access;
      public:
        //! The type of the occupation
        typedef T_TYPE t_Type;

        //! Tags to freeze cell coordinates.
        enum t_FreezeAtom
        {
          FREEZE_NONE =  0, //!< Freeze no atomic coordinate.
          FREEZE_X   =  1,  //!< Freeze x coordinate.
          FREEZE_Y   =  2,  //!< Freeze y coordinate.
          FREEZE_Z   =  4,  //!< Freeze z coordinate.
          FREEZE_T   =  8,  //!< Freeze type.
          FREEZE_CARTESIANS  =  7,  //!< Freeze cartesians coordinates.
          FREEZE_ALL =  15,  //!< Freeze all.
        };

      public:
        //! The frozen status
        types::t_unsigned freeze;
        
        //! Constructor
        Atom() : AtomBase<std::string>(), freeze(FREEZE_NONE) {};
        //! Constructor and Initializer
        explicit  Atom_Type( const math::rVector3d &_pos, t_Type _type) 
                    : AtomBase<std::string>(_pos, _type), freeze(FREEZE_NONE) {};
        //! Copy Constructor
        Atom(const Atom &_c ) : AtomBase<std::string>(_c), freeze(_c.freeze) {};

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
        _ar & bs::base_object< AtomBase<std::string> >(*this);
        _ar & freeze;
      }
#   ifdef LADA_WITH_LNS
      //! To load and save to xml-like input.
      template<class T_ARCHIVE>
        bool Atom :: lns_access(T_ARCHIVE &_ar, const unsigned int _version) 
        {
          namespace lns = LaDa :: load_n_save;
          std::map<std::string, LaDa::types::t_unsigned> freeze_map;
          freeze_map["none"] = FREEZE_NONE;
          freeze_map["x"] = FREEZE_X;
          freeze_map["y"] = FREEZE_Y;
          freeze_map["z"] = FREEZE_Z;
          freeze_map["t"] = FREEZE_T;
          freeze_map["cartesian"] = FREEZE_CARTESIANS;
          freeze_map["all"] = FREEZE_ALL;
          return _ar & 
                 ( 
                   lns::base< AtomBase<std::string> >(*_this) 
                    << lns::option("freeze", lns::action=lns::enum_(freeze, freeze_map),
                                        lns::default_=FREEZE_NONE)
                 );
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

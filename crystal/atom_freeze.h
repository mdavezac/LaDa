#ifndef LADA_CRYSTAL_ATOM_FREEZE_H
#define LADA_CRYSTAL_ATOM_FREEZE_H

#include "LaDaConfig.h"

# ifdef LADA_WITH_LNS
#  include <load_n_save/xpr/utilities.h>
#  include <load_n_save/action/enum.h>
# endif

namespace LaDa
{
  namespace crystal
  {
    //! \brief Mixin providing atoms with freezing variable.
    //! \details The frozen status variable indicate which, if any, coordinate
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
    struct AtomFreezeMixin
    {
      //! Namespace containing enums for freezing coordinates.
      struct frozen
      {
        //! Tags to freeze cell coordinates.
        enum type
        {
          NONE       =  0,  //!< Freeze no atomic coordinate.
          X          =  1,  //!< Freeze x coordinate.
          Y          =  2,  //!< Freeze y coordinate.
          Z          =  4,  //!< Freeze z coordinate.
          T          =  8,  //!< Freeze type.
          CARTESIANS =  7,  //!< Freeze cartesians coordinates.
          ALL        =  15  //!< Freeze all.
        };
      };

      //! The frozen status
      types::t_unsigned freeze;

      //! Constructor.
      AtomFreezeMixin(types::t_unsigned _in = frozen::NONE) : freeze(_in) {}
      //! Copy constructor.
      AtomFreezeMixin(AtomFreezeMixin const &_c) : freeze(_c.freeze) {}

      //! Serializes the status.
      template<class ARCHIVE> void serialize(ARCHIVE & _ar, const unsigned int _version)
        { _ar & freeze; }
#       ifdef LADA_WITH_LNS
          //! To load and save to xml-like input.
          template<class T_ARCHIVE> bool lns_access(T_ARCHIVE &_ar, const unsigned int _version)
          {
            namespace lns = LaDa :: load_n_save;
            std::map<std::string, LaDa::types::t_unsigned> freeze_map;
            freeze_map["none"]      = frozen::NONE;
            freeze_map["x"]         = frozen::X;
            freeze_map["y"]         = frozen::Y;
            freeze_map["z"]         = frozen::Z;
            freeze_map["t"]         = frozen::T;
            freeze_map["cartesian"] = frozen::CARTESIANS;
            freeze_map["all"]       = frozen::ALL;
            return _ar & ( lns::section("AtomFreeze") 
                << lns::option( "freeze",
                                lns::action=lns::enum_(freeze, freeze_map),
                                lns::default_=frozen::NONE ) );
          }
#       endif
      
    };
  }
}

#endif

#ifndef LADA_MATH_SERIALIZE_H
#define LADA_MATH_SERIALIZE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "eigen.h"

namespace boost {
  namespace serialization {

    //! Serializes eigen real vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::rVector3d & g, const unsigned int version)
     { ar & g.x(); ar & g.y(); ar & g.z(); }
     //! Serializes eigen integer vectors.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::iVector3d & g, const unsigned int version)
     { ar & g.x(); ar & g.y(); ar & g.z(); }
    //! Serializes eigen real matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::rMatrix3d & g, const unsigned int version)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          ar & g(i,j);
    }
    //! Serializes eigen integer matrices.
    template<class Archive>
    void serialize(Archive & ar, LaDa::math::iMatrix3d & g, const unsigned int version)
    {
      for(size_t i(0); i < 3; ++i)
        for(size_t j(0); j < 3; ++j)
          ar & g(i,j);
    }
  }
}

#endif

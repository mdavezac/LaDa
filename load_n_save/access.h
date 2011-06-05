#ifndef LADA_LNS_ACCESS_H
#define LADA_LNS_ACCESS_H

#include "LaDaConfig.h"

#include <boost/serialization/strong_typedef.hpp>
#include <boost/serialization/pfto.hpp>

namespace LaDa 
{
  namespace load_n_save
  {
    // Defines version type for overloading.
    typedef unsigned int version_type;
    typedef BOOST_PFTO version_type long_version_type;
    BOOST_STRONG_TYPEDEF(version_type, pfto_version_type);
    

    //! Intrusive external access. 
    struct access
    {
      template<class T_ARCHIVE, class T> 
        static inline bool lns_access(T_ARCHIVE &_ar, T &_var, version_type const _version)
          { return _var.lns_access(_ar, _version); }
    };
    
    // default implementation - call the member function "serialize"
    template<class Archive, class T>
      inline bool lns_access( Archive & ar, T & t, long_version_type const file_version )
        { access::lns_access(ar, t, static_cast<unsigned int>(file_version)); }

    // Taken from boost/serialization/serialization.hpp
    // layer 3 - move call into serialization namespace so that ADL will function
    // in the manner we desire.
    //
    // on compilers which don't implement ADL. only the current namespace
    // i.e. boost::serialization will be searched.
    // 
    // on compilers which DO implement ADL
    // serialize overrides can be in any of the following
    // 
    // 1) same namepace as Archive
    // 2) same namespace as T
    // 3) boost::serialization
    //
    // Due to Martin Ecker
    template<class Archive, class T>
      inline bool lns_access_adl( Archive & ar, T & t, version_type const file_version )
      { 
        // note usage of function overloading to delay final resolution
        // until the point of instantiation.  This works around the two-phase
        // lookup "feature" which inhibits redefintion of a default function
        // template implementation. Due to Robert Ramey
        //
        // Note that this trick generates problems for compiles which don't support
        // PFTO, suppress it here.  As far as we know, there are no compilers
        // which fail to support PFTO while supporting two-phase lookup.
        #if ! defined(BOOST_NO_FUNCTION_TEMPLATE_ORDERING)
            pfto_version_type const v(file_version);
            lns_access(ar, t, v);
        #else
            lns_access(ar, t, file_version);
        #endif
      }

  } // namespace load_n_save
} // namespace LaDa

#endif

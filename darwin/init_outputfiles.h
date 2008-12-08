//
//  Version: $Id: print_callbacks.h 871 2008-12-02 05:24:15Z davezac $
//
#ifndef _LADA_GA_START_XMG_OUT_H_
#define _LADA_GA_START_XMG_OUT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <boost/filesystem/path.hpp>

namespace LaDa
{

  namespace GA
  {
    //! Reads filenames, synchronizes between procs.
    void init_outputfiles( const boost::filesystem::path &_input );
  } // namespace GA
} // namespace LaDa

#endif

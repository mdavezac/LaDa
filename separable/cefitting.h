//
//  Version: $Id$
//
#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <istringstream>

#include <boost/filesystem.hpp>

#include <opt/types.h>
#include <boost/random/uniform_real.hpp>

#include "ce.h"

namespace Fitting
{

  class SepCeInterface
  {
    public:
      //! Constructor
      SepCeInterface() : generator( 42u ), uni_dist(0,1), rng(generator, uni_dist)
        { random_engine.seed( static_cast<unsigned int>(std::time(0)) ); }
      
      //! Reads the structures from input.
      void read( CE::SymSeparables &_symseps,
                 std::string _dir,
                 std::string _ldasdat = "LDAs.dat" );

      //! Fits a separable function to the complete training set.
      template< class T_ALLSQ >
      void fittoall(  T_ALLSQ &_llsq, CE::Separables &_sep ) const;

      //! Fits a separable function to the complete training set.
      template< class T_ALLSQ >
      void fitexcept( T_ALLSQ &_llsq, CE::Separables &_sep,
                      std::vector< types::t_unsigned > &_excpt ) const;

    protected:
      //! Reads structure names and energies.
      void read_ldasdat( boost::filesystem::path &_path,
                         const std::string _ldasdat );
      //! Reads structure names and energies.
      void read_structure( CE::SymSeparables &_symseps,
                           boost::filesystem::path &_path, 
                           const std::string _filename );


      //! The type of vector holding all equvivalent configurations.
      typedef CE::SymSeparables::t_Configurations t_Configurations;

      //! An optional weight vector.
      std::vector< types::t_real > weight;
      //! A vector of target values.
      std::vector< types::t_real > targets;
      //! A vector of structures.
      std::vector< Crystal::Structure > structures;
      //! The complete training set.
      std::vector< t_Configurations > training;
      //! The names of the structures.
      std::vectors< std::string > names;
      //! A random number generator.
      boost::random::mt11213b generator;
      //! A uniform distribution.
      boost::uniform_real<> uni_dist;
      //! A random number generator
      boost::variate_generator< booos::random::mt11213b&, 
                                boost::uniform_real<> > uni;
  };

}

#include "cefitting.impl.h"


#endif

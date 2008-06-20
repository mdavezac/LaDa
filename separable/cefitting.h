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
      //! Structures excluded from fit.
      std::vector< types::t_unsigned > exclude;
      //! A random number generator
      boost::variate_generator< booos::random::mt11213b&, 
                                boost::uniform_real<> > uni;

      //! Constructor
      SepCeInterface() : generator( 42u ), uni_dist(0,1), rng(generator, uni_dist)
        { random_engine.seed( static_cast<unsigned int>(std::time(0)) ); }
      
      //! Reads the structures from input.
      void read( CE::SymSeparables &_symseps,
                 std::string _dir,
                 std::string _ldasdat = "LDAs.dat" );

      //! Fits a separable function to the complete training set.
      template< class T_ALLSQ, class T_COLLAPSE>
      void fit( T_ALLSQ &_allsq, T_COLLAPSE &_collapse ) const;


      //! Check training convergence.
      std::pair< types::t_real, types::t_real>
        check_training( CE::Separables &_sep, bool _verbose = false ) const
        { return check( false, _verbose ); } 
      //! Check predictions.
      std::pair< types::t_real, types::t_real>
        check_predictions( CE::Separables &_sep, bool _verbose = false ) const
        { return check( true, _verbose ); } 

      //! Number of structures in training set.
      types::t_unsigned training_set_size() const
        { return training.size(); }


    protected:
      //! Reads structure names and energies.
      void read_ldasdat( boost::filesystem::path &_path,
                         const std::string _ldasdat );
      //! Reads structure names and energies.
      void read_structure( CE::SymSeparables &_symseps,
                           boost::filesystem::path &_path, 
                           const std::string _filename );
      //! \brief Check convergence.
      //! \param _which = true checks training convergence.
      //!               = false checks prediction.
      //! \param _verbose whether to print error for each structure.
      //! \return an std::pair where the first value is the average error,
      //!         and the second value the maximum error.
      std::pair< types::t_real, types::t_real>
        check( CE::Separables &_sep, bool _which , bool _verbose = false ) const;


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
  };

  template< class T_ALLSQ, class T_COLLAPSE>
    std::pair<types::t_real, types::t_real> 
      leave_one_out( SepCeInterface &_interface,
                    T_ALLSQ &_allsq, 
                    T_COLLAPSE &_collapse, 
                    bool _verbose = false );

}

#include "cefitting.impl.h"


#endif

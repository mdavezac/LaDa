//
//  Version: $Id$
//
#ifndef _SEPARABLE_CEFITTING_H_
#define _SEPARABLE_CEFITTING_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

#include <opt/types.h>

#include "ce.h"

namespace Fitting
{

  class SepCeInterface
  {
    public:
      typedef std::pair< types::t_real, types::t_real > t_PairErrors;
      //! Structures excluded from fit.
      std::vector< types::t_unsigned > exclude;
      //! Constructor
      SepCeInterface() : offset( 0e0 ) {}
      
      //! Reads the structures from input.
      void read( CE::SymSeparables &_symseps,
                 const std::string &_dir,
                 const std::string &_ldasdat = "LDAs.dat",
                 bool _verbose = false );

      //! Fits a separable function to the complete training set.
      template< class T_ALLSQ, class T_COLLAPSE>
      void fit( T_ALLSQ &_allsq, T_COLLAPSE &_collapse ) const;


      //! Check training convergence.
      std::pair< types::t_real, types::t_real>
        check_training( const CE::Separables &_sep, bool _verbose = false ) const
        { return check( _sep, false, _verbose ); } 
      //! Check predictions.
      std::pair< types::t_real, types::t_real>
        check_predictions( const CE::Separables &_sep, bool _verbose = false ) const
        { return check( _sep, true, _verbose ); } 

      //! Number of structures in training set.
      types::t_unsigned training_set_size() const
        { return training.size(); }

      //! Sets the offset. 
      void set_offset( types::t_real _offset ) { offset = _offset; }

    protected:
      //! Reads structure names and energies.
      void read_ldasdat( const boost::filesystem::path &_path,
                         const std::string& _ldasdat );
      //! Reads structure names and energies.
      void read_structure( const CE::SymSeparables &_symseps,
                           const boost::filesystem::path &_path, 
                           const std::string& _filename );
      //! \brief Check convergence.
      //! \param _which = true checks training convergence.
      //!               = false checks prediction.
      //! \param _verbose whether to print error for each structure.
      //! \return an std::pair where the first value is the average error,
      //!         and the second value the maximum error.
      std::pair< types::t_real, types::t_real>
        check( const CE::Separables &_sep,
               bool _which , bool _verbose = false ) const;


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
      std::vector< std::string > names;
      //! An energy offset to add to the structures.
      types::t_real offset;
  };

  template< class T_ALLSQ, class T_COLLAPSE>
    std::pair< SepCeInterface::t_PairErrors, SepCeInterface::t_PairErrors> 
      leave_one_out( SepCeInterface &_interface,
                    T_ALLSQ &_allsq, 
                    T_COLLAPSE &_collapse, 
                    bool _verbose = false );

}

#include "cefitting.impl.h"


#endif

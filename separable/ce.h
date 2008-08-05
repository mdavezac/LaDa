//
//  Version: $Id: separable.h 626 2008-06-18 22:05:05Z davezac $
//
#ifndef _SEPARABLE_CE_H_
#define _SEPARABLE_CE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <vector>

#include <opt/types.h>
#include <atat/vectmac.h>
#include <crystal/lattice.h>
#include <crystal/structure.h>

#include "separable.h"
#include "boolean.h"

namespace CE
{
  //! A separable function for a fixed-lattice incorporating all symmetry
  //! operations.
  class SymSeparables;
  //! A separable function for fixed-lattice, fixed orientation, fixed-origin.
  class Separables : public ::Separable::Function< 
                               std::vector< 
                                ::Separable::Summand<
                                   ::Separable::BooleanBasis > > >
  {
      friend class SymSeparables;
      //! Type of the base.
      typedef ::Separable::Function< std::vector< ::Separable::Summand<
                                     ::Separable::BooleanBasis > > > t_Base;
    public:

      //! Constructor.
      Separables( types::t_unsigned _rank = 2,
                  types::t_unsigned _size = 3, 
                  const std::string &basis_type = "cube");
      //! Constructor.
      Separables( types::t_unsigned _rank,
                  const atat::rMatrix3d &_cell, 
                  const std::string &basis_type = "parallelepiped");

      //! Sets the rank of the separable function.
      void set_rank( types::t_unsigned _rank );
      //! Sets the basis of the separable function.
      void set_basis( types::t_unsigned _size, const std::string &_type = "cube" );
      //! Sets the basis of the separable function.
      void set_basis( const atat::rMatrix3d& _cell, const std::string &_type = "cell" );
      //! Returns the rank of the separable function.
      types::t_unsigned rank() const  { return basis.size(); }
      //! Returns the type of the basis.
      std::string type() const  { return basis_type; }
      //! Returns the size of the basis.
      types::t_unsigned size() const  { return basis_size; }

    protected:
      //! Renames functions in sep. basis.
      void rename();
      //! Size of the basis.
      types::t_unsigned basis_size;
      //! Type of the basis.
      std::string basis_type;
      //! A vector of positions used in creating the basis.
      std::vector< atat::rVector3d > positions;
      //! Name of this function.
      std::string name;
  };

  // A separable function for a fixed-lattice incorporating all symmetry
  // operations.
  class SymSeparables
  {
    public:
      //! Type of the container of positions.
      typedef std::vector< atat::rVector3d > t_Basis;
      //! Type of the container of symmetry operations.
      typedef std::vector< atat::rMatrix3d > t_SymOps;
      //! Type of the pure bitset representing a configuration for a set symmetry.
      typedef std::vector< bool > t_Bitset;
      //! Type  containing a pure bitset and an attached coefficient.
      typedef std::pair< t_Bitset, types::t_real > t_CoefBitset;
      //! Type representing a set of configurations with all possible symmetries.
      typedef std::vector< t_CoefBitset > t_Configurations;

      //! Constructor.
      SymSeparables ( Separables &_sep );
      //! Constructor.
      SymSeparables   ( t_Basis &_poss, Crystal::Lattice &_lat )
                    : basis( _poss ) { init_syms( _lat ); }

      //! Creates all necessary configurations for a given structure.
      void configurations( const Crystal :: Structure &_structure,
                           SymSeparables :: t_Configurations& _confs ) const;
      //! Evaluates a set of configurations. 
      types::t_real operator()( t_Configurations &_conf,
                                const Separables &_func ) const;

    protected:
      //! Initializes list of symmetry operations.
      void init_syms ( Crystal::Lattice &_lat );
      
      //! A reference to the positions.
      t_Basis &basis;
      //! The syemmetry operations.
      t_SymOps syms;
  };

} // end of CE namespace

#endif //  _SEPARABLE_CE_H_

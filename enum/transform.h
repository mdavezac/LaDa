//
//  Version: $Id$
//
#ifndef LADA_ENUM_TRANSFORM_H_
#define LADA_ENUM_TRANSFORM_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <vector>
#include <utility>

#include <crystal/smith.h>
#include <crystal/symmetry_operator.h>

#include "numeric_type.h"

namespace LaDa
{
  namespace enumeration
  {
    //! \brief Symmetry operation of the lattice operating on an integer structure.
    //! \see Transform an integer structure to another. Appendix of PRB 80, 014120.
    class Transform : public Crystal::SymmetryOperator
    {
      //! Type of the container of supercell-independent elements.
      typedef std::pair<types::t_int, atat::rVector3d> t_Independent;
      //! Type of the container of supercell-independent elements.
      typedef std::vector<t_Independent> t_Independents;
      public:
        //! Copy Constructor
        Transform   (Crystal::SymmetryOperator const &_c,
                     Crystal::Lattice const &_lat);
        //! Copy Constructor
        Transform    (Transform const &_c)
                   : Crystal::SymmetryOperator(_c),
                     permutations_(_c.permutations_), 
                     independents_(_c.independents_), 
                     nsites_(_c.nsites_),
                     card_(_c.card_), is_trivial_(_c.is_trivial_) {}
        //! Initializes transform for specific supercell.
        void init(atat::rMatrix3d const &_left, atat::iVector3d const &_smith);
        //! Initializes transform for specific supercell.
        void init( Crystal::t_SmithTransform const &_t )
          { return init(boost::tuples::get<0>(_t), boost::tuples::get<1>(_t)); }
        //! Performs transformation.
        t_uint operator()(t_uint _x, FlavorBase const &_flavorbase) const;
        //! Is trivial.
        bool is_trivial() const { return is_trivial_; };

      private:
        //! Permutation.
        std::vector<size_t> permutations_;
        //! Independent elements d_{N,d} and t_{N,d} 
        t_Independents independents_;
        //! Number of sites in the unit lattice-cell.
        size_t nsites_;
        //! Total number of sites in the supercell.
        size_t card_;
        //! is trivial.
        bool is_trivial_;
    };

    boost::shared_ptr< std::vector<Transform> >  create_transforms( Crystal::Lattice const &_lat );
  }
}

#endif

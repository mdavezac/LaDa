//
//  Version: $Id$
//
#ifndef LADA_ENUM_TRANSLATION_H_
#define LADA_ENUM_TRANSLATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>

#include <vector>

#include "numeric_type.h"

namespace LaDa
{
  namespace enumeration
  {
    //! \brief Transform an integer structure to another. 
    //! \see Appendix of PRB 80, 014120.
    class Translation
    {
      public:
        //! Constructor
        Translation   (atat::iVector3d const& _trans, atat::iVector3d const &_smith, size_t _nsites)
                  : smith_(_smith), nsites_(_nsites),
                    card_(_nsites*_smith(0)*_smith(1)*_smith(2)), 
                    translation_(_trans) {}
                   
        //! Copy Constructor
        Translation   (Translation const &_c) 
                    : smith_(_c.smith_), nsites_(_c.nsites_),
                      card_(_c.card_), translation_(_c.translation_)  {}

        //! Returns the transformed integer.
        t_uint operator()(t_uint _x, FlavorBase const &_flavorbase) const;

        //! Comparison operator.
        bool operator==(Translation const& _c) const
        { return smith_ == _c.smith_ and nsites_ == _c.nsites_ and translation_ == _c.translation_; }

      private:
        //! Smith id.
        atat::iVector3d smith_;
        //! Number of sites.
        size_t nsites_;
        //! Cardinality.
        size_t card_;
        //! Translation within the Smith normal form.
        atat::iVector3d translation_;
    };

    boost::shared_ptr< std::vector<Translation> > 
      create_translations(atat::iVector3d const &_smith, size_t _nsites);
  }
}

#endif

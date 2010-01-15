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
        Translation   (Eigen::Vector3i const &_smith, size_t _nsites);
        //! Copy Constructor
        Translation   (Translation const &_c) 
                    : card_(_c.card_), permutations_(_c.permutations_),
                      i_trans_(_c.i_trans_), i_trans_end_(_c.i_trans_end_)  {}

        //! Returns the transformed integer.
        t_uint operator()(t_uint _x, FlavorBase const &_flavorbase) const;
        //! Iterates to next translation.
        bool operator++() 
        {
          ++i_trans_;
          if( i_trans_ != i_trans_end_ ) return true;
          i_trans_ = permutations_.rbegin();
          return false;
        }
        //! Returns number of translations.
        size_t size() const { return permutations_.size(); }
        
      private:
        //! Cardinality.
        size_t card_;
        //! Permutations.
        std::vector< std::vector<size_t> > permutations_;
        //! Iterator over permutation.
        std::vector< std::vector<size_t> > :: const_reverse_iterator i_trans_;
        //! Iterator over permutation.                      
        std::vector< std::vector<size_t> > :: const_reverse_iterator i_trans_end_;
    };
  }
}

#endif

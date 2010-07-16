#ifndef LADA_ENUM_LABEL_EXCHANGE_H_
#define LADA_ENUM_LABEL_EXCHANGE_H_

#include "LaDaConfig.h"

#include <iostream>

#include <vector>

#include "numeric_type.h"

namespace LaDa
{
  namespace enumeration
  {
    //! \brief Transform an integer structure to another. 
    //! \see Appendix of PRB 80, 014120.
    class LabelExchange
    {
      public:
        //! Constructor
        LabelExchange   (size_t _card, size_t _nflavors)
                      : card_(_card)
          { for(size_t i(0); i < _nflavors; ++i) permutation_.push_back(i); }
                   
        //! Copy Constructor
        LabelExchange   (LabelExchange const &_c) 
                      : card_(_c.card_), permutation_(_c.permutation_) {}

        //! Returns the transformed integer.
        t_uint operator()(t_uint _x, FlavorBase const &_flavorbase) const;
        //! Performs next permutation.
        bool operator++()
          { return std::next_permutation(permutation_.begin(), permutation_.end()); }

      private:
        //! Cardinality.
        size_t const card_;
        //! LabelExchange within the Smith normal form.
        std::vector<size_t> permutation_;
    };
  }
}

#endif

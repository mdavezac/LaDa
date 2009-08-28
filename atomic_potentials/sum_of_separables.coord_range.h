//
//  Version: $Id$
//
#ifdef LADA_WITH_CONST
#  define LADA_CONST_(a) const_ ## a
#  define LADA_CONST const 
#else
#  define LADA_CONST_(a) a
#  define LADA_CONST
#endif

 //! Inverts iteration, such that it goes over variables first, and ranks second.
 class LADA_CONST_(coord_range)
 {
#  ifdef LADA_VMRASSERT
#    error LADA_VMRASSERT already defined.
#  endif
#  ifdef LADA_DEBUG
#    define LADA_VMRASSERT LADA_ASSERT( index_ < max_index_, "Range out of range.\n") 
#  else
#    define LADA_VMRASSERT 
#  endif
   public:
     //! Type of the dereferenced value.
     typedef size_t value_type;

     
#    ifndef LADA_WITH_CONST
       //! Iterator over ranks.
       class rank_range;
#      include "coord_range.rank_range.h"

#      define LADA_WITH_CONST
#      undef LADA_CONST_
#      undef LADA_CONST
#      define LADA_CONST_(a) const_ ## a
#      define LADA_CONST const 
       //! Iterator over ranks.
       class const_rank_range;
#      include "coord_range.rank_range.h"
#      undef LADA_CONST_
#      undef LADA_CONST
#      define LADA_CONST_(a) a
#      define LADA_CONST
#      undef LADA_WITH_CONST
#    else 
       //! Iterator over ranks.
       class const_rank_range;
#      include "coord_range.rank_range.h"
#    endif
     //! Constructor.
     LADA_CONST_(coord_range)   (SumOfSeparables LADA_CONST &_sumofseps )
                              : sumofseps_(_sumofseps), index_(0),
                                max_index_(_sumofseps.nb_coordinates() ) {}
     //! Copy Constructor.
     LADA_CONST_(coord_range)   ( LADA_CONST_(coord_range) const &_c)
                              : sumofseps_(_c.sumofseps_),
                                index_(_c.index_), max_index_(_c.max_index_) {}

#    ifndef LADA_WITH_CONST
       //! Range over ranks.
       rank_range range()
       {
         LADA_VMRASSERT; 
         return rank_range( sumofseps_.begin(), sumofseps_.end(), index_ ); 
       }
#    endif
 
     //! Range over ranks.
     const_rank_range range() const
     {
       LADA_VMRASSERT; 
       return const_rank_range( sumofseps_.begin(), sumofseps_.end(), index_ ); 
     }

     //! Increments.
     bool operator++() { ++index_; return this->operator bool(); }

     //! Continue iterating or not.
     operator bool() const { return index_ < max_index_; }

     //! Returns size of the range.
     size_t size() const { return max_index_; }

     value_type operator*() const { return index_; }

   private:
     //! Type of the numeric values.
     SumOfSeparables LADA_CONST &sumofseps_;
     //! Current index.
     size_t index_;
     //! Maximum index.
     size_t max_index_;

 };

#undef LADA_VMRASSERT
#undef LADA_WITH_CONST
#undef LADA_CONST
#undef LADA_CONST_

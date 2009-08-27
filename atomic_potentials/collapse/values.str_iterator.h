//
//  Version: $Id$
//
#ifdef LADA_CONST_
#  error LADA_CONST_ already defined
#endif
#ifdef LADA_WITH_CONST
#  define LADA_CONST_(a) const_ ## a
#else
#  define LADA_CONST_(a) a
#endif

class LADA_CONST_(str_iterator)
{
    friend class Values;
  public:
    //! Iterator over representations.
    class LADA_CONST_(rep_iterator);
#   include "values.rep_iterator.h"

    //! Constructor.
    LADA_CONST_(str_iterator)   (size_t _i, t_CoordRankValues const &_cr)
                              : i_(_i), n_(0), coord_rank_values_(_cr) {}
    //! Copy Constructor.
    LADA_CONST_(str_iterator)   (LADA_CONST_(str_iterator) const &_c) 
                              : i_(_c.i_), n_(_c.n_), coord_rank_values_(_c.coord_rank_values_),
                                i_coord_rank_value_(_c.i_coord_rank_value_),
                                i_rank_value_(_c.i_rank_value_),
                                i_function_values_(_c.i_function_values_) {}

    //! Increments iterator. 
    void operator++() 
      { ++i_coord_rank_value_; ++i_rank_value_; ++i_function_values_; ++n_; }
    //! Decrements iterator.
    void operator--() 
      { --i_coord_rank_value_; --i_rank_value_; --i_function_values_; --n_; }
    
    //! Returns iterator over representations.
    LADA_CONST_(rep_iterator) begin() const
    {
      LADA_CONST_(rep_iterator) result(i_, n_, coord_rank_values_);
      result.i_coord_rank_value_ = i_coord_rank_value_->begin();
      result.i_rank_value_ = i_rank_value_->begin();
      result.i_function_values_ = i_function_values_->begin();
      return result;
    }
    //! Returns iterator over representations.
    LADA_CONST_(rep_iterator) end() const
    {
      LADA_CONST_(rep_iterator) result(i_, n_, coord_rank_values_);
      result.i_coord_rank_value_ = i_coord_rank_value_->end();
      // result.i_rank_value_ = i_rank_value_->end();
      // result.i_function_values_ = function_values_->end();
      return result;
    }

    //  True if iterators are at same position.
    bool operator==(LADA_CONST_(str_iterator) const& _b) const
    {
      LADA_ASSERT( i_ == _b.i_, "Inequivalent iterators.\n");
      return i_coord_rank_value_ != _b.i_coord_rank_value_; 
    }
    //! False if iterators are at same position.
    bool operator!=(LADA_CONST_(str_iterator) const& _b) const { return not operator==(_b); }

  protected:
    //! Current variable.
    size_t i_;
    //! Current variable.
    size_t n_;
    //! Values aggregate per coordinate and rank.
    t_CoordRankValues const& coord_rank_values_;
    //! Iterator over other factors.
    t_CoordRankValues::value_type::LADA_CONST_(iterator) i_coord_rank_value_;
    //! Iterator over other factors.
    t_RankValues::LADA_CONST_(iterator) i_rank_value_;
    //! Iterator over function values.
    t_FunctionValues::value_type::LADA_CONST_(iterator) i_function_values_;
};

# undef LADA_CONST_
# undef LADA_WITH_CONST

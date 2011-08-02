class LADA_CONST_(rank_iterator)
{
    friend class LADA_CONST_(rep_iterator);
  public:
    //! Iterator over function values.
    typedef t_FunctionValues::value_type::value_type::value_type
                            ::value_type::LADA_CONST_(iterator) function_iterator;
    //! Constructor.
    LADA_CONST_(rank_iterator)   (size_t _i, size_t _nstr, size_t _nrep,
                                  t_CoordRankValues const& _cr)
                               : i_(_i), n_str_(_nstr), n_rep_(_nrep), n_(0),
                                 coord_rank_values_(_cr) {}
    //! Copy Constructor.
    LADA_CONST_(rank_iterator)   (LADA_CONST_(rank_iterator) const &_c) 
                               : i_(_c.i_), n_str_(_c.n_str_), n_rep_(_c.n_rep_), n_(_c.n_),
                                 coord_rank_values_(_c.coord_rank_values_),
                                 i_coord_rank_value_(_c.i_coord_rank_value_),
                                 i_rank_value_(_c.i_rank_value_),
                                 i_function_values_(_c.i_function_values_) {}

    //! Returns factors (i!=j) for alternating least-square-fit.
    numeric_type other() const
    {
      numeric_type const coord_rank_value( *i_coord_rank_value_ );
      if( math::is_null(coord_rank_value) ) 
      {
        numeric_type result(1);
        t_CoordRankValues::const_iterator i_cr = coord_rank_values_.begin();
        t_CoordRankValues::const_iterator const i_cr_end = coord_rank_values_.end();
        for( size_t i(0); i_cr != i_cr_end; ++i_cr, ++i)
          if( i != i_ ) result *= (*i_cr)[n_str_][n_rep_][n_];
        return result;
      }
      return (*i_rank_value_) / coord_rank_value;
    }
    //! Returns rank values.
    numeric_type all() const { return *i_rank_value_; }
    //! Returns iterator to function values.
    function_iterator begin() const { return i_function_values_->begin(); }
    //! Returns iterator to function values.
    function_iterator end() const { return i_function_values_->end(); }
    //! Increments operator.
    void operator++() 
      { ++i_coord_rank_value_; ++i_rank_value_; ++i_function_values_; ++n_; }
    //! Decrements operator.
    void operator--() 
      { --i_coord_rank_value_; --i_rank_value_; --i_function_values_; --n_; }


   //  True if iterators are at same position.
   bool operator==(LADA_CONST_(rank_iterator) const& _b ) const
     { return i_coord_rank_value_ == _b.i_coord_rank_value_;  }
   //! False if iterators are at same position.
   bool operator!=(LADA_CONST_(rank_iterator) const& _b ) const { return not operator==(_b); }

  protected:
    //! Current variable.
    size_t i_;
    //! Current structure.
    size_t n_str_;
    //! Current representation.
    size_t n_rep_;
    //! Current rank
    size_t n_;
    //! Values aggregate per coordinate and rank.
    t_CoordRankValues const& coord_rank_values_;
    //! Iterator over other factors.
    t_CoordRankValues::value_type::value_type::value_type::LADA_CONST_(iterator) i_coord_rank_value_;
    //! Iterator over other factors.
    t_RankValues::value_type::value_type::LADA_CONST_(iterator) i_rank_value_;
    //! Iterator over function values.
    t_FunctionValues::value_type::value_type::value_type::LADA_CONST_(iterator) i_function_values_;
};

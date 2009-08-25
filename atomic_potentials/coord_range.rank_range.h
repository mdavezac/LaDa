class LADA_CONST_(rank_range)
{
# ifdef LADA_WITH_CONST
    friend class const_coord_range;
# endif
  friend class coord_range;
  //! Type of the var iterator.
  typedef SumOfSeparables::t_Function::const_iterator const_var_iterator;
  //! Type of the var iterator.
  typedef SumOfSeparables::t_Function::iterator var_iterator;
  public:
    //! value type
    typedef SumOfSeparables::t_Function::LADA_CONST_(iterator)::value_type value_type;
    //! reference type
    typedef SumOfSeparables::t_Function::LADA_CONST_(iterator)::reference reference;
    //! pointer type
    typedef SumOfSeparables::t_Function::LADA_CONST_(iterator)::pointer pointer;
    
    //! Iterator over functions and coefficients for set rank and variable.
    typedef SumOfSeparables::t_Function::t_Function::const_iterator const_iterator;
# ifndef LADA_WITH_CONST
    //! Iterator over functions and coefficients for set rank and variable.
    typedef SumOfSeparables::t_Function::t_Function::iterator iterator;
# endif

    //! Copy Constructor.
    LADA_CONST_(rank_range)   (LADA_CONST_(rank_range) const &_c)
                            : i_rank_(_c.i_rank_), i_rank_end_(_c.i_rank_end_),\
                              index_(_c.index_) {}

#   ifndef LADA_WITH_CONST
      //! Iterator over ranks.
      iterator begin() { return get_i_var()->begin(); }
      //! Iterator over ranks.
      iterator end() { return get_i_var()->end(); }
      //! Iterator over ranks.
      const_iterator cbegin() const { return const_get_i_var()->begin(); }
      //! Iterator over ranks.
      const_iterator cend() const { return const_get_i_var()->end(); }
#   else
      //! Iterator over ranks.
      const_iterator begin() const { return const_get_i_var()->begin(); }
      //! Iterator over ranks.
      const_iterator end() const { return const_get_i_var()->end(); }
#   endif 

    //! Increments.
    void operator++()
    {
      if( i_rank_ != i_rank_end_ ) return;
      do 
      { 
        ++i_rank_; 
        if( i_rank_ == i_rank_end_ ) break; 
      } while( i_rank_->get<0>().size() <= index_  ); 
    }

    //! Deref.
    value_type operator*() const { return LADA_CONST_(get_i_var)().operator*(); }
    //! Deref.
    pointer operator->() const { return LADA_CONST_(get_i_var)().operator->(); }
    //! True if still iteratable.
    operator bool() const { return i_rank_ != i_rank_end_; }

  private:
    //! Constructor.
    LADA_CONST_(rank_range)   (SumOfSeparables::LADA_CONST_(iterator) const &_irank,
                               SumOfSeparables::LADA_CONST_(iterator) const &_irank_end,
                               size_t _i)
                            : i_rank_(_irank), i_rank_end_(_irank_end), index_(_i) 
    {
      try
      {
        while( i_rank_->get<0>().size() <= index_  )
        {
          ++i_rank_;
          if( i_rank_ == i_rank_end_ ) break;
        } 
      }
      catch(...)
      {
        i_rank_ == i_rank_end_;
      }
    }

    //! Gets correct variable iterator.
    const_var_iterator const_get_i_var() const
    { 
      LADA_ASSERT( i_rank_ != i_rank_end_, "Range not iteratable.\n" ) 
      const_var_iterator result( i_rank_->get<0>().begin() ); 
      for(size_t i(0); i < index_; ++i, ++result); 
      return result; 
    };

#   ifndef LADA_WITH_CONST
      //! Gets correct variable iterator.
      var_iterator get_i_var() const
      { 
        LADA_ASSERT( i_rank_ != i_rank_end_, "Range not iteratable.\n" ) 
        var_iterator result( i_rank_->get<0>().begin() ); 
        for(size_t i(0); i < index_; ++i, ++result); 
        return result; 
      };
#   endif

    //! Current rank positions.
    SumOfSeparables::LADA_CONST_(iterator) i_rank_;
    //! Last rank position.
    SumOfSeparables::LADA_CONST_(iterator) i_rank_end_;
    //! Variable index.
    size_t index_;
};



#ifndef LADA_ITERATOR_NAME
#  error LADA_ITERATOR_NAME must be defined.
#endif
#ifndef LADA_WITH_CONST
#  error LADA_WITH_CONST must be defined.
#endif
class LADA_ITERATOR_NAME
{
  public:
    //! Dereferenced value.
    class reference
    {
      friend class LADA_ITERATOR_NAME;
      public:
        //! Constructor.
        reference() {}
        //! Constructor.
        reference   ( t_Coefficients::LADA_ITERATOR_NAME const& _c,
                      t_Functions::LADA_ITERATOR_NAME const& _f )
                  : i_coef_(_c), i_func_(_f) {}
        //! CopyConstructor.
        reference   ( reference const &_c)
                  : i_coef_(_c.i_coef_), i_func_(_c.i_func_) {}
        //! Returns a reference to that component.
        numeric_type LADA_WITH_CONST& operator[](size_t _i) const
        {
          LADA_ASSERT( _i < Functions::N, "Index out of range." )
          return *(i_coef_+_i);
        }
        //! Reference to the function.
        t_Functions::value_type LADA_WITH_CONST& function() const
          { return *i_func_; }
        //! Calls function with coefficient == 1.
        Functions::result_type operator()( Functions::arg_type::first_type _a ) const
          { return (*i_func_)(_a);  }
        //! Calls function.
        Functions::result_type operator()( Functions::arg_type const &_a ) const
          { return (*i_func_)(_a.first) * operator[](_a.second); }

      private:
        //! Current position.
        t_Coefficients::LADA_ITERATOR_NAME i_coef_;
        //! Current position.
        t_Functions::LADA_ITERATOR_NAME i_func_;
    };

    //! Pointer.
    typedef reference* pointer;
    //! Null value type.
    typedef void value_type;

    //! Constructor.
    LADA_ITERATOR_NAME   (t_Functions::LADA_ITERATOR_NAME const &_f, 
                          t_Coefficients::LADA_ITERATOR_NAME const &_c)
                       : i_func_(_f), i_coef_(_c), is_set_(false) {};
    //! Copy Constructor.
    LADA_ITERATOR_NAME   (LADA_ITERATOR_NAME const& _c)
             : i_func_(_c.i_func_), i_coef_(_c.i_coef_), 
               is_set_(_c.is_set_), value_(_c.value_) {};

    //! Increments iterator.
    LADA_ITERATOR_NAME operator++(int) { return operator++(); }
    //! Increments iterator.
    LADA_ITERATOR_NAME& operator++()
      { is_set_ = false; ++i_func_; i_coef_ += Functions::N; return *this; }
    //! Decrements iterator.
    LADA_ITERATOR_NAME operator--(int) { return operator--(); }
    //! Decrements iterator.
    LADA_ITERATOR_NAME& operator--()
      { is_set_ = false; --i_func_; i_coef_ -= Functions::N; return *this; }
    //! Dereference iterator.
    reference operator*() const { set_(); return value_; }
    //! Dereference iterator.
    pointer operator->() const { set_(); return &value_; }
    //! Returns true if iterator are at same position.
    bool operator==(LADA_ITERATOR_NAME const& _b) const
      { return i_func_ == _b.i_func_; }
    //! Returns false if iterator are at same position.
    bool operator!=(LADA_ITERATOR_NAME const& _b) const
      { return i_func_ != _b.i_func_; }

  private:
    //! Sets the value.
    void set_() const
    {
      if(is_set_) return;
      value_.i_coef_ = i_coef_;
      value_.i_func_ = i_func_;
      is_set_ = true;
    }
    //! Current function position.
    t_Functions::LADA_ITERATOR_NAME i_func_;
    //! Current coefficient position.
    t_Coefficients::LADA_ITERATOR_NAME i_coef_;
    //! Whether the value has been set.
    mutable bool is_set_;
    //! A tuple for pointer stuff.
    mutable reference value_;
};

#undef LADA_ITERATOR_NAME
#undef LADA_WITH_CONST

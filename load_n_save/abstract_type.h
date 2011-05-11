
template< class T1, class T2 > struct cloner;

template<>
  struct cloner<void, void> 
  {
    virtual Base* clone( Base* _a1, Base *_a2 )
    {
      _a1.clone( _a2 );
    }
    virtual Base* clone( Base *_a2 );
    virtual Base* clone();
    template< class T2 >
      void another( T2 *_a2 ) 
      {
      };
  } 
template< class T1 >
  struct cloner<T1, void> : public cloner<void, void>
  {
    virtual Base* clone( Base *_a2 )
    {
      return a2.clone<T1>(a1_);
    }
  };

class Base
{
  virtual Base* clone( Base *_a1, Base* _a2 ) const
  {
    _a1->clone( _a2 );
  }
  virtual Base* clone( Base* _a2 ) const;
  template< class T1 > 
    Base* self( T1* _a1 ) const;
    {
      class dd : public T1
      {
      };
    }
}

template< class T1 > class Derived
  {
    virtual Base* clone( Base* _a2 ) const;
    {
      a2->self< Derived<T1> >( this );
    }
  }

template< class T1, class T2>
  Base* cloner( Base* _a1, Base* _a2 )
  {
    return new Double<T1, T2>( _a1, _a2 );
  }


template< class T1 >
  Base* cloner_left( Base* _a1, Base* _a2 )
  {
    return _a2.clone( _a1, &cloner<T1> )
  }

template< class T2 >
  Base* cloner_right( Base* _a1, Base* _a2 )
  {
    return _a1.clone( _a2, &cloner<T2> )
  }

#ifndef _DARWIN_FUNCTORS_H_
#define _DARWIN_FUNCTORS_H_

#include <eo/eoOp.h>
#include <eo/eoContinue.h>
#include <eo/utils/eoRNG.h>
#include <eo/utils/eoState.h>

#include <opt/types.h>
#include <lamarck/structure.h>

namespace darwin 
{
  template <class A1, class R>
 class const_eoUF : public eoFunctorBase, public std::unary_function<A1, R>
 {
   public:
     typedef A1 t_Argument;
     typedef R  t_Return;

   public:
     virtual ~const_eoUF() {}
   
     virtual t_Return operator()( t_Argument ) const = 0;
   
     static eoFunctorBase::unary_function_tag functor_category()
       { return eoFunctorBase::unary_function_tag(); }
 };

  // generic class for converting member function to binary genetic operators
  template<class T_CLASS, class T_OBJECT, class T_ARG >
  class mem_binop_t : public eoBinOp<T_OBJECT> 
  {
    public:
      typedef T_CLASS t_Class;
      typedef T_OBJECT t_Object;
      typedef T_ARG t_Arg; 
      typedef bool ( t_Class::*t_Function )(t_Object &, const t_Object&, t_Arg);

    private:
      t_Class &class_obj;
      t_Arg arg;
      t_Function class_func;
      std::string class_name;

    public:
      explicit
        mem_binop_t   ( t_Class &_co, t_Function _func,
                        const std::string &_cn, t_Arg _arg )
                    : class_obj(_co), arg(_arg), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Object &_object1, const t_Object &_object2) 
        {  return ( (class_obj.*class_func) )( _object1, _object2, arg ); }

  }; 

  // generic class for converting member function to binary operators
  template<class T_CLASS, class T_OBJECT>
  class mem_binop_t< T_CLASS, T_OBJECT, void > : public eoBinOp<T_OBJECT> 
  {
    public:
      typedef T_CLASS t_Class;
      typedef T_OBJECT t_Object;
      typedef bool ( t_Class::*t_Function )(t_Object &, const t_Object&);

    private:
      t_Class &class_obj;
      bool ( t_Class::*class_func )(t_Object &, const t_Object&);
      std::string class_name;

    public:
      explicit
        mem_binop_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Object &_object1, const t_Object &_object2) 
        {  return ( (class_obj.*class_func) )( _object1, _object2 ); }

  }; 

  template<class T_INDIVIDUAL>
  class ObjectOp_To_BinGenOp : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef typename t_Individual::t_Object t_Object;

    private:
      eoBinOp<t_Object> &class_object;

    public:
      explicit
        ObjectOp_To_BinGenOp ( eoBinOp<t_Object> &_class ) : class_object(_class) {};

      unsigned max_production(void) { return 1; } 
   
      void apply(eoPopulator<t_Individual>& _pop)
      {
        t_Individual& a = *_pop;
        const t_Individual& b = _pop.select();
  
        if ( class_object( (t_Object&) a, (const t_Object&)b ) )
          a.invalidate();
      }
      virtual std::string className() const {return class_object.className();}

  }; 

  // generic class for converting member function to binary operators
  template<class T_CLASS, class T_OBJECT>
  class mem_monop_t : public eoMonOp<T_OBJECT> 
  {
    public:
      typedef T_CLASS t_Class; 
      typedef T_OBJECT t_Object;
      typedef bool ( t_Class::*t_Function )(t_Object &);

    private:
      t_Class &class_obj;
      t_Function class_func;
      std::string class_name;

    public:
      explicit
        mem_monop_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Object &_object) 
        {  return ( (class_obj.*class_func) )( _object); }

  }; 
  template<class T_CLASS, class T_OBJECT>
  class const_mem_monop_t : public eoMonOp<const T_OBJECT> 
  {
    public:
      typedef T_CLASS t_Class; 
      typedef T_OBJECT t_Object;
      typedef bool ( t_Class::*t_Function )(const t_Object &);

    private:
      t_Class &class_obj;
      t_Function class_func;
      std::string class_name;

    public:
      explicit
        const_mem_monop_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (const t_Object &_object) 
        {  return ( (class_obj.*class_func) )( _object); }

  }; 
  // generic class for converting member function to monary operators
  template<class T_CLASS, class T_INDIVIDUAL>
  class mem_monop_indiv_t : public eoMonOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_CLASS t_Class; 
      typedef T_INDIVIDUAL t_Individual;
      typedef typename t_Individual::t_Object t_Object;
      typedef bool ( t_Class::*t_Function )(t_Object &);

    private:
      t_Class &class_obj;
      t_Function class_func;
      std::string class_name;

    public:
      explicit
        mem_monop_indiv_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                          : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Individual &_indiv) 
        { return ( (class_obj.*class_func) )( (t_Object &) _indiv ); }

  }; 
  // generic class for converting member function to zero operators
  template<class T_CLASS, class T_OBJECT>
  class mem_zerop_t : public eoF<bool>
  {
    public:
      typedef T_CLASS t_Class; 
      typedef bool ( t_Class::*t_Function )();

    private:
      t_Class &class_obj;
      t_Function class_func;
      std::string class_name;

    public:
      explicit
        mem_zerop_t   ( t_Class &_co, t_Function _func, const std::string &_cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() () 
        {  return ( (class_obj.*class_func) )(); }

  }; 
  template<class T_INDIVIDUAL>
  class ObjectOp_To_MonGenOp : public eoGenOp<T_INDIVIDUAL> 
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef typename t_Individual::t_Object t_Object;

    private:
      eoMonOp<t_Object> &class_object;

    public:
      explicit
        ObjectOp_To_MonGenOp ( eoMonOp<t_Object> &_class ) : class_object(_class) {};

      unsigned max_production(void) { return 1; } 
   
      void apply(eoPopulator<t_Individual>& _pop)
      {
        if ( class_object( (t_Object &)(*_pop) ) )
          (*_pop).invalidate();
      }
      virtual std::string className() const {return class_object.className();}
  }; 

  // translates an eoOp<t_Object> into a genetic operator
  // eoGenOp<t_Individual>, where trait t_Individual::t_Object should
  // be defined. Note that _class is stored at this point for
  // convenience (see note in darwin.impl.h to this effect) 
  template<class T_INDIVIDUAL>
  eoGenOp<T_INDIVIDUAL>* Wrap_ObjectOp_To_GenOp( eoOp<typename T_INDIVIDUAL::t_Object> *_class, eoState &_state)
  {
    typedef T_INDIVIDUAL t_Individual;
    typedef typename t_Individual::t_Object t_Object;
    switch(_class->getType())
    {
      case eoOp<t_Object>::unary: 
        _state.storeFunctor( (eoMonOp<t_Object>*)(_class) );
        return & _state.storeFunctor(
            new ObjectOp_To_MonGenOp<t_Individual>( *((eoMonOp<t_Object>*)_class)) );
      case eoOp<T_INDIVIDUAL>::binary:
        _state.storeFunctor( (eoBinOp<t_Object>*)(_class) );
        return & _state.storeFunctor(
            new ObjectOp_To_BinGenOp<t_Individual>(*( (eoBinOp<t_Object>*)_class) ) );
      case eoOp<T_INDIVIDUAL>::quadratic:
      case eoOp<T_INDIVIDUAL>::general:
        std::cerr << "No Implementation for quadratic and general operators "
                  << "in Wrap_ObjectOp_To_GenOp yet!"
                  << std::endl;
        exit(0);
    }
    return NULL; // just to please gcc
  }

  template< class T_OBJECT, class T_CONTAINER = typename T_OBJECT::t_Container >
  class Crossover : public eoBinOp<T_OBJECT>
  {
    protected:
      typedef T_OBJECT t_Object;
      typedef T_CONTAINER t_Container;

    public:
      types::t_real probability;

    public:
      Crossover(types::t_real _prob=0.0) : probability(_prob) {}
      ~Crossover() {}

      bool operator() ( t_Object &_obj1, const t_Object &_obj2 )
      {
        typename t_Container :: iterator i_var1 = _obj1.begin();
        typename t_Container :: const_iterator i_var2 = _obj2.begin();
        typename t_Container :: const_iterator i_var2_end = _obj2.end();
        for(; i_var2 != i_var2_end; ++i_var1, ++i_var2)
          if ( rng.uniform() < probability ) 
            *i_var1 = *i_var2;
        return true;
      }
  };
  template<> 
  class Crossover<Ising_CE::Structure, Ising_CE::Structure::t_Atoms> : public eoBinOp<Ising_CE::Structure>
  {
    protected:
      typedef Ising_CE::Structure t_Object;
      typedef Ising_CE::Structure::t_Atoms t_Container;

    public:
      types::t_real probability;

    public:
      Crossover(types::t_real _prob=0.0) : probability(_prob) {}
      ~Crossover() {}

      bool operator() ( t_Object &_obj1, const t_Object &_obj2 )
      {
        t_Container :: iterator i_var1 = _obj1.atoms.begin();
        t_Container :: const_iterator i_var2 = _obj2.atoms.begin();
        t_Container :: const_iterator i_var2_end = _obj2.atoms.end();
        for(; i_var2 != i_var2_end; ++i_var1, ++i_var2)
          if ( rng.uniform() < probability ) 
            i_var1->type = i_var2->type;
        return true;
      }
  };


  // a dummy operator which does nothing 
  template< class T_OBJECT >
  class DummyOp : public eoMonOp<T_OBJECT>
  {
    protected:
      typedef T_OBJECT t_Object;

    public:
      DummyOp() {};
      ~DummyOp() {}

      bool operator() ( t_Object &_obj1 )
      { return false; } // do nothing!
  };

  template< class T_INDIVIDUAL  >
  class Continuator : public eoContinue< T_INDIVIDUAL >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;

    protected:
      eoF<bool> &op;

    public:
      Continuator( eoF<bool> &_op ) : op(_op) {};
      ~Continuator() {}

      bool operator()(const eoPop<t_Individual> &_pop )
      {
        return op();
      }
  };


}
#endif // _DARWIN_FUNCTORS_H_

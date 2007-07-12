#ifndef _GENERIC_OPS_H_
#define _GENERIC_OPS_H_

#include <eo/eoOp.h>


namespace LaDa 
{
  // generic class for converting member function to binary genetic operators
  template<class _Tp, class t_Object, class t_Arg >
  class mem_binop_t : public eoBinOp<t_Object> 
  {
    private:
      _Tp &class_obj;
      t_Arg arg;
      bool ( _Tp::*class_func )(t_Object &, const t_Object&, t_Arg);
      std::string class_name;

    public:
      explicit
        mem_binop_t   ( _Tp &_co, bool (_Tp::*_func)(t_Object&, const t_Object&, t_Arg),
                        std::string _cn, t_Arg _arg )
                    : class_obj(_co), arg(_arg), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Object &_object1, const t_Object &_object2) 
        {  return ( (class_obj.*class_func) )( _object1, _object2, arg ); }

  }; 

  // generic class for converting member function to binary genetic operators
  template<class _Tp, class t_Object>
  class mem_binop_t< _Tp, t_Object, void > : public eoBinOp<t_Object> 
  {
    private:
      _Tp &class_obj;
      bool ( _Tp::*class_func )(t_Object &, const t_Object&);
      std::string class_name;

    public:
      explicit
        mem_binop_t   ( _Tp &_co, bool (_Tp::*_func)(t_Object&, const t_Object&),
                            std::string _cn )
                    : class_obj(_co), class_func(_func), class_name(_cn) {};

      void set_className( std::string &_cn) { class_name = _cn; }
      virtual std::string className() const { return class_name; }

      bool operator() (t_Object &_object1, const t_Object &_object2) 
        {  return ( (class_obj.*class_func) )( _object1, _object2 ); }

  }; 

} // endif LaDa

#endif

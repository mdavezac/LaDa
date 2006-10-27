#ifndef _GENCOUNT_H_
#define _GENCOUNT_H_

namespace LaDa 
{
  // generation counter
  class GenCount
  {
    protected:
      t_unsigned age;
    public:
      GenCount( const GenCount &_gc) : age(_gc.age) {};
      GenCount( t_unsigned _age) : age(_age) {};
      GenCount() : age(0) {};
      void operator ++() 
        { ++age; }
      const t_unsigned operator()() const
        { return age; }
  };
} // namespace LaDa

#endif

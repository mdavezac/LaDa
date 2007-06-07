#ifndef _INDIVIDUAL_H_
#define _INDIVIDUAL_H_

#include <eo/eoScalarFitness.h>
#include <opt/opt_function_base.h>
#include <exception>
#include <iostream>
#include <stdexcept>       // std::runtime_error
#include <eo/eoObject.h>      // eoObject
#include <eo/eoPersistent.h>  // eoPersistent
#include <eo/eoOpContainer.h>  

#include<algorithm>
#include<functional>

#include <opt/types.h>
using namespace types;

// functional should be some class derived from 
// opt::Function<t_Type, t_Container>, where t_Container supports random
// access iterators
//

namespace LaDa 
{
  template<class FITNESS = eoMinimizingFitness, class FUNCTIONAL = opt::Function<> >
  class Individual : public FUNCTIONAL, public eoObject, public eoPersistent
  {

    public:
      typedef FUNCTIONAL t_Functional;
      typedef typename t_Functional::t_Type t_Type;
      typedef typename t_Functional::t_Container t_Container;
      typedef typename t_Functional::iterator iterator;
      typedef typename t_Functional::const_iterator const_iterator;
      typedef t_Type AtomType;
      typedef FITNESS t_Fitness;
      typedef FITNESS Fitness; // for eo

    public:
      using t_Functional :: variables;
      using t_Functional :: set_variables;

    private:
      t_Fitness repFitness;
      t_int MotU_ref; // ref to master of the universe
      static bool baseline_is_valid;
      bool quantity_is_valid;
      t_Type quantity, baseline;
      t_unsigned age;
      t_Container *phenotype;

    public:
      static bool is_using_phenotype;
      

    public: 
      Individual() : MotU_ref(0), quantity_is_valid (false), age(0)
      {
        variables = new t_Container;
        phenotype = NULL;
        if ( is_using_phenotype )
          phenotype = new t_Container;
        else 
          phenotype = variables;
      };
      
      ~Individual()
      {
        if (variables)
          delete variables;  
        variables = NULL; 
        if (is_using_phenotype)
          delete phenotype;  
        phenotype = NULL;
      };

      Individual(const Individual<FITNESS, t_Functional> &_indiv ) :
              repFitness(_indiv.repFitness),
              MotU_ref(_indiv.MotU_ref), 
              quantity_is_valid (_indiv.quantity_is_valid),
              quantity(_indiv.quantity), baseline(_indiv.baseline), 
              age(_indiv.age), phenotype(NULL)
      {
        variables = new t_Container;
        *variables = *_indiv.variables;
        typename t_Container :: iterator i_var = variables->begin();
        typename t_Container :: iterator i_end = variables->end();
        for(; i_var != i_end; ++i_var )
          *i_var = ( *i_var > 0.0 ) ? 1.0 : -1.0;
        phenotype = variables;
        if ( is_using_phenotype )
        {
          phenotype = new t_Container;
          *phenotype = *_indiv.phenotype;
        }
      };

      virtual t_Container* get_variables() const
      { return phenotype; }
      virtual t_Type get_concentration() const
      { return ( std::accumulate(phenotype->begin(), phenotype->end(), 0.0) 
                  / ( (t_real) phenotype->size() ) ); }

      void invalidate() { quantity_is_valid = false; }
      static void invalidate_baseline();
      static void validate_baseline();
      static bool is_baseline_valid();
      bool invalid() const
        { return not quantity_is_valid; }


      void set_quantity( t_Type _q )
      {
        quantity = _q;
        quantity_is_valid = true;
      }
      t_Type get_quantity() const
        { return quantity; }
      void set_baseline( t_Type _b )
        { baseline = _b; }
      t_Type get_baseline() const
        { return baseline; }

      void  set_age( t_unsigned _age )
        { age = _age; }
      t_unsigned  get_age() const
        { return age; }
      void  set_MotU_ref( t_int _m )
        { MotU_ref = _m; }
      t_int get_MotU_ref() const
        { return MotU_ref; }

      bool operator<(const Individual<FITNESS, t_Functional>& _eo2) const
        { return fitness() < _eo2.fitness(); }
      bool operator>(const Individual<FITNESS, t_Functional>& _eo2) const
        { return !(fitness() <= _eo2.fitness()); }
      void operator=( const Individual<FITNESS, t_Functional> &_indiv )
      {
        repFitness = _indiv.repFitness;
        MotU_ref   = _indiv.MotU_ref;
        quantity_is_valid = _indiv.quantity_is_valid;
        quantity = _indiv.quantity;
        baseline = _indiv.baseline;
        *variables = *_indiv.variables;
        typename t_Container :: iterator i_var = variables->begin();
        typename t_Container :: iterator i_end = variables->end();
        for(; i_var != i_end; ++i_var )
          *i_var = ( *i_var > 0.0 ) ? 1.0 : -1.0;
        age = _indiv.age;
        if ( is_using_phenotype )
          *phenotype = *_indiv.phenotype;
        else
          phenotype = variables;
      }
      bool operator==( const Individual<FITNESS, t_Functional> &_indiv ) const
      {
        t_unsigned size = variables->size();
        if ( size != _indiv.variables->size() )
          return false;
        if ( size == 0 )
          return true;
        if ( not std::equal( variables->begin(), variables->end(),
                             _indiv.variables->begin() ) )
          return false;
        if ( not is_using_phenotype )
          return true;

        if ( not std::equal( variables->begin(), variables->end(),
                             _indiv.phenotype->begin() ) )
          return false;
        if ( not std::equal( phenotype->begin(), phenotype->end(),
                             _indiv.variables->begin() ) )
          return false;
        return std::equal( phenotype->begin(), phenotype->end(),
                           _indiv.phenotype->begin() );
      }
        
      t_Fitness fitness() const
      {
         if (not quantity_is_valid)
           throw std::runtime_error("wtf invalid fitness");
         return repFitness;
      }
      t_Type value() const
        { return quantity-baseline; }
      void set_fitness()
      {
        repFitness = quantity - baseline;
        quantity_is_valid = true;
      }

      // eo stuff
      virtual std::string className() const
        { return "Individual<...>"; }
      void printOn(std::ostream &_os) const
      {
        { // EO stuff
          if (invalid()) {
              _os << "INVALID ";
          }
          else
          {
              _os << repFitness << ' ';
          }
        } // EO stuff
        // prints ref
        _os << MotU_ref << " " << variables->size() << " ";
        const_iterator i_var  = variables->begin();
        const_iterator i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _os << *i_var << " ";
      }
      void readFrom(std::istream &_is)
      {

        { // EO stuff
          std::string fitness_str;
          t_int pos = _is.tellg();
          _is >> fitness_str;

          if (fitness_str == "INVALID")
          {
             quantity_is_valid = false;
          }
          else
          {
             quantity_is_valid = true;
             _is.seekg(pos); // rewind
             _is >> repFitness;
          }
        } // EO stuff
        t_unsigned var_size;
        _is >> var_size;
        variables->resize(var_size);

        iterator i_var = variables->begin();
        iterator i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _is >> *i_var;
      }

      void print_out( std::ostream &_stream ) const
      {
        iterator i_var = variables->begin();
        iterator i_last = variables->end();
        for( ; i_var != i_last; ++i_var )
          _stream << t_int(*i_var) << " ";
      }

      void set_genotype_to_phenotype()
        { std::copy( phenotype->begin(), phenotype->end(), variables->begin() ); }
      void set_phenotype_to_genotype()
        { std::copy( variables->begin(), variables->end(), phenotype->begin() ); }
 
      void resize( t_unsigned n )
      {
        t_Functional::resize( n );
        if ( is_using_phenotype )
          phenotype->resize( n );
      }

  };


  // redefines operators < and > with respect to FITNESS
  template <class FUNCTIONAL, class FITNESS>
  bool operator<(const Individual<FUNCTIONAL, FITNESS>& _eo1,
                 const Individual<FUNCTIONAL, FITNESS>& _eo2)
  {
      return _eo1.operator<(_eo2);
  }
  template <class FUNCTIONAL, class FITNESS>
  bool operator>(const Individual<FUNCTIONAL, FITNESS>& _eo1,
                 const Individual<FUNCTIONAL, FITNESS>& _eo2)
  {
      return _eo2.operator<(_eo1);
  }

  template<class FITNESS, class FUNCTIONAL> 
  bool Individual<FITNESS, FUNCTIONAL> :: baseline_is_valid = true; 

  template<class FITNESS, class FUNCTIONAL> 
  bool Individual<FITNESS, FUNCTIONAL> :: is_using_phenotype = false; 

  
  template<class FITNESS, class FUNCTIONAL> 
  void Individual<FITNESS, FUNCTIONAL> :: invalidate_baseline()
    { baseline_is_valid = false; }
  template<class FITNESS, class FUNCTIONAL> 
  void Individual<FITNESS, FUNCTIONAL> :: validate_baseline()
    { baseline_is_valid = true; } 
  template<class FITNESS, class FUNCTIONAL> 
  bool Individual<FITNESS, FUNCTIONAL> :: is_baseline_valid()
    { return baseline_is_valid; }

} // endif LaDa

#endif

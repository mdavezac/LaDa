//
//  Version: $Id$
//
#ifndef _DARWIN_TABOO_H_
#define _DARWIN_TABOO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ext/functional>
#include <list>
#include <algorithm>

#include <eo/eoFunctor.h>
#include <eo/eoOp.h>
#include <eo/eoGenOp.h>

#include "opt/types.h"

#include "functors.h"
#include "gatraits.h"

namespace darwin
{

  // taboo base class, declares virtual stuff
  template<class T_INDIVIDUAL>
  class Taboo_Base : public const_eoUF<const T_INDIVIDUAL&, bool>
  {
    protected:
      typedef T_INDIVIDUAL t_Individual;

    public:
      Taboo_Base() {}
      Taboo_Base( const Taboo_Base<t_Individual> &_taboo ) {}
      virtual ~Taboo_Base() {};

      types::t_unsigned max_production(void) { return 1; }

      virtual bool is_problematic() {return false;}
      virtual void set_problematic(bool _p = false) {return; }

      virtual void print_out( std::ostream &str ) const {return;}
  };

  template< class T_INDIVIDUAL, class T_CONTAINER = std::list<T_INDIVIDUAL> >
  class Taboo : public Taboo_Base<T_INDIVIDUAL>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_CONTAINER t_Container;
      typedef t_Individual value_type;

    protected:
      bool problematic;
      bool owns_pointer;
      t_Container *taboo_list;

    public:
      Taboo   ( t_Container *_list )
            : problematic(false), owns_pointer( false ),
              taboo_list(_list) {};
      Taboo() : problematic(false), owns_pointer(true)
        { taboo_list = new t_Container; }
      Taboo   ( const Taboo<t_Individual, t_Container> & _taboo )
            : owns_pointer( false ),
              taboo_list(_taboo.taboo_list) {};
      virtual ~Taboo()
      {
        if (owns_pointer and taboo_list)
          delete taboo_list;
        taboo_list = NULL;
      }

      // returns true if _indiv is in taboo_list
      virtual bool operator()( const t_Individual& _indiv ) const
      {
        typename t_Container :: const_iterator i_end = taboo_list->end();
        typename t_Container :: const_iterator i_begin = taboo_list->begin();

        return i_end != std::find( i_begin, i_end, _indiv);
      }
      

      void add( const t_Individual &_indiv, bool add_fast = true ) 
      {
        if ( not owns_pointer )
          return;
        if ( add_fast )
        {
          taboo_list->push_back( _indiv );
          return;
        }
        if (  not operator()(_indiv) )
        {
          taboo_list->push_back( _indiv );
          problematic = true;
        }
      }
      void push_back( const t_Individual &_indiv )
        { add( _indiv, false ); }

      virtual void print_out( std::ostream &str ) const
      {
        typename t_Container :: const_iterator i_pop = taboo_list->begin();
        typename t_Container :: const_iterator i_end = taboo_list->end();

        str << "Taboo Population" << std::endl;
        for(types::t_unsigned i=0 ; i_pop != i_end; ++i, ++i_pop )
        {
          str << "  Indiv " << i << " -- ";
          i_pop->print_out(str);
          str << std::endl;
        }
      };

      virtual bool is_problematic() const
        { return problematic; }
      virtual void set_problematic( bool _p = false ) 
        { problematic = _p; }

      template<class tt_Container>
      void append( const tt_Container &_pop )
      {
        if ( not owns_pointer )
          return;
        types::t_unsigned size = _pop.size();
        if ( size < 1 )  // nothing to append
          return; 
        typename tt_Container :: const_iterator i_indiv = _pop.begin();
        typename tt_Container :: const_iterator i_end = _pop.end();
        for(; i_indiv != i_end; ++i_indiv)
          if ( not operator()(*i_indiv) )
            taboo_list->push_back(*i_indiv);
      }
      typename t_Container :: const_iterator begin() const
        { return taboo_list->begin(); } 
      typename t_Container :: iterator begin() 
        { return taboo_list->begin(); } 
      typename t_Container :: const_iterator end() const
        { return taboo_list->end(); } 
      typename t_Container :: iterator end() 
        { return taboo_list->end(); } 
      types::t_unsigned size() const { return taboo_list->size(); }
      void clear() { taboo_list->clear(); }
  };

  template<class T_INDIVIDUAL, class T_INDIVTRAITS = Traits::Indiv<T_INDIVIDUAL> >
  class OffspringTaboo : public Taboo<T_INDIVIDUAL, typename T_INDIVTRAITS :: t_Population >
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIVTRAITS t_IndivTraits;
      typedef t_Individual value_type;
    protected:
      typedef typename t_IndivTraits :: t_Population t_Population;
      using Taboo<t_Individual, t_Population> :: taboo_list;

    public:
      OffspringTaboo ( t_Population *_list ) : Taboo<t_Individual, t_Population>( _list ) {}
      OffspringTaboo () : Taboo<t_Individual, t_Population>() {}
      virtual ~OffspringTaboo() {};
       
      // returns true if _indiv is in taboo_list 
      virtual bool operator()( const t_Individual& _indiv ) const
      {
        typename t_Population :: const_iterator i_end = taboo_list->end();
        typename t_Population :: const_iterator i_begin = taboo_list->begin();
        if ( i_begin == i_end )
          return false;
        --i_end; // last is current
        return i_end != std::find( i_begin, i_end, _indiv);
      }
  };

  template<class T_INDIVIDUAL, class T_CONTAINER = std::list<T_INDIVIDUAL> >
  class History : public Taboo<T_INDIVIDUAL, T_CONTAINER>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_CONTAINER t_Container;
      typedef t_Individual value_type;
    private:
      typedef typename t_Individual::t_Object t_Object;
    protected:
      using Taboo<t_Individual, t_Container> :: taboo_list;
      using Taboo<t_Individual, t_Container> :: owns_pointer;
      using Taboo<t_Individual, t_Container> :: problematic;
#ifdef _MPI
      t_Container new_taboos;
#endif 

    public:
      History() : Taboo<t_Individual, t_Container>() {}
      virtual ~History() {};

      virtual bool clone(t_Individual &_indiv)
      {
        typename t_Container :: const_iterator i_end = taboo_list->end();
        typename t_Container :: const_iterator i_indiv = taboo_list->begin();
        if ( i_indiv == i_end )
          return false;
        // t_Object since we do not want to compare fitness,
        // quantity, validity, etc...
        // but only wether these are truly different individual 
        // in terms of t_Object
        i_indiv = std::find( i_indiv, i_end, _indiv);
        if ( i_end == i_indiv )
          return false;
        _indiv.quantities() = i_indiv->quantities();
        _indiv.fitness() = i_indiv->fitness();
        return true;
      }
#ifdef _MPI
      void add( const t_Individual &_indiv, bool add_fast = true ) 
      {
        if ( not owns_pointer )
          return;
        if ( add_fast )
        {
          taboo_list->push_back( _indiv );
          new_taboos.push_back( _indiv );
          return;
        }
        if (  not operator()(_indiv) )
        {
          taboo_list->push_back( _indiv );
          new_taboos.push_back( _indiv );
          problematic = true;
        }
      }
      void synchronize()
      {
        mpi::AllGather allgather(mpi::main);
        typename t_Container :: iterator i_indiv = new_taboos.begin();
        typename t_Container :: iterator i_end = new_taboos.end();
        for(; i_indiv != i_end; ++i_indiv)
          i_indiv->broadcast(allgather);
        allgather.allocate_buffers();
        i_indiv = new_taboos.begin();
        for(; i_indiv != i_end; ++i_indiv)
          i_indiv->broadcast(allgather);
        allgather();
        new_taboos.clear();
        t_Individual indiv;
        while( indiv.broadcast(allgather) )
          { Taboo<t_Individual, t_Container>::add( indiv ); }
        new_taboos.clear();
      }

      bool broadcast( mpi::BroadCast &_bc )
      {
        types::t_int n = taboo_list->size();
        if( not _bc.serialize(n) ) return false;
        if( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
          taboo_list->resize(n);
        typename t_Container :: iterator i_indiv = taboo_list->begin();
        typename t_Container :: iterator i_indiv_end = taboo_list->end();
        for(; i_indiv != i_indiv_end; ++i_indiv )
          if ( not i_indiv->broadcast(_bc) ) return false;
        return true;
      }
#endif
  };

  // a class with a number of taboos
  template<class T_INDIVIDUAL>
  class Taboos : public Taboo_Base<T_INDIVIDUAL>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;

    protected:
      typedef Taboo_Base<t_Individual> t_Type;
      typedef std::list< t_Type* > t_Container;

    protected: 
      t_Container taboos;

    public:
      Taboos() {};
      Taboos( const Taboos<t_Individual> &_taboo ) : taboos( _taboo.taboos ) {};
      virtual ~Taboos(){};

      types::t_unsigned size() const { return taboos.size(); }
      t_Type* front() { return taboos.front(); }

      void add( Taboo_Base<t_Individual> * _taboo )
      {
        if ( _taboo == NULL )
          return;
        taboos.push_back( _taboo );
      }
      void clear()
        { taboos.clear(); } 

      // as soon as one taboo operator returns true,
      // function exits with true as well
      virtual bool operator()( const t_Individual &_indiv ) const
      {
        if ( not taboos.empty() )
        {
          typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
          typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
          for ( ; i_taboo != i_end; ++i_taboo )
            if ( (*i_taboo)->operator()( _indiv ) )
              return true;
        }

        return false;
      } 

      bool is_problematic() const
      {
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          if ( (*i_taboo)->is_problematic() )
            return true;
        return false;
      }
      void set_problematic( bool _p = false ) 
      {
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          (*i_taboo)->set_problematic( _p );
      }
      virtual void print_out( std::ostream &str ) const
      {
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          (*i_taboo)->print_out( str );
      };
  };
  
  // a class which taboos a whole list of pops
  template<class T_INDIVIDUAL, class T_INDIVTRAITS = Traits::Indiv<T_INDIVIDUAL> > 
  class IslandsTaboos : public Taboo_Base<T_INDIVIDUAL>
  {
    public:
      typedef T_INDIVIDUAL t_Individual;
      typedef T_INDIVTRAITS t_IndivTraits;
    private:
      typedef IslandsTaboos<t_Individual, t_IndivTraits>  t_Base;
      typedef typename t_IndivTraits :: t_Population  t_Container;
      typedef typename t_IndivTraits :: t_Islands     t_Islands;
      typedef typename t_IndivTraits :: t_Object t_Object;

    protected: 
      bool problematic;
      t_Islands &populations;

    public:
      IslandsTaboos   ( t_Islands &_islands )
                    : problematic(false), 
                      populations( _islands ) {};
      IslandsTaboos   ( const t_Base &_taboos )
                    : problematic(_taboos.is_problematic), 
                      populations( _taboos.populations ) {};
      virtual ~IslandsTaboos(){};

      // as soon as one taboo operator returns true,
      // function exits with true as well
      virtual bool operator()( const t_Individual &_indiv ) const
      {
        if ( populations.empty() )
          return false;

        typename t_Islands :: const_iterator i_pop = populations.begin();
        typename t_Islands :: const_iterator i_pop_end = populations.end();
        typename t_Container :: const_iterator i_end;
        for ( ; i_pop != i_pop_end; ++i_pop )
        {
          i_end = i_pop->end();
          if ( i_end != std::find( i_pop->begin(), i_end, _indiv) )
            return true;
        }

        return false;
      } 
      virtual void print_out( std::ostream &str ) const {}
      virtual bool is_problematic() const
        { return problematic; }
      virtual void set_problematic( bool _p = false ) 
        { problematic = _p; }
  };

  template<class T_INDIVIDUAL>
  class TabooOp : public eoGenOp<T_INDIVIDUAL>
  {
    protected:
      typedef T_INDIVIDUAL t_Individual; 

    protected:
      Taboo_Base<t_Individual> &taboo;
      eoMonOp<t_Individual> &utterrandom;
      types::t_unsigned max;
      eoGenOp<t_Individual> &op;

    public:
      TabooOp   ( eoGenOp<t_Individual> &_op, 
                  Taboo_Base<t_Individual> &_taboo,
                  types::t_unsigned _max,
                  eoMonOp<t_Individual> &_ur )
              : taboo(_taboo), utterrandom(_ur), max(_max), op(_op) {}

      virtual types::t_unsigned max_production()
        { return op.max_production(); }

      // tries to create an untaboo object on applying _op
      // after max tries, creates a random untaboo object
      virtual void apply( eoPopulator<t_Individual> &_indiv ) 
      {
        types::t_unsigned  i = 0;
        do
        {
          ++i;
          op( _indiv );
          if ( not taboo( *_indiv ) )
            return;
          *_indiv = _indiv.select();
        } while ( i < max );

        std::cerr << "Could not find original individual in this crowd" << std::endl
                  << "Trying random initialization" << std::endl;
        taboo.set_problematic();

        (*_indiv).invalidate(); // utterrandom is NOT a GenOp, does not invalidate!!
        do 
        {
          ++i;
          utterrandom( *_indiv ); // _indiv is a populator
          if ( not taboo( *_indiv ) )
            return;
        } while ( i < UINT_MAX );

        std::cerr << "Could not find original individual in this crowd" << std::endl
                  << "Not even from re-initialization generation" << std::endl;
        throw "";
      }

      virtual std::string className () const { return "Darwin :: TabooOp"; }

  };

  template< class T_EVALUATOR >
  class TabooFunction : public Taboo_Base< typename T_EVALUATOR::t_Individual >
  {
    public:
      typedef T_EVALUATOR t_Evaluator;
      typedef bool ( t_Evaluator::*t_Function )(const typename t_Evaluator::t_Individual &);
    protected:
      typedef typename t_Evaluator::t_Individual t_Individual;

    protected:
      t_Evaluator &evaluator;
      t_Function member_func;
      std::string class_name;

    public:
      explicit
        TabooFunction   ( t_Evaluator &_eval, t_Function _func, const std::string &_cn )
                      : evaluator(_eval), member_func(_func), class_name(_cn) {};
      ~TabooFunction() {}

      std::string className() const { return class_name; }

      bool operator()( const t_Individual& _indiv ) const
        { return ( (evaluator.*member_func) )( _indiv); }
  };

} // namespace LaDa
#endif

#ifndef _TABOO_H_
#define _TABOO_H_


#include <eo/eoOp.h>
#include <eo/eoGenOp.h>
#include <list>
#include <algorithm>

#include <opt/types.h>
using namespace types;
#include <eo/eotypes.h>

namespace LaDa
{
  // taboo base class, declares virtual stuff
  template<class t_Object>
  class Taboo_Base : public eoUF<const t_Object&, bool>
  {
    public:
      Taboo_Base() {}
      Taboo_Base( const Taboo_Base<t_Object> &_taboo ) {}
      virtual ~Taboo_Base() {};

      eotypes::t_unsigned max_production(void) { return 1; }

      virtual bool is_problematic() const = 0;
      virtual void set_problematic(bool _p = false) = 0;

      virtual void print_out( std::ostream &str ) const = 0;
  };

  template<class t_Object, class t_Container = eoPop<t_Object> >
  class Taboo : public Taboo_Base<t_Object>
  {
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
      Taboo   ( const Taboo<t_Object, t_Container> & _taboo )
            : owns_pointer( false ),
              taboo_list(_taboo.taboo_list) {};
      virtual ~Taboo()
      {
        if (owns_pointer and taboo_list)
          delete taboo_list;
        taboo_list = NULL;
      }

      // returns true if _object is in taboo_list
      virtual bool operator()( const t_Object& _object ) 
      {
        typename t_Container :: iterator i_end = taboo_list->end();
        return not ( i_end == std::find( taboo_list->begin(), i_end, _object ) );
      }
      

      void add( const t_Object &_object, bool add_fast = true ) 
      {
        if ( not owns_pointer )
          return;
        if ( add_fast )
        {
          taboo_list->push_back( _object );
          return;
        }
        if (  not operator()(_object) )
        {
          taboo_list->push_back( _object );
          problematic = true;
        }
      }

      virtual void print_out( std::ostream &str ) const
      {
        typename t_Container :: const_iterator i_pop = taboo_list->begin();
        typename t_Container :: const_iterator i_end = taboo_list->end();

        str << "Taboo Population" << std::endl;
        for(t_unsigned i=0 ; i_pop != i_end; ++i, ++i_pop )
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
        t_unsigned size = _pop.size();
        if ( size < 1 )  // nothing to append
          return; 
        typename tt_Container :: const_iterator i_indiv = _pop.begin();
        typename tt_Container :: const_iterator i_end = _pop.end();
        for(; i_indiv != i_end; ++i_indiv)
          if ( not operator()(*i_indiv) )
            taboo_list->push_back(*i_indiv);
      }
  };

  template<class t_Object, class t_Container = eoPop<t_Object> >
  class OffspringTaboo : public Taboo<t_Object, t_Container>
  {
    protected:
      using Taboo<t_Object, t_Container> :: taboo_list;

    public:
      OffspringTaboo ( t_Container *_list ) : Taboo<t_Object, t_Container>( _list ) {}
      OffspringTaboo () : Taboo<t_Object, t_Container>() {}
      virtual ~OffspringTaboo() {};
       
      // returns true if _object is in taboo_list 
      virtual bool operator()( const t_Object& _object ) 
      {
        typename t_Container :: iterator i_end = taboo_list->end();
        typename t_Container :: iterator i_begin = taboo_list->begin();
        if ( i_begin == i_end )
          return false;
        --i_end; // last is current
        return not ( i_end == std::find( i_begin, i_end, _object ) );
      }
  };

  template<class t_Object, class t_Container = eoPop<t_Object> >
  class History : public Taboo<t_Object, t_Container>
  {
    protected:
      using Taboo<t_Object, t_Container> :: taboo_list;

    public:
      History() : Taboo<t_Object, t_Container>() {}
      virtual ~History() {};

      // sets quantity in object if object is in history list
      virtual bool set_quantity(t_Object &_object)
      {
        typename t_Container :: iterator i_end = taboo_list->end();
        typename t_Container :: iterator i_indiv = taboo_list->begin();
        if ( i_indiv == i_end )
          return false;
        i_indiv = std::find( i_indiv, i_end, _object );
        if ( i_end == i_indiv )
          return false;
        _object.set_quantity( i_indiv->get_quantity() );
        _object.set_fitness();
        return true;
      }
  };

  // a class with a number of taboos
  template<class t_Object>
  class Taboos : public Taboo_Base<t_Object>
  {
    protected: 
      std::list< Taboo_Base<t_Object>* > taboos;

    public:
      Taboos() {};
      Taboos( const Taboos<t_Object> &_taboo ) : taboos( _taboo.taboos ) {};
      virtual ~Taboos(){};

      void add( Taboo_Base<t_Object> * _taboo )
      {
        if ( _taboo == NULL )
          return;
        taboos.push_back( _taboo );
      }
      void clear()
        { taboos.clear(); } 

      // as soon as one taboo operator returns true,
      // function exits with true as well
      virtual bool operator()( const t_Object &_object ) 
      {
        if ( not taboos.empty() )
        {
          typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_taboo = taboos.begin();
          typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_end = taboos.end();
          for ( ; i_taboo != i_end; ++i_taboo )
            if ( (*i_taboo)->operator()( _object ) )
              return true;
        }

        return false;
      } 

      bool is_problematic() const
      {
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          if ( (*i_taboo)->is_problematic() )
            return true;
        return false;
      }
      void set_problematic( bool _p = false ) 
      {
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          (*i_taboo)->set_problematic( _p );
      }
      virtual void print_out( std::ostream &str ) const
      {
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Object> * > :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          (*i_taboo)->print_out( str );
      };
  };
  
  // a class which taboos a whole list of pops
  template<class t_Object, class t_Container = eoPop<t_Object>, class t_Islands = std::list< t_Container > > 
  class IslandsTaboos : public Taboo_Base<t_Object>
  {
    protected: 
      bool problematic;
      t_Islands &populations;

    public:
      IslandsTaboos   ( t_Islands &_islands )
                    : problematic(false), 
                      populations( _islands ) {};
      IslandsTaboos   ( const IslandsTaboos<t_Object, t_Container, t_Islands> &_taboos )
                    : problematic(_taboos.is_problematic), 
                      populations( _taboos.populations ) {};
      virtual ~IslandsTaboos(){};

      // as soon as one taboo operator returns true,
      // function exits with true as well
      virtual bool operator()( const t_Object &_object ) 
      {
        if ( populations.empty() )
          return false;

        typename t_Islands :: iterator i_pop = populations.begin();
        typename t_Islands :: iterator i_pop_end = populations.end();
        typename t_Container :: iterator i_end;
        for ( ; i_pop != i_pop_end; ++i_pop )
        {
          i_end = i_pop->end();
          if  ( not ( i_end == std::find( i_pop->begin(), i_end, _object ) ) )
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

  template<class t_Object>
  class TabooOp : public eoGenOp<t_Object>
  {
    protected:
      Taboo_Base<t_Object> &taboo;
      UtterRandom<t_Object> utterrandom;
      t_unsigned max;
      eoGenOp<t_Object> *op;

    public:
      TabooOp   ( eoOp<t_Object> &_op, 
                  Taboo_Base<t_Object> &_taboo,
                  t_unsigned _max,
                  eoFunctorStore &_store )
              : taboo(_taboo), utterrandom(), max(_max)
        { op = &wrap_op<t_Object>( _op, _store ); }

      virtual eotypes::t_unsigned max_production()
        { return op->max_production(); }

      // tries to create an untaboo object on applying _op
      // after max tries, creates a random untaboo object
      virtual void apply( eoPopulator<t_Object> &_object ) 
      {
        t_unsigned  i = 0;
        do
        {
          eotypes::t_unsigned pos = _object.tellp();
          ++i;
          (*op)( _object );
          _object.seekp(pos);
          if ( not taboo( *_object ) )
            return;
          *_object = _object.select();
        } while ( i < max );

        std::cerr << "Could not find original individual in this crowd" << std::endl
                  << "Going UtterRandom" << std::endl;
        taboo.set_problematic();

        do 
        {
          ++i;
          utterrandom( *_object ); // _object is a populator
          if ( not taboo( *_object ) )
            return;
        } while ( i < UINT_MAX );

        std::cerr << "Could not find original individual in this crowd" << std::endl
                  << "Not even from UtterRandom generation" << std::endl;
        throw "";
      }

      virtual std::string className () const { return "LaDa :: TabooOp"; }

  };

} // namespace LaDa
#endif

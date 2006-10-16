#ifndef _TABOO_H_
#define _TABOO_H_


#include <eo/eoOp.h>
#include <eo/eoGenOp.h>
#include <list>
#include <algorithm>

namespace LaDa
{
  // taboo base class, declares virtual stuff
  template<class t_Object>
  class Taboo_Base : public eoUF<t_Object&, bool>
  {
    public:
      Taboo_Base() {}
      Taboo_Base( const Taboo_Base<t_Object> &_taboo ) {}
      virtual ~Taboo_Base() {};

      unsigned max_production(void) { return 1; }

      virtual bool is_problematic() const = 0;
      virtual void set_problematic(bool _p = false) = 0;

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
      virtual bool operator()( t_Object& _object ) 
      {
        typename t_Container :: iterator i_end = taboo_list->end();
        return not ( i_end == std::find( taboo_list->begin(), i_end, _object ) );
      }
      

      void add( const t_Object &_object ) 
      {
        if ( owns_pointer )
        {
          taboo_list->push_back( _object );
          problematic = true;
        }
      }

      void print_out( std::ostream &str ) const
      {
        typename t_Container :: const_iterator i_pop = taboo_list->begin();
        typename t_Container :: const_iterator i_end = taboo_list->end();

        str << "Taboo Population" << std::endl;
        for(unsigned i=0 ; i_pop != i_end; ++i, ++i_pop )
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
      virtual bool operator()( t_Object &_object ) 
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
  };

  template<class t_Object>
  class TabooOp : public eoGenOp<t_Object>
  {
    protected:
      Taboo_Base<t_Object> &taboo;
      UtterRandom<t_Object> utterrandom;
      unsigned max;
      eoGenOp<t_Object> *op;

    public:
      TabooOp   ( eoOp<t_Object> &_op, 
                  Taboo_Base<t_Object> &_taboo,
                  unsigned _max,
                  eoFunctorStore &_store )
              : taboo(_taboo), utterrandom(), max(_max)
        { op = &wrap_op<t_Object>( _op, _store ); }

      virtual unsigned max_production()
        { return op->max_production(); }

      // tries to create an untaboo object on applying _op
      // after max tries, creates a random untaboo object
      virtual void apply( eoPopulator<t_Object> &_object ) 
      {
        unsigned  i = 0;
        do
        {
          unsigned pos = _object.tellp();
          ++i;
          (*op)( _object );
          _object.seekp(pos);
          if ( not taboo( *_object ) )
            return;
        } while ( i < max );

        do 
        {
          ++i;
          utterrandom( *_object ); // _object is a populator
          if ( not taboo( *_object ) )
            return;
        } while ( i < UINT_MAX );

        std::cerr << "Could not find original individual in this crowd" << std::endl;
        throw "";
      }

      virtual std::string className () const { return "LaDa :: TabooOp"; }

  };

} // namespace LaDa
#endif

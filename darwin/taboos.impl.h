//
//  Version: $Id$
//

#ifndef _DARWIN_TABOO_IMPL_H_
#define _DARWIN_TABOO_IMPL_H_

#include "debug.h"

namespace GA 
{
  template< class T_INDIVIDUAL, class T_CONTAINER >
    Taboo<T_INDIVIDUAL, T_CONTAINER> :: ~Taboo()
    {
      if (owns_pointer and taboo_list)
        delete taboo_list;
      taboo_list = NULL;
    }
  template< class T_INDIVIDUAL, class T_CONTAINER >
    inline void Taboo<T_INDIVIDUAL, T_CONTAINER> :: 
      add( const t_Individual &_indiv, bool add_fast) 
      {
        if ( not owns_pointer ) return;
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
  
  template< class T_INDIVIDUAL, class T_CONTAINER >
    inline bool Taboo<T_INDIVIDUAL, T_CONTAINER> ::  
      operator()( const t_Individual& _indiv ) const
      {
        typename t_Container :: const_iterator i_end = taboo_list->end();
        typename t_Container :: const_iterator i_begin = taboo_list->begin();

        return i_end != std::find( i_begin, i_end, _indiv);
      }

  template< class T_INDIVIDUAL, class T_CONTAINER >
    void Taboo<T_INDIVIDUAL, T_CONTAINER> :: print_out( std::ostream &str ) const
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

  template< class T_INDIVIDUAL, class T_CONTAINER > template<class tt_Container>
    void Taboo<T_INDIVIDUAL, T_CONTAINER> ::
      append( const tt_Container &_pop )
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




  template<class T_GATRAITS>
    inline bool OffspringTaboo<T_GATRAITS> :: 
      operator()( const t_Individual& _indiv ) const
      {
        typename t_Population :: const_iterator i_end = taboo_list->end();
        typename t_Population :: const_iterator i_begin = taboo_list->begin();
        if ( i_begin == i_end )
          return false;
        --i_end; // last is current
        return i_end != std::find( i_begin, i_end, _indiv);
      }







  template<class T_INDIVIDUAL, class T_CONTAINER >
    bool History<T_INDIVIDUAL, T_CONTAINER> :: clone(t_Individual &_indiv)
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
      _indiv.quantities() = i_indiv->const_quantities();
      _indiv.fitness() = i_indiv->fitness();
      return true;
    }

#ifdef _MPI
  template<class T_INDIVIDUAL, class T_CONTAINER >
    void History<T_INDIVIDUAL, T_CONTAINER> ::
      add( const t_Individual &_indiv, bool add_fast ) 
      {
        if ( not owns_pointer ) return;
        if ( add_fast )
        {
          taboo_list->push_back( _indiv );
          new_taboos.push_back( _indiv );
          return;
        }
        if ( not operator()(_indiv) )
        {
          taboo_list->push_back( _indiv );
          new_taboos.push_back( _indiv );
          problematic = true;
        }
      }
  template<class T_INDIVIDUAL, class T_CONTAINER >
    void History<T_INDIVIDUAL, T_CONTAINER> :: synchronize()
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

  template<class T_INDIVIDUAL, class T_CONTAINER >
    bool History<T_INDIVIDUAL, T_CONTAINER> :: broadcast( mpi::BroadCast &_bc )
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





  template<class T_INDIVIDUAL>
    inline void Taboos<T_INDIVIDUAL> :: add( Taboo_Base<t_Individual> * _taboo )
    {
      if ( _taboo == NULL )
        return;
      taboos.push_back( _taboo );
    }


  template<class T_INDIVIDUAL>
    bool Taboos<T_INDIVIDUAL> :: operator()( const t_Individual &_indiv ) const
    {
      if ( not taboos.empty() )
      {
        typename std::list< Taboo_Base<t_Individual> * >
            :: const_iterator i_taboo = taboos.begin();
        typename std::list< Taboo_Base<t_Individual> * >
            :: const_iterator i_end = taboos.end();
        for ( ; i_taboo != i_end; ++i_taboo )
          if ( (*i_taboo)->operator()( _indiv ) )
            return true;
      }

      return false;
    } 

  template<class T_INDIVIDUAL>
    bool Taboos<T_INDIVIDUAL> :: is_problematic() const
    {
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
      for ( ; i_taboo != i_end; ++i_taboo )
        if ( (*i_taboo)->is_problematic() )
          return true;
      return false;
    }
  template<class T_INDIVIDUAL>
    void Taboos<T_INDIVIDUAL> :: set_problematic( bool _p )
    {
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
      for ( ; i_taboo != i_end; ++i_taboo )
        (*i_taboo)->set_problematic( _p );
    }

  template<class T_INDIVIDUAL>
    void Taboos<T_INDIVIDUAL> :: print_out( std::ostream &str ) const
    {
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_taboo = taboos.begin();
      typename std::list< Taboo_Base<t_Individual> * > :: const_iterator i_end = taboos.end();
      for ( ; i_taboo != i_end; ++i_taboo )
        (*i_taboo)->print_out( str );
    };





  template<class T_GATRAITS>
    bool IslandsTaboos<T_GATRAITS> :: operator()( const t_Individual &_indiv ) const
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







  template<class T_INDIVIDUAL>
    void TabooOp<T_INDIVIDUAL> :: apply( eoPopulator<t_Individual> &_indiv ) 
    {
      types::t_unsigned  i = 0;
      do
      {
        ++i;
        op( _indiv );
        if ( not taboo( *_indiv ) ) return;
        *_indiv = _indiv.select();
      } while ( i < max );

      taboo.set_problematic();

      (*_indiv).invalidate(); // utterrandom is NOT a GenOp, does not invalidate!!
      do 
      {
        ++i;
        utterrandom( *_indiv ); // _indiv is a populator
        if ( not taboo( *_indiv ) )
          return;
      } while ( i < UINT_MAX );

      std::ostringstream sstr;
      sstr << __LINE__ << ", line: " << __LINE__ << "\n"
           << "Could not create a non-taboo individual\n";
      throw std::runtime_error( sstr.str() );
    }



} // namespace GA

#endif

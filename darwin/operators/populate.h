//
//  Version: $Id$
//

#ifndef _LADA_GA_POPULATE_
#define _LADA_GA_POPULATE_

#include <boost/lambda/bind.hpp>

namespace LaDa
{
  namespace GA
  {
    namespace Operators
    {
      //! \brief Grows a population to a set size \a _popsize using operator \a
      //!        _operator to create individuals.

      template< class T_OPERATOR, class T_TABOOS >
        struct Populate
        {
          typedef void result_type;
          T_OPERATOR const& operator_;
          T_TABOOS const& taboos_;
          Populate   ( T_OPERATOR const& _op, T_TABOOS const &_tab )
                   : operator_(_op), taboos_(_tab) {}
          Populate   ( Populate const& _c )
                   : operator_(_c.operator_), taboos_(_c.taboos_) {}
                  

          template< class T_POPULATION >
            result_type operator()( T_POPULATION& _pop,  size_t _popsize )
            {
              if( _pop.size() == _popsize ) return;
              typedef typename T_POPULATION :: value_type t_Individual;
              types::t_unsigned i = 0, j = 0;
              types::t_unsigned _maxiter = 50 * _popsize;
              while ( _pop.size() < _popsize and i < _maxiter)
              {
                t_Individual indiv;
                operator_( indiv );
                if ( not taboos_( indiv ) )
                {
                  _pop.push_back(indiv);
                    ++j;
                }
                ++i;
              } // while ( i_pop->size() < target )
              __DOASSERT( _pop.size() < _popsize,
                             "Created " << j << " individuals in " << i << " iterations.\n"
                          << "Are taboos/concentration constraints to restrictive?\n" )
              Print::out << "Created " << j << " new individuals.\n";
              _pop.resize( _popsize );
            }
        };

      template< class T_OPERATOR, class T_TABOOS >
        Populate<T_OPERATOR, T_TABOOS> populator( T_OPERATOR const& _op, T_TABOOS &_tab )
          { return Populate<T_OPERATOR, T_TABOOS>( _op, _tab); }

      //! Could use a lambda... but can't get it right...
      template< class T_INDIVIDUAL >
        void call_object_mask( T_INDIVIDUAL& _indiv, size_t _first, size_t _end )
          { _indiv.Object().mask( _first, _end ); }

      
      //! \brief Grows a population: (i) by creating an indivdual using \a _initop
      //!        and (ii) adding "masked" individual using \a _maskop.
      template< class T_INITOP, class T_MASKOP, class T_TABOOS, class T_POPULATION >
        void partition_populate( const T_INITOP& _initop, const T_MASKOP& _maskop,
                                 const T_TABOOS& _taboos, T_POPULATION& _pop, 
                                 size_t _popsize, size_t _maxiter )
        {
          if( _pop.size() == _popsize ) return;
          typedef typename T_POPULATION :: value_type t_Individual;
          size_t i = 0, popstart( _pop.size() );
          while ( _pop.size() < _popsize and i < _maxiter )
          {
            ++i;
            // first random object
            t_Individual indiv;
            _initop( indiv ); 
            types::t_unsigned start = 0;  
            types::t_unsigned end = indiv.Object().size();  
          
            // indiv: ----
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
            
            _maskop( indiv, start, end );         // indiv: ****
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
            
            _maskop( indiv, start, end/2 );       // indiv: --**
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: **--
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end/4 );       // indiv: -*--
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: *-**
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, end/2, end*3/4 );     // indiv: *--*
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: -**-
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, end/4, end/2 );       // indiv: --*-
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: **-*
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end/2 );       // indiv: ---*
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: *---
            if ( not _taboos(indiv) ) _pop.push_back(indiv);
          }
          __DOASSERT( _pop.size() < _popsize,
                         "Created " << _pop.size() - popstart << " individuals in "
                      << i << " iterations.\n"
                         "Are taboos/concentration constraints to restrictive?\n" )
          _pop.resize( _popsize );
        }

    } // namespace Operators
  } // namespace GA
} // namespace LaDa

#endif

#ifndef LADA_LNS_PARSED_TREE_BASE_H
#define LADA_LNS_PARSED_TREE_BASE_H

#include "LaDaConfig.h"

#include <iostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include "../string_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace tree
    {
      namespace details
      {
        template<class T_SECTION, class T_FILTER> class Base;

        //! Filters according to name string.
        template<class T_SECTION>
          class NameFilter
          {
              typedef T_SECTION t_Section;
            public:
              //! Constructor.
              NameFilter   ( t_String const & _filter = "" ) 
                         : with_filter_( _filter.size() > 0 ), filter_(_filter) {}
              //! Copy Constructor.
              NameFilter   ( NameFilter const & _c )
                         : with_filter_(_c.with_filter_), filter_(_c.filter_) {}
              //! Predicate.
              bool operator()( t_Section const &_section ) const
                { return with_filter_ ? filter_ == _section.name: true; }

            private:
              //! Whether this iterator shall be filtering or not.
              bool with_filter_;
              //! The filter.
              t_String filter_;
          };

        //! Head node of the xml tree.
        template<class T_SECTION, class T_FILTER = NameFilter<T_SECTION> >
          class Base 
          {
            public:
              //! Type of the subsections.
              typedef T_SECTION t_Section; 
              //! Type of the subsections.
              typedef T_FILTER t_Filter; 
            public:
              //! Type of option tree.
              typedef std::vector<t_Section> t_Subsections;
              //! Iterators
              struct iterator
              {
                //! Type of the subsections iterator.
                typedef boost::filter_iterator
                        <
                          t_Filter, 
                          typename t_Subsections::iterator
                        > subsection;
              };
              //! Iterators
              struct const_iterator
              {
                //! Type of the subsections iterator.
                typedef boost::filter_iterator
                        <
                          t_Filter, 
                          typename t_Subsections::const_iterator
                        > subsection;
              };
      
              //! Constructor.
              Base() {}
              //! Copy Constructor.
              Base( Base const& _c ) : subsections_( _c.subsections_ ) {}
      
              //! Subsection iterator.
              typename iterator::subsection subsections_begin();
              //! Subsection iterator.
              typename const_iterator::subsection subsections_begin() const;
              //! Subsection iterator over a subset of subsections.
              typename iterator::subsection subsections_begin( t_String const& _op );
              //! Subsection iterator.
              typename const_iterator::subsection subsections_begin( t_String const& _op ) const;
              //! Subsection iterator.
              typename iterator::subsection subsections_end();
              //! Subsection iterator.
              typename const_iterator::subsection subsections_end() const;
              //! Subsection iterator over a subset of subsections.
              typename iterator::subsection subsections_end( t_String const& _op );
              //! Subsection iterator.
              typename const_iterator::subsection subsections_end( t_String const& _op ) const;
      
              //! Erase a subsection.
              void erase( typename iterator::subsection const& _it ) 
                { if( subsections_ ) subsections_->erase( _it.base() ); }
      
              //! Insert a subsection.
              void push_back( const t_Section &_value );
      
              //! Copies a section.
              void copy_to( Base& _c ) const;

              //! Swaps to objects.
              void swap( Base& _c ) { subsections_.swap(_c.subsections_); }
      
            private:
              //! A map of sections.
              boost::shared_ptr< t_Subsections > subsections_;
          };

      //! Thin wrapper to change routine names.
      template<class T_OPTION>
        class Options 
        {
            //! Holds implementation.
            details::Base<T_OPTION> impl_;
          public:
            //! Iterator.
            typedef typename Base<T_OPTION> :: iterator :: subsection iterator;
            //! Iterator.
            typedef typename Base<T_OPTION> :: const_iterator :: subsection const_iterator;
          
            //! Constructor.
            Options() : impl_() {}
            //! Copy Constructor.
            Options( Options const& _c ) : impl_( _c.impl_ ) {}
            
            //! Option iterator.
            iterator options_begin()
              { return impl_.subsections_begin(); }
            //! Option iterator.
            const_iterator options_begin() const
              { return impl_.subsections_begin(); };
            //! Option iterator over a subset of options.
            iterator options_begin( t_String const& _op )
              { return impl_.subsections_begin(_op); }
            //! Option iterator.
            const_iterator options_begin( t_String const& _op ) const
              { return impl_.subsections_begin(_op); }
            //! Option iterator.
            iterator options_end()
              { return impl_.subsections_end(); }
            //! Option iterator.
            const_iterator options_end() const
              { return impl_.subsections_end(); }
            //! Option iterator over a subset of options.
            iterator options_end( t_String const& _op )
              { return impl_.subsections_end(_op); }
            //! Option iterator.
            const_iterator options_end( t_String const& _op ) const
              { return impl_.subsections_end(_op); }
          
            //! Erase an option.
            void erase( iterator const& _it ) 
              { return impl_.erase(_it); }
            //! Insert an option.
            void push_back( T_OPTION const &_op )
              { return impl_.push_back( _op); }
            //! Insert an option.
            void copy_to( Options &_op ) const
              { return impl_.copy_to( _op.impl_ ); }
            //! Swaps to objects.
            void swap( Options& _c ) { impl_.swap(_c.impl_); }
        };

      template<class T_SECTION, class T_FILTER>
        typename Base<T_SECTION, T_FILTER>::iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_begin() 
          { 
            if( not subsections_ ) 
            {
              static typename iterator::subsection a;
              return a;
            }
            return typename iterator::subsection( subsections_->begin(), subsections_->end() );
          }

      template<class T_SECTION, class T_FILTER>
        typename Base<T_SECTION, T_FILTER>::const_iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_begin() const
          { 
            if( not subsections_ ) 
            {
              static typename const_iterator::subsection a;
              return a;
            }
            return typename const_iterator::subsection( subsections_->begin(), subsections_->end() );
          }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_begin( t_String const& _op )
          {
            return subsections_ ? 
                     typename iterator::subsection
                     (
                       t_Filter(_op), 
                       subsections_->begin(), 
                       subsections_->end()
                     )
                   : subsections_begin();
          }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::const_iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_begin( t_String const& _op ) const
          {
            return subsections_ ? 
                     typename const_iterator::subsection
                     (
                        t_Filter(_op), 
                        subsections_->begin(), 
                        subsections_->end()
                     )
                   : subsections_begin();
          }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_end( t_String const& _op )
          {
            return subsections_ ? 
                     typename iterator::subsection
                     (
                       t_Filter(_op), 
                       subsections_->end(), 
                       subsections_->end()
                     )
                   : subsections_begin();
          }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::const_iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_end( t_String const& _op ) const
          {
            return subsections_ ? 
                     typename const_iterator::subsection
                     (
                       t_Filter(_op), 
                       subsections_->end(), 
                       subsections_->end()
                     )
                   : subsections_begin();
          }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::iterator::subsection Base<T_SECTION, T_FILTER>::subsections_end()
        {
          return subsections_ ? 
                   typename iterator::subsection
                   (
                     subsections_->end(), 
                     subsections_->end()
                   )
                 : subsections_begin();
        }
      template<class T_SECTION, class T_FILTER>
        inline typename Base<T_SECTION, T_FILTER>::const_iterator::subsection
          Base<T_SECTION, T_FILTER>::subsections_end() const
          {
            return subsections_ ? 
                     typename const_iterator::subsection
                     (
                       subsections_->end(), 
                       subsections_->end()
                     )
                   : subsections_begin();
          }
        
      template<class T_SECTION, class T_FILTER>
        void Base<T_SECTION, T_FILTER>::push_back( T_SECTION const &_value )
        {
          if( not subsections_ ) subsections_.reset( new t_Subsections );
          subsections_->push_back(_value);
        }
        
      template<class T_SECTION, class T_FILTER>
        void Base<T_SECTION, T_FILTER>::copy_to( Base<T_SECTION, T_FILTER>& _c ) const
        {
          if( subsections_ ) _c.subsections_.reset( new t_Subsections( *subsections_ ) );
          else if( _c.subsections_ ) 
          {
            boost::shared_ptr<t_Subsections> s;
            _c.subsections_.swap(s);
          }
        }
      } // end of namespace details.

    } // namespace parser
  } // namespace xml
} // namespace LaDa

#endif

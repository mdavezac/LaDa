//
//  Version: $Id: subsections.h 1250 2009-07-26 21:04:07Z davezac $
//

#ifndef _LADA_LNS_XPR_SUBSECTIONS_H_
#define _LADA_LNS_XPR_SUBSECTIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <string>
#include <list>

#include <boost/shared_ptr.hpp>

#include "../tree/base.h"
#include "../string_type.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      //! Head node of the xml tree.
      template<class T_SECTION>
        class Subsections
        {
          public:
            //! Type of the subsections.
            typedef T_SECTION t_Section; 
          public:
            //! Type of option tree.
            typedef std::list<t_Section> t_Subsections;
            //! Iterators
            typedef typename t_Subsections::iterator iterator;
            //! Iterators
            typedef typename t_Subsections::const_iterator const_iterator;
            //! Constructor.
            Subsections() {}
            //! Copy Constructor.
            Subsections( Subsections const& _c ) : subsections_( _c.subsections_ ) {}
    
            //! Subsection iterator.
            iterator begin();
            //! Subsection iterator.
            const_iterator begin() const;
            //! Subsection iterator.
            iterator end()
              { return subsections_ ?  subsections_->end(): begin(); }
            //! Subsection iterator.
            const_iterator end() const
              { return subsections_ ?  subsections_->end(): begin(); }
    
            //! Erase a subsection.
            void erase( iterator const& _it ) 
              { if( subsections_ ) subsections_->erase( _it ); }
    
            //! Insert a subsection.
            void push_back( const t_Section &_value );
    
            //! Copies a section.
            void copy_to( Subsections& _c ) const;

            //! Swaps to objects.
            void swap( Subsections& _c ) { subsections_.swap(_c.subsections_); }
            
            //! Prints to string.
            t_String print(t_String const& _depth, t_String const& _tab ) const;
    
          private:
            //! A map of sections.
            boost::shared_ptr< t_Subsections > subsections_;
        };

      template<class T_SECTION>
        typename Subsections<T_SECTION>::iterator Subsections<T_SECTION>::begin() 
        { 
          if( not subsections_ ) 
          {
            static iterator a;
            return a;
          }
          return subsections_->begin();
        }
      
      template<class T_SECTION>
        typename Subsections<T_SECTION>::const_iterator Subsections<T_SECTION>::begin() const
        { 
          if( not subsections_ ) 
          {
            static const_iterator a;
            return a;
          }
          return subsections_->begin();
        }
      
        
      template<class T_SECTION>
        void Subsections<T_SECTION>::push_back( T_SECTION const &_value )
        {
          if( not subsections_ ) subsections_.reset( new t_Subsections );
          subsections_->push_back(_value);
        }
        
      template<class T_SECTION>
        void Subsections<T_SECTION>::copy_to( Subsections<T_SECTION>& _c ) const
        {
          if( subsections_ ) _c.subsections_.reset( new t_Subsections( *subsections_ ) );
          else if( _c.subsections_ ) 
          {
            boost::shared_ptr<t_Subsections> s;
            _c.subsections_.swap(s);
          }
        }
      
      template<class T_SECTION>
        t_String Subsections<T_SECTION>::print(t_String const& _depth, t_String const& _tab ) const 
        {
          if( (not subsections_) or subsections_->size() == 0 )
            return " none.\n";
          t_String result ("\n");
          for( const_iterator i_first( begin() ), i_end( end() ); 
               i_first != i_end; ++i_first )
               result += i_first->print(_depth, _tab);
          return result;
        }
    } // end of namespace details.

  } // namespace xml
} // namespace LaDa

#endif

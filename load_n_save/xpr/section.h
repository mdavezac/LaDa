#ifndef LADA_LNS_XPR_SECTION_H
#define LADA_LNS_XPR_SECTION_H

#include "LaDaConfig.h"

#include <boost/xpressive/detail/utility/tracking_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include "../string_type.h"
#include "../tags.h"
#include "../parser_base.h"
#include "../sequencer/binary.h"
#include "option.h"
#include "subsections.h"
#include "section_data.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      namespace details
      {
        namespace tracking = boost::xpressive::detail;

        //! Holds section implementation.
        template< class T_SECTION >
          struct SectionImpl : public tracking::enable_reference_tracking< SectionImpl<T_SECTION> >
          {
            //! Tracking base.
            typedef tracking::enable_reference_tracking< SectionImpl<T_SECTION> > t_TrackingBase;
            //! Subsections container type.
            typedef Subsections<T_SECTION> t_Subsections;
            //! Options container type.
            typedef tree::details::Options<Option> t_Options;
            //! Iterators.
            struct iterator
            {
              //! Subsection iterator.
              typedef typename t_Subsections::iterator subsection;
              //! Option iterator.
              typedef t_Options::iterator option;
            };
            //! Iterators.
            struct const_iterator
            {
              //! Subsection iterator.
              typedef typename t_Subsections::const_iterator subsection;
              //! Option iterator.
              typedef t_Options::const_iterator option;
            };
  
            //! Constructor.
            SectionImpl() : t_TrackingBase(), 
                            subsections_( new t_Subsections),
                            options_( new t_Options),
                            sequence_(new sequencer::Binary) {}
            //! Copy Constructor.
            SectionImpl   (const SectionImpl &_c)
                        : t_TrackingBase(_c),
                          subsections_(_c.subsections_), options_(_c.options_),
                          data_(_c.data_), sequence_(_c.sequence_)  {}

            void swap(SectionImpl &_b)
            {
              t_TrackingBase::swap(_b);
              subsections_.swap( _b.subsections_ );
              options_.swap( _b.options_ );
              sequence_.swap(_b.sequence_);
              data_.swap(_b.data_);
            }

            //! Sets the data.
            template<class T_DATA> 
              void set_data( T_DATA &_data )
               { data_.reset( new SectionData<T_SECTION, T_DATA>( _data ) ); }

            //! Prints out this section.
            t_String print(t_String const& _depth, t_String const& _tab) const;

            //! Parsing double dispatch.
            bool parse( parser_base::Section const& _sec,
                        T_SECTION const& _this, version_type _version) const
              { return data_ ? data_->parse(_sec, _this, _version): _sec.content(_this); }

            //! Acquires \a _sec as subsections.
            void push_back( Section const& _sec, SectionImpl & _impl )
            {
              track_reference( _impl );
              return subsections_->push_back(_sec); 
            }

            //! Subsections.
            boost::shared_ptr< t_Subsections > subsections_;
            //! Options.
            boost::shared_ptr< t_Options > options_;
            //! Holds data of section.
            boost::shared_ptr< SectionDataBase<T_SECTION> > data_;
            //! Tree of options and subsections.
            boost::shared_ptr<sequencer::Binary> sequence_;
            //! Can't copy directly.
            void operator=( SectionImpl<T_SECTION> const& );
          };

        template<class T_SECTION>
          t_String SectionImpl<T_SECTION>::print(t_String const& _depth, t_String const &_tab ) const
          {
            std::ostringstream sstr;
            if( not data_ ) sstr << "Base:\n";
            else sstr << data_->print( _depth, _tab );
            t_String newdepth = _depth + _tab;
            sstr << newdepth << "subsections: "
                 << this->subsections_->print( newdepth + _tab, _tab );
            return sstr.str();
          }
  
      } // namespace details.
 
      //! \cond
      class Section;
      //! \endcond 

      //! Two sections of equivalent depth.
      Section operator&&( Section _a, Section _b );
      //! Two options of equivalent depth.
      Section operator&&( Option const& _a, Option const &_b );
      //! Content and an option of equivalent depth.
      Section operator&&( Section _a, Option const &_b );
      //! Content and an option of equivalent depth.
      Section operator&&( Option const& _b, Section _a );

      //! Two options.
      Section operator||( Option const &_a, Option const &_b );
      //! Content and an option of equivalent depth
      Section operator||( Section _a, Option const &_b );
      //! Content and an option of equivalent depth
      Section operator||( Option const &_b, Section _a );
      //! Two sections of equivalent depth.
      Section operator||( Section _a, Section _b );
       
    
      //! Insert subsection.
      Section operator<<( Section _a, Section const& _b );
      //! Insert option.
      Section operator<<( Section _a, Option const& _b );

      //! Operator node.
      class Section
      {
          //! Type of the implementation.
          typedef details::SectionImpl<Section> t_Impl;
          //! Type of the tracking ptr.
          typedef details::tracking::tracking_ptr< t_Impl > t_Tracking;
          //! Type of the options implementation.
          typedef t_Impl::t_Options t_OptionsImpl;
          //! Type of the subsections implementation.
          typedef t_Impl::t_Subsections t_SubsectionsImpl;

        public:
          //! Iterators.
          struct iterator
          {
            //! Options iterator.
            typedef t_Impl::iterator::option option;
            //! Sections iterator.
            typedef t_Impl::iterator::subsection subsection;
          };
          //! Iterators.
          struct const_iterator
          {
            //! Options iterator.
            typedef t_Impl::const_iterator::option option;
            //! Sections iterator.
            typedef t_Impl::const_iterator::subsection subsection;
          };


          //! Constructor.
          Section() : impl_() {}
          //! Copy Constructor.
          Section   ( Section const& _c ) : impl_(_c.impl_) {}

          //! Option iterator.
          iterator::option options_begin()
            { return impl_->options_->options_begin(); }
          //! Option iterator.
          const_iterator::option options_begin() const
            { return impl_->options_->options_begin(); };
          //! Option iterator over a subset of options.
          iterator::option options_begin( t_String const& _op )
            { return impl_->options_->options_begin(_op); }
          //! Option iterator.
          const_iterator::option options_begin( t_String const& _op ) const
            { return impl_->options_->options_begin(_op); }
          //! Option iterator.
          iterator::option options_end()
            { return impl_->options_->options_end(); }
          //! Option iterator.
          const_iterator::option options_end() const
            { return impl_->options_->options_end(); }
          //! Option iterator over a subset of options.
          iterator::option options_end( t_String const& _op )
            { return impl_->options_->options_end(_op); }
          //! Option iterator.
          const_iterator::option options_end( t_String const& _op ) const
            { return impl_->options_->options_end(_op); }
          //! Erase an option.
          void erase( iterator::option const& _it ) 
            { return impl_.get()->options_->erase(_it); }
          //! Insert an option.
          void push_back( Option const&_op )
            { return impl_.get()->options_->push_back( _op); }

          //! Subsection iterator.
          iterator::subsection subsections_begin()
            { return impl_->subsections_->begin(); }
          //! Subsection iterator.
          const_iterator::subsection subsections_begin() const
            { return impl_->subsections_->begin(); };
          //! Subsection iterator.
          iterator::subsection subsections_end()
            { return impl_->subsections_->end(); }
          //! Subsection iterator.
          const_iterator::subsection subsections_end() const
            { return impl_->subsections_->end(); }
          //! Erase a subsection.
          void erase( iterator::subsection const& _it ) 
            { return impl_->subsections_->t_SubsectionsImpl::erase(_it); }
          //! Insert a subsection.
          void push_back( Section const& _sec ) { impl_.get()->push_back(_sec, *_sec.impl_.get() ); }

          //! Sets the data.
          template<class T_DATA> void set_data(T_DATA &_data) { impl_.get()->set_data(_data); }

          //! Prints out to a string.
          t_String print( t_String const& _depth = "", t_String const &_tab = "  " ) const
            { return impl_->print(_depth, _tab); }

          //! Whether a section is complete with name, tag... or is a content holder.
          bool incomplete() const { return not impl_->data_; }

          //! Parsing double dispatch.
          bool parse(parser_base::Section const& _sec, version_type _version) const
            { return impl_->parse(_sec, *this, _version); }

          //! Returns true if has data.
          bool has_data() const { return impl_->data_; }

          //! Returns reference to sequencer list.
          sequencer::Binary& sequence() const { return *impl_->sequence_; }

          //! Returns tag.
          tags tag() const { return impl_->data_->tag(); }

          //! Merges input's subsections and options into self.
          Section merge(Section const &_in);

        private:
          //! Handle to the body.
          t_Tracking impl_;
      };


    } // namespace xpr
  } // namespace load_n_save
} // namespace LaDa

namespace boost
{
  namespace xpressive
  {
    namespace detail
    {
      template< class T_SECTION >
        void swap( LaDa::load_n_save::xpr::details::SectionImpl<T_SECTION> &_a,
                   LaDa::load_n_save::xpr::details::SectionImpl<T_SECTION> &_b  )
          { _a.swap(_b); }
    }
  }
}

#endif

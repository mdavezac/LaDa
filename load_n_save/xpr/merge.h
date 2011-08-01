#ifndef LADA_LNS_XPR_MERGE_H
#define LADA_LNS_XPR_MERGE_H

#include "LaDaConfig.h"

#include <boost/exception/get_error_info.hpp>
#include "../exceptions.h"
#include "../access.h"
#include "utilities.h"
#include "section.h"


namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      namespace details
      {
        //! \brief Merges a section into the first section.
        //! \details Useful to add new objects from a derived class into the
        //! section declared by its base. 
        template<class T0, class T1>
          struct Merge
          {
            struct CatchSection : public parser_base::Section
            {
              CatchSection(bool _loading) : loading_(loading_), version_(0u) {};
              virtual bool operator()(xpr::Section const &_object, const unsigned int _version) const
                { version_ = _version; return _object.parse(*this, _version); }
              virtual bool operator&(xpr::Section const &_object) const 
                { return operator()(_object, version_); }

              bool regular(xpr::Section const &_object, regular_data const &) const
                { return content(_object); }
              bool content(xpr::Section const &_object, t_String const &_name = "") const
                { result =  _object; return true;  }
              bool is_loading() const { return loading_; }
              bool grammar_only() const { return true; }
              
              //! Starts recurrence;
              virtual void start_recurrence() const
                {  BOOST_THROW_EXCEPTION(error::cannot_merge_recurrent()); }
              
              //! Increments recurrence.
              virtual void step_recurrence() const 
                {  BOOST_THROW_EXCEPTION(error::cannot_merge_recurrent()); }

              //! Starts recurrence;
              virtual void stop_recurrence() const 
                {  BOOST_THROW_EXCEPTION(error::cannot_merge_recurrent()); }

              mutable bool loading_;
              mutable unsigned long version_;
              mutable xpr::Section result; 
            };
            
            Merge(T0 &_a, T1 &_b) : a(_a), b(_b) {}
            Merge(Merge const &_c) : a(_c.a_), b(_c.b_) {}

            template<class T_ARCHIVE>
              bool lns_access(T_ARCHIVE &_ar, const unsigned int _version)
              {
                CatchSection catchsection(_ar.is_loading());
                catchsection(ext(a), _version);
                Section sec0 = catchsection.result; 
                catchsection(ext(b), _version);
                Section sec1 = catchsection.result; 
                return _ar & sec0.merge(sec1); 
              }
            T0 &a;
            T1 &b;
          };
      } // namespace details
    } // namespace xpr

    //! \brief Merges a section into the first section.
    //! \details Useful to add new objects from a derived class into the
    //! section declared by its base. 
    template<class T0, class T1> xpr::Section merge(T0 &_merger, T1 &_mergee)
    {
      boost::shared_ptr< xpr::details::Merge<T0, T1> >
        ptr(new xpr::details::Merge<T0, T1>(_merger, _mergee) );
      return ext(ptr);
    }
    //! \brief Merges a section into the first section.
    //! \details Useful to add new objects from a derived class into the
    //! section declared by its base. 
    template<class T0, class T1> xpr::Section merge(T0 *_merger, T1 *_mergee)
      { return merge(*_merger, *_mergee); }
  } // namespace load_n_save
} // namespace LaDa
#endif


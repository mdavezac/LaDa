#ifndef LADA_LNS_XPR_PUSHBACK_H
#define LADA_LNS_XPR_PUSHBACK_H

#include "LaDaConfig.h"

#include <boost/exception/get_error_info.hpp>
#include "../exceptions.h"
#include "utilities.h"

namespace LaDa 
{
  namespace load_n_save
  {
    namespace xpr
    {
      namespace details
      {
        //! LNS for containers of external types, eg vector of atoms.
        template<class T>
          struct PushBack
          {
            PushBack(T &_container, int _required = 0) : container_(_container), required_(_required) {}
            PushBack(PushBack<T> const &_c) : container_(_c.container_), required_(_c.required_) {}
            template<class T_ARCHIVE>
              bool lns_access(T_ARCHIVE &_ar, const unsigned int)
              {
                if(_ar.is_loading())
                {
                  typename T::value_type value;
                  _ar.start_recurrence();
                  size_t i(0);
                  try
                  { 
                    while( _ar & ext(value) )
                    {
                      _ar.step_recurrence();
                      ++i;
                      if(not _ar.grammar_only()) container_.push_back(value); 
                    }
                  }
                  catch( error::section_not_found &_e)
                  { 
                    _ar.stop_recurrence();
                    if(required_ == 0u) return true;
                    if(i >= required_) return true;
                    std::string const * section_name 
                      = boost::get_error_info<error::section_name>(_e);
                    BOOST_THROW_EXCEPTION
                    ( 
                      error::too_few_sections() << error::section_name(*section_name) 
                    );
                  }
                  _ar.stop_recurrence();
                  if(i < required_)
                    BOOST_THROW_EXCEPTION( error::too_few_sections() );
                  return true;
                }
                else 
                {
                  typename T::iterator i_first = container_.begin();
                  typename T::iterator i_end = container_.end();
                  if(i_first == i_end) return true;
                  xpr::Section section(ext(*i_first));
                  for(++i_first; i_first != i_end; ++i_first)
                    section = section && ext(*i_first);  
                  return _ar & section;
                }
              }
            T &container_;
            size_t required_;
          };
//       //! Saves or loops containers with begin/end, and pushback method.
//       template<class T>
//         struct PushBack<T const>
//         {
//           PushBack(T const &_container) : container_(_container) {}
//           PushBack(PushBack<T> const &_c) : container_(_c.container_) {}
//           template<class T_ARCHIVE>
//             bool lns_access(T_ARCHIVE &_ar, const unsigned int)
//             {
//               LADA_ASSERT(not _ar.is_loading(), "Cannot load into constant container.")
//               typename  T::value_type value;
//               while(_ar & ext(value)) container_.push_back(value);
//               return true;
//             }
//           T const &container_;
//         };
      } // namespace details
    } // namespace xpr
//   template<class T> xpr::Section push_back(T const &_container)
//   {
//     boost::shared_ptr< xpr::details::PushBack<T const> >
//       ptr(new xpr::details::PushBack<T const>(_container) );
//     return ext(ptr);
//   }
    template<class T> xpr::Section push_back(T &_container, size_t _required = 0)
    {
        boost::shared_ptr< xpr::details::PushBack<T> >
          ptr(new xpr::details::PushBack<T>(_container, _required) );
        return ext(ptr);
    }
  } // namespace load_n_save
} // namespace LaDa
#endif

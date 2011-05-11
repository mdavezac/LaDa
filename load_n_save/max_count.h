//
//  Version: $Id: max_count.h 1226 2009-07-13 06:28:01Z davezac $
//

#ifndef _LADA_LOADNSAVE_MAX_COUNT_ACTION_H_
#define _LADA_LOADNSAVE_MAX_COUNT_ACTION_H_


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

namespace LaDa 
{
  namespace load_n_save
  {
    //! Pushes back into a vector.
    template<class T_LAMBDA >
      class MaxCount
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! Constructor. 
          MaxCount(size_t _i, T_LAMBDA &_op) : count_(0), max_(_i), op_(_op) {}
          //! Copy
          MaxCount(MaxCount const &_c) : count_(_c.count_), max_(_i), op_(_op) {}
          //! Does bitwise operation.
          result_type operator()( std::string const &_x ) const
          { 
            if( count_ > max_ )
            {
              std::cerr << "Too many sections or options of a kind.\n";
              return false;
            }
            if( not op_(_x) ) return false;
            ++count_;
          }

        protected:
          //! Holds reference to container.
          T_CONTAINER &container_;
      };

    //! Pushes back into a vector.
    template<class T_CONTAINER, class T_CONVERTER>
      class PushBackWithConverter : PushBack<T_CONTAINER>
      {
        public:
          //! A trait wich says this is an action.
          typedef void action;
          //! The return type is bool.
          typedef bool result_type;
          //! Constructor. 
          PushBackWithConverter   ( T_CONTAINER &_cont )
                                : PushBack<T_CONTAINER>(_cont), converter_(_converter_) {} 
          //! CopyConstructor. 
          PushBackWithConverter   ( PushBackWithConverter const &_c )
                                : PushBack<T_CONTAINER>(_c.container_), converter_(_c.converter_) {} 
          //! Does bitwise operation.
          result_type operator()( std::string const &_x ) const
          { 
            if( not initializer::parse_value( _x, converter_ ) ) return false;
            container_.push_back( converter_ );
            return true; 
          }

        protected:
          //! Holds reference to an action which converts to the container value type.
          T_CONVERTER &converter_;
      };

    //! Returns an action to push_back items into a container.
    template< class T_CONTAINER >
      PushBack<T_CONTAINER> push_back(T_CONTAINER &_cont)
        { return PushBack<T_CONTAINER>(_cont); }

    //! Returns an action to push_back items into a container.
    template< class T_CONTAINER, class T_CONVERTER >
      PushBackWithConverter<T_CONTAINER, T_CONVERTER> push_back( T_CONTAINER &_cont, T_CONVERTER& _conv )
        { return PushBackWithConverter<T_CONTAINER, T_CONVERTER>(_cont, _conv); }

  }
}

#endif

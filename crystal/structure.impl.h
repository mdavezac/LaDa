#include <boost/lambda/lambda.hpp>

namespace LaDa
{
  namespace crystal
  {
    template< class T_TYPE >
      std::ostream& operator<<(std::ostream &_stream, Structure<T_TYPE> const &_str)
        { return _stream << *_str.impl_; }

    template<class TYPE> template<class T_ARCHIVE>
      bool Structure :: lns_access(T_ARCHIVE &_ar, load_n_save::version_type const _version) 
      {
        if( _ar.is_loading() )
        {
          boost::shared_ptr< StructureData<TYPE> > dummy(impl_);
          impl_.reset(new StructureData<TYPE>());
        }
        return _ar & *impl_;
      }
    


  }
} // namespace LaDa

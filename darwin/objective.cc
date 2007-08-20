#include  "objective.h"

std::ostream & operator<<( std::ostream &_os, darwin::Fitness &_fit )
  {  return _os << (types::t_real ) _fit; }


namespace darwin
{
  bool Fitness :: Load( const TiXmlElement & _node )
  {
    if ( not _node.Attribute( "fitness" ) )
    {
      is_valid = false;
      return false; 
    }
    double d; _node.Attribute( "fitness", &d );
    quantity = (t_Quantity) d;
    is_valid = true;
    return true;
  }
  bool Fitness :: Save( TiXmlElement & _node ) const
  {
    double d = (double) quantity;
    _node.SetDoubleAttribute("fitness", d);
    return true;
  }
}

#ifdef _MPI
namespace mpi
{
  template<>
  bool BroadCast::serialize<darwin::Fitness>( darwin::Fitness &_fit )
  {
     bool result = serialize( _fit.is_valid )
            and serialize( _fit.quantity );
     if ( stage == COPYING_TO_HERE )
       std::cout <<  "Copying to here " << _fit.is_valid << " " << _fit.quantity << std::endl;
     else if ( stage == COPYING_FROM_HERE )
       std::cout <<  "Copying from here " << _fit.is_valid << " " << _fit.quantity << std::endl;
     return result;
  }
//   if ( stage == GETTING_SIZE )
//   {
//     ++buffer_size[0];
//     if ( _fit.is_valid )  ++buffer_size[2];
//     return true;
//   }
//   else if(    end_real_buff - cur_real_buff < 1 
//            or end_int_buff  - cur_int_buff < 1 )
//     return false;
//   else if ( stage == COPYING_TO_HERE )
//   {
//     *cur_int_buff = ( _fit.is_valid ) ? 1 : 0;
//     ++cur_int_buff;
//     if ( not _fit.is_valid ) return true;
//     *cur_real_buff = _fit.quantity;
//     ++cur_real_buff;
//     return true;
//   }
//   else if ( not stage == COPYING_FROM_HERE )
//     return false;
//
//   // COPYING FROM HERE
//   _fit.is_valid = ( *cur_int_buff != 0 );
//   ++cur_int_buff;
//   if ( not _fit.is_valid ) return true;
//   _fit.quantity = *cur_real_buff;
//   ++cur_real_buff;
//    
//   return true;   
// }
} // namespace mpi
#endif 


#include "gatraits.h"

namespace Traits 
{
    inline template<> void zero_out<types::t_real>( types::t_real &_cont )
      { _cont = types::t_real(0);  }
    inline template<> void zero_out<types::t_int>( types::t_int &_cont )
      { _cont = types::t_int(0);  }
    inline template<> void zero_out<types::t_unsigned>( types::t_unsigned &_cont )
      { _cont = types::t_unsigned(0);  }
    inline template<> void zero_out<bool>( bool &_cont )
      { _cont = false;  }
    inline template<> void zero_out<std::string>( std::string &_cont )
      { _cont = "";  }
}

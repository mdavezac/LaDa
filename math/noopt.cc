#include "PyladaConfig.h"

#include "math.h"

namespace Pylada
{
  namespace math
  {
    namespace details
    {
      bool no_opt_change_test(types::t_real _new, types::t_real _last)
      {
        types::t_real const m_new = 16e0 * _new;
        types::t_real const diff = _new - _last;
        types::t_real const m_new_plus_diff = m_new + diff;
        types::t_real const m_new_plus_diff_minus_m_new = m_new_plus_diff - m_new;
        return m_new_plus_diff_minus_m_new != 0;
      }
    }
  }
}

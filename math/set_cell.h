#ifndef LADA_MATH_SETCELL_H
#define LADA_MATH_SETCELL_H
#include "LaDaConfig.h"

#include <boost/mpl/int.hpp>
#include "eigen.h"

namespace LaDa
{
  namespace math
  {
    namespace details
    {
      //! Easily sets the cell parameters of a StructureData and TemplateData.
      template<class D = boost::mpl::int_<0> >
        class SetCell
        {
          public:
            SetCell(math::rMatrix3d &_cell) : cell_(_cell) {}
            SetCell<typename boost::mpl::next<D>::type>
              operator()(types::t_real _x, types::t_real _y, types::t_real _z)
              {
                cell_(D::value, 0) = _x;
                cell_(D::value, 1) = _y;
                cell_(D::value, 2) = _z;
                return SetCell<typename boost::mpl::next<D>::type>(cell_);
              }
            SetCell<typename boost::mpl::next<D>::type>
              operator()(math::rVector3d const &_pos)
              {
                cell_(D::value, 0) = _pos[0];
                cell_(D::value, 1) = _pos[1];
                cell_(D::value, 2) = _pos[2];
                return SetCell<typename boost::mpl::next<D>::type>(cell_);
              }
          private:
            math::rMatrix3d &cell_;
        };
      //! \brief Easily sets the cell parameters of a StructureData and TemplateData.
      //! \details This one returns void, so that only three instances of
      //!          SetCell can be called.
      template<> class SetCell< boost::mpl::int_<2> >
      {
        public:
          SetCell(math::rMatrix3d &_cell) : cell_(_cell) {}
          void operator()(types::t_real _x, types::t_real _y, types::t_real _z)
          {
            cell_(boost::mpl::int_<2>::value, 0) = _x;
            cell_(boost::mpl::int_<2>::value, 1) = _y;
            cell_(boost::mpl::int_<2>::value, 2) = _z;
          }
          void operator()(math::rVector3d const &_pos)
          {
            cell_(boost::mpl::int_<2>::value, 0) = _pos[0];
            cell_(boost::mpl::int_<2>::value, 1) = _pos[1];
            cell_(boost::mpl::int_<2>::value, 2) = _pos[2];
          }
        private:
          math::rMatrix3d &cell_;
      };
    }
  }
}
#endif 

      double precision function df1dim (func,x,dx)
!
! Evaluate function func and its directional derivative dx at
! a point specified by the point pcom (in common), the direction
! xicom (in common), and the parameter x, which specifies the
! distance along xicom from pcom.  func must calculate its
! gradient as a side effect.
!
! $Id: df1dim.fpp,v 1.1.1.1 2004/02/13 21:44:24 wjones Exp $
!
! Input:
      double precision func,x
      external func
! Output:
      double precision dx
! Parameters:
      integer mxparm
      parameter (mxparm = 135)
! Common:
      integer ncom
      double precision pcom,xicom
      common /f1com/ ncom,pcom(mxparm),xicom(mxparm)
! Local:
      integer j
      double precision xt(mxparm),df(mxparm)
      common/df1dim_com/xt,df

! Constants:
      double precision zero
      parameter (zero = 0.0d0)
!
! Find the n-dimensional point:
!

      do 11 j = 1,ncom
         xt(j) = pcom(j) + x * xicom(j)
   11 continue
!
! Evaluate the function and its gradient:
!
      df1dim = func (xt,df)
!
! Determine directional derivative:
!
      dx = zero
      do 12 j = 1,ncom
         dx = dx + df(j) * xicom(j)
   12 continue
!
      return
      end

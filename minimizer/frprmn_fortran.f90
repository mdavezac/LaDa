      subroutine frprmn (p,n,ftol,ltol,zeps,iter,fret,func,itmax,rtol)
!
! From William H. Press, Brian P. Flannery, Saul A. Teukolsky, and
! William T. Vetterling, Numerical Recipes (Cambridge, New York, 1986),
! pp. 305-306.
!
! Given a starting point p that is a vector of length n, Fletcher-
! Reeves-Polak_Ribiere minimization is performed on a function func,
! using its gradient.  The convergence tolerance on the function value
! is input as ftol.  Returned quantities are p (the location of the
! minimum), iter (the number of iterations that were performed), and
! fret (the minimum value of the function).  The routine linmin is
! called to perform line minimizations.
!
! $Id: frprmn.fpp,v 1.1.1.1 2004/02/13 21:44:24 wjones Exp $
!
! Modifications:
! * func is now a dummy function name (used here and passed to linmin)
! * Caller specifies maximum number of iterations (itmax), and is
!   responsible for verifying that convergence was achieved
! * Achieved tolerance is returned (rtol)
! * Put ltol and zeps (for dbrent) into arg list
!
! Input:
      integer n,itmax
      double precision ftol,ltol,zeps,func
      external func
! Output:
      integer iter
      double precision fret,rtol
! Input and output:
      double precision p(n)
! Parameters:
      double precision eps
      parameter (eps = 1.0d-10)
! Local:
      integer j,its
      double precision g(n),h(n),xi(n),fp,fpl,gg,dgg,gam
!     common/frpmn_com/g,h,xi
! Constants:
      double precision zero,two
      parameter (zero = 0.0d0, two = 2.0d0)
      integer uout
      common /io/ uout
      double precision xmin
      common /myxmin/xmin

      xmin = 0d0

!
! Initializations:
!

      fp = func (p,xi)
      do 11 j = 1,n
         g(j) = -xi(j)
         h(j) = g(j)
         xi(j) = h(j)
   11 continue
!
! Loop over iterations:
!
      do 14 its = 1,itmax
         iter = its
         call linmin (p,xi,n,ltol,zeps,func,fret)
!PK	 write(*,*)"iter residual criterion",
!PK     $        iter,two*abs(fret-fp),ftol*(abs(fret)+abs(fp)+eps)
         if (two*abs(fret-fp) .le. ftol*(abs(fret)+abs(fp)+eps)) then
            rtol = two * abs(fret-fp) / (abs(fret)+abs(fp)+eps)
            return
         endif
         fpl = fp
         fp = func (p,xi)
         gg = zero
         dgg = zero
         do 12 j = 1,n
            gg = gg + g(j)**2
! This statement for Fletcher-Reeves:
!           dgg = dgg + xi(j)**2
! This statement for Polak-Ribiere:
            dgg = dgg + (xi(j)+g(j)) * xi(j)
   12    continue
         if (gg .eq. zero) then
           return
         endif
         gam = dgg / gg
         do 13 j = 1,n
            g(j) = -xi(j)
            h(j) = g(j) + gam * h(j)
            xi(j) = h(j)
   13    continue
   14 continue
!
! Return achieved tolerance upon convergence failure:
!
      rtol = two * abs(fret-fpl) / (abs(fret)+abs(fpl)+eps)
!
      return
      end


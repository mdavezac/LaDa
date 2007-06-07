      double precision function dbrent (ax,bx,cx,f,tol,zeps,func,xmin)
!
! From William H. Press, Brian P. Flannery, Saul A. Teukolsky, and
! William T. Vetterling, Numerical Recipes (Cambridge, New York, 1986),
! pp. 287-289.
!
! Given a function f that calculates its derivative as a side effect,
! and given a bracketing triplet of abscissas ax, bx, cx [such that bx
! is between ax and cx, and f(bx) is less than both f(ax) and f(cx)],
! this routine isolates the minimum to a fractional precision of
! about tol using a modification of Brent's method that uses
! derivatives.  The abscissa of the minimum is returned as xmin, and
! the minimum function value is returned as dbrent, the returned
! function value.
!
! $Id: dbrent.fpp,v 1.1.1.1 2004/02/13 21:44:24 wjones Exp $
!
! Modifications:
! * func is a dummy function that is passed to f -- it does the
!   n-dimensional evaluations needed by f to determine the function
!   value and directional derivative
! * Eliminated the function df -- the derivative is computed by f
! * Put zeps into arg list
!
! Input:
      double precision ax,bx,cx,f,tol,zeps,func
      external f,func
! Output:
      double precision xmin
! Parameters:
      integer itmax
      parameter (itmax = 100)
! Local:
      integer iter
      double precision a,b,v,w,x,e,fx,fv,fw,dx,dv,dw,xm,tol1,tol2
      double precision d1,d2,u1,u2,olde,d,u,fu,du
      logical ok1,ok2
#ifdef SYSVAT
      double precision xtmp
#endif
! Constants:
      double precision zero,half,two
      parameter (zero = 0.0d0, half = 0.5d0, two = 2.0d0)
!
      a = min(ax,cx)
      b = max(ax,cx)
      v = bx
      w = v
      x = v
      e = zero
      fx = f(func,x,dx)
      fv = fx
      fw = fx
      dv = dx
      dw = dx
      do 11 iter = 1,itmax
         xm = half * (a+b)
         tol1 = tol*abs(x) + zeps
         tol2 = two * tol1
#ifndef SYSVAT
         if (abs(x-xm) .le. (tol2-half*(b-a))) goto 3
#else
         xtmp = tol2 - half*(b-a)
         if (abs(x-xm) .le. xtmp) goto 3
#endif
         if (abs(e) .gt. tol1) then
            d1 = two * (b-a)
            d2 = d1
            if (dw .ne. dx) d1 = (w-x) * dx / (dx-dw)
            if (dv .ne. dx) d2 = (v-x) * dx / (dx-dv)
            u1 = x + d1
            u2 = x + d2
            ok1 = ((a-u1)*(u1-b) .gt. zero) .and. (dx*d1 .le. zero)
            ok2 = ((a-u2)*(u2-b) .gt. zero) .and. (dx*d2 .le. zero)
            olde = e
            e = d
            if (.not. (ok1 .or. ok2)) then
               goto 1
            else if (ok1 .and. ok2) then
               if (abs(d1) .lt. abs(d2)) then
                  d = d1
               else
                  d = d2
               endif
            else if (ok1) then
               d = d1
            else
               d = d2
            endif
            if (abs(d) .gt. abs(half*olde)) goto 1
            u = x + d
            if ((u-a .lt. tol2) .or. (b-u .lt. tol2)) d = sign(tol1,xm-x)
            goto 2
         endif
    1    if (dx .ge. zero) then
            e = a - x
         else
            e = b - x
         endif
         d = half * e
    2    if (abs(d) .ge. tol1) then
            u = x + d
            fu = f(func,u,du)
         else
            u = x + sign(tol1,d)
            fu = f(func,u,du)
            if (fu .gt. fx) goto 3
         endif
         if (fu .le. fx) then
            if (u .ge. x) then
               a = x
            else
               b = x
            endif
            v = w
            fv = fw
            dv = dw
            w = x
            fw = fx
            dw = dx
            x = u
            fx = fu
            dx = du
         else
            if (u .lt. x) then
               a = u
            else
               b = u
            endif
            if ((fu .le. fw) .or. (w .eq. x)) then
               v = w
               fv = fw
               dv = dw
               w = u
               fw = fu
               dw = du
            else if ((fu .le. fv) .or. (v .eq. x) .or. (v .eq. w)) then
               v = u
               fv = fu
               dv = du
            endif
         endif
   11 continue
      pause 'dbrent exceeded maximum iterations.'
    3 xmin = x
      dbrent = fx
!
      return
      end

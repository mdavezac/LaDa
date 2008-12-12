      subroutine linmin (p,xi,n,tol,zeps,func,fret)
!
!c  rewritten by LW Wang, don't use the stupid mnbrack.x and dbrent.x routines
!c
! From William H. Press, Brian P. Flannery, Saul A. Teukolsky, and
! William T. Vetterling, Numerical Recipes (Cambridge, New York, 1986),
! pp. 300-301.
!
! Given an n dimensional point p and an n dimensional direction xi,
! moves and resets p to where the function func(p) takes on a minimum
! along the direction xi from p, and replaces xi by the actual vector
! displacement that p was moved.  Also returns as fret the value of
! func at the returned location p.  This is actually all accomplished
! by calling the routines mnbrak and brent (or dbrent, if using
! derivatives).
!
! $Id: linmin.fpp,v 1.1.1.1 2004/02/13 21:44:24 wjones Exp $
!
! Modifications:
! * func is a dummy function that is passed to mnbrak and dbrent to
!   do the n-dimensional evaluations needed for the line minimization
! * Put tol (for dbrent) into arg list
! * Put zeps (for dbrent) into arg list
!
! Input:
      integer n
      double precision tol,zeps,func
      double precision x,dx,x1,dx1,x2,dx2,a,b
      integer i
      external func
!      save xmin
! Input and output:
      double precision p(n),xi(n)
! Output:
      double precision fret,fret1,fret2,fret3
! Parameters:
      integer mxparm
      parameter (mxparm = 300)
! Local:
      integer j
!      double precision ax,xx,bx,fa,fx,fb,xmin
      double precision ax,xx,bx,fa,fx,fb
! External:
      double precision df1dim
      external df1dim
! Common:
	double precision xmin
	common /myxmin/xmin
      integer ncom
      double precision pcom,xicom
      common /f1com/ ncom,pcom(mxparm),xicom(mxparm)
! Constants:
      double precision zero,one,two
      parameter (zero = 0.0d0, one = 1.0d0, two = 2.0d0)
      double precision shit
      integer uout
      common /io/ uout
!

!cccccccc pcom(j) is the current position
!cccccccc xicom(j) is the searching direction

      ncom = n
      do 11 j = 1,n
         pcom(j) = p(j)
         xicom(j) = xi(j)
   11 continue
!
! Initial guess for brackets:
!
!cccccccccccccc old version, it is very expensive, and unstable
!c      ax = zero
!c      xx = one    ! far too larg, break down for higher order effects 
!c      bx = two
!c      xx = 0.000005
!c      bx = 0.00001
!c      call mnbrak (ax,xx,bx,fa,fx,fb,df1dim,func)
!c      fret = dbrent (ax,xx,bx,df1dim,tol,zeps,func,xmin)
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Construct the vector results to return:
!

!ccccc using simple quadratic interpolation.
!cccc  E=E0+b*x+a*x**2

!cccc xmin is not assigned values at the first iter, 
!cccc so this might break down for
!cccc other machines.

!       do i=-200,200
!       x1=i/200.d0*0.01
!       fret1=df1dim(func,x1,dx1)
!       write(6,*) x1,fret1
!       enddo
!       stop

       x1=0.d00
       fret1=df1dim(func,x1,dx1)
       x2=xmin
       if(x2.lt.1.D-10.or.x2.gt.0.1) x2=0.0001      ! starting value
       fret2=df1dim(func,x2,dx2)
       b=dx1
       a=(fret2-fret1-b*x2)/x2**2
! Protect against constant expressions.
! Check on exact ZERO - no numeric tolerances
! PK 19991018
       if (a.eq.0.) then
          xmin=fret1
          fret3=fret1
          fret=df1dim(func,xmin,dx)
       else	
          xmin=-b/(2*a)
          fret3=fret1+b*xmin+a*xmin**2
          fret=df1dim(func,xmin,dx)
       end if
!PK       write(6,*) "xmin,E_pred,E",xmin,fret3,fret


      do 12 j = 1,n
         shit = p(j)
         xi(j) = xmin * xi(j)
         p(j) = p(j) + xi(j)
   12 continue
!
      return
      end

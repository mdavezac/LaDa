MODULE rational_mathematics
  implicit none
  private
  public gcd, SmithNormalFormInternal

  ! Overloaded procedure for computing the greatest common denominator
  INTERFACE gcd
     MODULE PROCEDURE gcd_2ints, gcd_rank1, gcd_3ints, gcd_4ints
  END INTERFACE

  CONTAINS
  !*******************************************************************************
  ! This routine takes an integer 3x3 matrix and computes its Smith Normal Form
  subroutine SmithNormalFormInternal(H,A,M,B)
    integer, intent(in) :: H(3,3) ! Input matrix
    ! Smith normal form (M), left (A) & right (B) transforms
    integer, intent(out), dimension(3,3) :: M, A, B 

    integer i, minm, maxidx, minidx, multiple, j, nondividx(2), check(9)
    logical Ldiv(2,2)

    A = 0; B = 0; M = H ! M starts out as H, the input matrix
    forall(i=1:3); A(i,i) = 1; B(i,i) = 1; end forall ! A & B = identity

    j = 1 ! Keep track of which row/col we are on
    do ! Keep doing steps (1) and (2) until all elements in j-th row and column are
       ! zero except the one on the diagonal. When j == 1, if the (1,1) element doesn't
       ! divide every other non-zero element left in the matrix at the end, then add the
       ! offending row to row 1 and try again.
    !print *,"J=",j
       !(1) Use row operations to zero out the column
       ! Divide the row with the smallest value into the largest
       !print *,"row operations"
       do while (count(M(:,j)/=0) > 1) ! Keep going until only 1 non-zero element in column j
          !call printAMB(A,M,B)
          call get_minmax_indices(M(:,j),minidx,maxidx)
          minm = M(minidx,j)
          !print *, "Min:",minm,"maxindex",maxidx,"minindex",minidx
          ! Subtract a multiple of the row containing the smallest element from
          ! the row containing the largest element
          multiple = M(maxidx,j)/minm
          M(maxidx,:) = M(maxidx,:) - multiple*M(minidx,:)
          A(maxidx,:) = A(maxidx,:) - multiple*A(minidx,:)
          if (any(matmul(matmul(A,H),B)/=M)) stop "ROW: Transformation matrices didn't work"
          !read(*,*)
       enddo ! End of step 1
       !call printAMB(A,M,B)
       !print *, "End of step 1"
          
       if (M(j,j) == 0) then; call swap_row(A,M,j)
          !print *, "row swap"
          !call printAMB(A,M,B);
       endif
       if (any(matmul(matmul(A,H),B)/=M)) stop "ROWSWAP: Transformation matrices didn't work"
       
       !print *,"column operations"
       !(2) Use column operations to zero out first row
       ! Divide the colum with the smallest value into the largest
       do while (count(M(j,:)/=0) > 1) ! Keep going until only 1 non-zero element in row 1
          !call printAMB(A,M,B)
          call get_minmax_indices(M(j,:),minidx,maxidx)
          minm = M(j,minidx)
          ! Subtract a multiple of the column containing the smallest element from
          ! the row containing the largest element
          multiple = M(j,maxidx)/minm ! Factor to multiply by
    !      !print *, "Min:",minm,"maxindex",maxidx,"minindex",minidx
          M(:,maxidx) = M(:,maxidx) - multiple*M(:,minidx)
          B(:,maxidx) = B(:,maxidx) - multiple*B(:,minidx)
          if (any(matmul(matmul(A,H),B)/=M)) stop "COLS: Transformation matrices didn't work"
          !read(*,*)
       enddo ! End of step 2
       !call printAMB(A,M,B)
       !print *, "End of step 2"
       if (M(j,j)<0) then ! Change signs
          M(:,j) = -M(:,j); B(:,j) = -B(:,j)
       elseif (M(j,j) == 0) then;call swap_column(M,B,j)
          !print *, "column swap"
          !call printAMB(A,M,B)
       endif
       if (count(M(j,:)/=0) > 1 .or. count(M(:,j)/=0) > 1) cycle

       if (any(matmul(matmul(A,H),B)/=M)) stop "COLSWAP: Transformation matrices didn't work"

       Ldiv = mod(M(2:,2:),M(1,1)) == 0
       if (j==1 .and. any(Ldiv .eqv. .false.)) then! Add the offending row to row 1 
          nondividx = maxloc(mod(M(2:,2:),M(1,1))) ! Find one of the elements that isn't 
          M(1,:) = M(1,:) + M(nondividx(1)+1,:)    ! divided by the diagonal element, and 
          A(1,:) = A(1,:) + A(nondividx(1)+1,:)    ! add the row it's in to row 1
          cycle ! Go back to the top of the outer do loop
       endif
       if (j==2) then 
          if (mod(M(3,3),M(2,2))/=0) then
             M(2,:) = M(2,:) + M(3,:)
             A(2,:) = A(2,:) + A(3,:)
             cycle
          endif
       else
          j = 2;cycle;endif ! Start row/col 2
       ! Try again if the matrix still isn't diagonal
       if (j == 2 .and. (M(3,2)/=0 .or. M(2,3)/=0)) then; cycle; endif
       exit ! We should be done if we hit this point
    enddo
    if (M(3,3)<0) then ! Change signs
       M(:,3) = -M(:,3); B(:,3) = -B(:,3);
       !call printAMB(A,M,B)
    endif
    if (any(matmul(matmul(A,H),B)/=M)) stop "END: Transformation matrices didn't work"
    check = reshape(M,(/9/))
    if (any(check((/2,3,4,6,7,8/))/=0)) stop "Not diagonal"
    if (mod(M(2,2),M(1,1))/=0 .or. mod(M(3,3),M(2,2))/=0) stop "SNF conditions not met"

    contains
      subroutine swap_row(A,M,k) ! Swap rows of M (and A)
        integer, intent(inout) :: M(3,3), A(3,3), k
        integer tmpRow(3), maxidx(1)
         
        maxidx = maxloc(abs(M(k:,k)))+k-1  ! find index of the non-zero element in col k
        tmpRow = A(k,:); A(k,:) = A(maxidx(1),:); A(maxidx(1),:) = tmpRow
        tmpRow = M(k,:); M(k,:) = M(maxidx(1),:); M(maxidx(1),:) = tmpRow
      endsubroutine swap_row

      subroutine swap_column(M,B,k) ! Swap columns of M (and B)
        integer, intent(inout) :: M(3,3), B(3,3),k
        integer tmpCol(3), maxidx(1)
         
        maxidx = maxloc(abs(M(k,k:)))+k-1 ! find index of the non-zero element in row k
        tmpCol = B(:,k); B(:,k) = B(:,maxidx(1)); B(:,maxidx(1)) = tmpCol
        tmpCol = M(:,k); M(:,k) = M(:,maxidx(1)); M(:,maxidx(1)) = tmpCol
      endsubroutine swap_column

      subroutine printAMB(A,M,B)
        integer, intent(in), dimension(3,3) :: A,M,B
        integer i
        do i = 1,3
           write(*,'(3(3x,3i4))') A(i,:),M(i,:),B(i,:);enddo;print *
      endsubroutine printAMB

      subroutine get_minmax_indices(invec,min,max)
        integer, intent(in) :: invec(3)
        integer, intent(out) :: min, max

        integer :: tmpmin(1), tmpmax(1), vec(3)
        vec = abs(invec)
        tmpmin = minloc(vec,vec>0)
        ! Search from the right for the max so it will be different from the minimum
        ! even if the min and max are the same value
        tmpmax = 4 - maxloc(vec(3:1:-1),vec(3:1:-1)>0)
        min = tmpmin(1)
        max = tmpmax(1)
      endsubroutine get_minmax_indices
  endsubroutine SmithNormalFormInternal

  !*****************************************************************************
  ! This function finds the greatest common denominator of several integers
  ! ****************************************************************************
  ! 
  ! This case works for two integers, given as separate arguments
  function gcd_2ints(x1, x2) result(divisor)
    integer, intent(in) :: x1, x2
    integer divisor

    integer a, b
    a = abs(x1); b = abs(x2) ! Make sure inputs are positive
    if (b>a) call swap(a,b)  

    do ! Keep dividing a by b, until one of them is zero
       if (b>a) call swap(a,b) ! Keep the bigger number in a's place
       if (b == 0) exit ! we're done when b == 0
       a = mod(a,b) ! Divide a by b and keep only the remainder
    enddo
    divisor = a

    contains
      subroutine swap(x,y) ! Swap two values
        integer x,y,tmp
        tmp = x; x = y; y = tmp
      endsubroutine swap
  end function gcd_2ints

  ! This case works on a list of integers (a rank-1 array)
  function gcd_rank1(x) result(divisor)
    integer, intent(in) :: x(:)
    integer divisor

    integer a(size(x)), N, indx(1), big2

    N = size(x); a = abs(x)
    if (any(a<1)) stop "GCD requires positive integers"
    do ! Divide the biggest number by the second biggest until 
       ! the second biggest is zero
       indx = maxloc(a)  ! Find the location of the biggest number
       if (all(a == a(indx(1)))) then ! check if all numbers are the same
          big2 = a(indx(1))
       else   ! The "real" around 'a' is a workaround for a problem in the Absoft compiler
          big2 = maxval(real(a),mask=a < a(indx(1))) ! Find the size of the 2nd biggest number
       endif
       if (big2 == 0) exit
       a(indx(1)) = mod(a(indx(1)),big2)
    enddo
    divisor = a(indx(1)) ! The divisor is the number left when every other member==0
  endfunction gcd_rank1

  ! This case works on 3 integers, not in an array
  function gcd_3ints(x1,x2,x3)
    integer, intent(in) :: x1,x2,x3
    integer gcd_3ints
    gcd_3ints = gcd_rank1((/x1,x2,x3/))
  end function gcd_3ints

  ! This case works on 4 integers, not in an array
  function gcd_4ints(x1,x2,x3,x4)
    integer, intent(in) :: x1,x2,x3,x4
    integer gcd_4ints
    gcd_4ints = gcd_rank1((/x1,x2,x3,x4/))
  end function gcd_4ints

END MODULE rational_mathematics

  subroutine SmithNormalForm(H,A,M,B)
    
    use rational_mathematics, only: SmithNormalFormInternal
    integer, intent(in) :: H(3,3) 
    integer, intent(out), dimension(3,3) :: M, A, B 

    call SmithNormalFormInternal( H, A, M, B )

  end subroutine

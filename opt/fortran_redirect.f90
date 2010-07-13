! Redirects stdout, stderr, or stdin
subroutine lada_redirect_open(which, N, path, isopen, doappend)
     use iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT, INPUT_UNIT
     integer, intent(in) :: which ! 0 => stderr, 5 => stdin, 6 => stdout
     integer, intent(in) :: N     ! size of the path.
     character(len=N), intent(in) :: path ! path name
     integer, intent(out) :: isopen ! 1 if no error.
     integer, intent(in) :: doappend ! 1 if append to file. 
     
     isopen = 1
     if(doappend .eq. 1) then
       if( which == 0 ) then ! should be using iso_fortran_env's OUTPUT_UNIT
         open(unit=ERROR_UNIT, file=path, action='WRITE', err=100, POSITION='APPEND')
       else if( which == 5 ) then ! should be using iso_fortran_env's INPUT_UNIT
         open(unit=INPUT_UNIT, file=path, action='READ', err=100, POSITION='APPEND')
       else if( which == 6 ) then ! should be using iso_fortran_env's ERROR_UNIT
         open(unit=OUTPUT_UNIT, file=path, action='WRITE', err=100, POSITION='APPEND')
       endif
     else
       if( which == 0 ) then ! should be using iso_fortran_env's OUTPUT_UNIT
         open(unit=ERROR_UNIT, file=path, action='WRITE', err=100, STATUS='REPLACE')
       else if( which == 5 ) then ! should be using iso_fortran_env's INPUT_UNIT
         open(unit=INPUT_UNIT, file=path, action='READ', err=100, STATUS='REPLACE')
       else if( which == 6 ) then ! should be using iso_fortran_env's ERROR_UNIT
         open(unit=OUTPUT_UNIT, file=path, action='WRITE', err=100, STATUS='REPLACE')
       endif
     endif

     if( isopen .eq. 1 ) return

100  isopen = 0
     
end subroutine

subroutine lada_redirect_close(which)
     use iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT, INPUT_UNIT
     integer, intent(in) :: which ! 0 => stderr, 5 => stdin, 6 => stdout
     
     if( which == 0 ) then 
       close(ERROR_UNIT)
       open(unit=ERROR_UNIT, file="/dev/stderr", POSITION='ASIS')
     else if( which == 5 ) then 
       close(INPUT_UNIT)
       open(unit=INPUT_UNIT, file="/dev/stdin" , POSITION='ASIS')
     else if( which == 6 ) then 
       close(ERROR_UNIT)
       open(unit=ERROR_UNIT, file="/dev/stdout", POSITION='ASIS')
     endif
     ! OpenVMS:
     ! Open (Unit=5, File="SYS$INPUT") ! Open Standard Input on Unit 5
     ! Open (Unit=6, File="SYS$OUTPUT") ! Open Standard Output on Unit 6
     ! Open (Unit=0, File="SYS$ERROR") ! Open Standard Error on Unit 0

     ! Windows
     ! Open (Unit=5, File="stdin") ! Open Standard Input on Unit 5
     ! Open (Unit=6, File="stdout") ! Open Standard Output on Unit 6
     ! Open (Unit=0, File="stderr") ! Open Standard Error on Unit 0
end subroutine

! permissions are incorrect when using aprun.
subroutine fucking_cray(which) 
     use iso_fortran_env, only: OUTPUT_UNIT, ERROR_UNIT, INPUT_UNIT
     integer, intent(in) :: which ! 0 => stderr, 5 => stdin, 6 => stdout
     
     if( which == 0 ) then 
       close(ERROR_UNIT)
     else if( which == 5 ) then 
       close(INPUT_UNIT)
     else if( which == 6 ) then 
       close(OUTPUT_UNIT)
     endif
end subroutine

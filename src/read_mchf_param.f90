!====================================================================
    SUBROUTINE read_mchf_param
!====================================================================
!   Read mchf_param if it is connected (exists)
!--------------------------------------------------------------------

    USE mchf_inout
    USE mchf_param
    USE cl_arguments

    IMPLICIT NONE
    INTEGER       :: ios
    CHARACTER(LEN=80) :: line, line1, line2
    CHARACTER(LEN=7)  :: var
    CHARACTER(LEN=20) :: value
    Logical           :: connected

    INQUIRE(FILE='mchf_param', OPENED=connected)
    !Print *, 'Is mchf_param connected?',connected
    If (.not. connected) return

    Read(iup,'(A)')  line1
    line1 = adjustl(line1)
    If (line1(1:4) /= 'MCHF') &
       Stop 'First line of mchf_param does not start with MCHF'
    Read(iup,'(A)')  line2
    
    Do 
      call get_next
      if (ios /= 0) EXIT
      select case (trim(var))
       case ('scf_tol'); Read(value,*) scf_tol
       case ('orb_tol'); Read(value,*) orb_tol
       case ('end_tol'); Read(value,*) end_tol 
       case ('acc_max'); Read(value,*) acc_max
       case ('varied1'); Read(value,*) varied1
       case ('varied2'); Read(value,*) varied2
       case ('list'   ); READ(value, FMT='(A)') vstring
       case ('etl_sum'); Read(value, *) etl_sum
      end select
    END DO
    !Print *, 'read mchf_param:', scf_tol, orb_tol, end_tol, acc_max, varied

    if (CL_scf_tol > 0    ) scf_tol = CL_scf_tol
    if (CL_orb_tol > 0    ) orb_tol = CL_orb_tol
    if (CL_end_tol > 0    ) end_tol = CL_end_tol
    if (CL_acc_max > 0    ) acc_max = CL_acc_max
    if (CL_varied1 /= '-1') varied1 = CL_varied1
    if (CL_varied2 /= '-1') varied2 = CL_varied2

   CONTAINS
!======================================================================
   SUBROUTINE get_next
!======================================================================
   IMPLICIT NONE

   INTEGER ::  i,ii

   Read(iup,'(A)',IOSTAT=ios) line
   IF (ios /= 0) Return
   i = INDEX(line,'=')  
   var =trim(adjustl( line(1:i-1)))
   value = adjustl(line(i+1:))
   
  END subroutine get_next

  END subroutine read_mchf_param
    
!====================================================================
    SUBROUTINE write_mchf_param
!====================================================================
!   Write mchf_param 
!--------------------------------------------------------------------
    USE spline_param
    USE mchf_inout
    USE mchf_param
    IMPLICIT NONE
    INTEGER       :: ios, nvar
    Logical           :: connected
    CHARACTER(LEN=3),DIMENSION(25) :: elc


    INQUIRE(FILE='mchf_param', OPENED=connected)
    If (.not. connected) then
      OPEN(UNIT=iup, FILE='mchf_param', STATUS='UNKNOWN', &
                     FORM='FORMATTED', ACTION='READWRITE')
    Else
        rewind(unit=iup)
    End if
    WRITE(UNIT=iup,FMT='(A)') 'MCHF_parameters'
    WRITE(UNIT=iup,FMT='(A)') '---------------'
    WRITE(UNIT=iup,FMT='(A,F12.6)') 'h       = ',h
    WRITE(UNIT=iup,FMT='(A,I12)') 'ks      = ',ks
    WRITE(UNIT=iup,FMT='(A,I12)') 'ns      = ',ns
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'scf_tol = ',scf_tol
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'orb_tol = ',orb_tol
    WRITE(UNIT=iup,FMT='(A,1PD12.2)') 'end_tol = ',end_tol
    WRITE(UNIT=iup,FMT='(A,I12)') 'acc_max = ',  1
    WRITE(UNIT=iup,FMT='(A,8X,A4)') 'varied1 = ', adjustr(varied1)
    WRITE(UNIT=iup,FMT='(A,8X,A4)') 'varied2 = ', adjustr(varied2)
   !If (varied == 'list') &
   !  WRITE(UNIT=iup,FMT='(A,A)')   'list    = ', trim(vstring)
    WRITE(UNIT=iup,FMT='(A,1PD16.9)') 'etl_sum = ', etl_sum
    
    
    END subroutine write_mchf_param 
!====================================================================

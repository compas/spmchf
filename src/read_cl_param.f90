!====================================================================
    SUBROUTINE  read_CL_param
!====================================================================
!   Read mchf_param if it is connected (exists)
!--------------------------------------------------------------------
    USE CL_arguments

    IMPLICIT NONE
    INTEGER       :: i, iarg, ieq
    CHARACTER(LEN=80) :: arg
    CHARACTER(LEN=7)  :: var
    CHARACTER(LEN=20) :: value

   ! initialize to -1
    CL_scf_tol=-1; CL_orb_tol=-1; CL_cfg_tol=-1; CL_end_tol=-1
    CL_Z      =-1; CL_acc_max=-1; CL_max_it1=-1; CL_max_it2=-1
    CL_h      =-1; CL_ns     =-1; CL_ks     =-1


    iarg = command_argument_count()
    If (iarg == 0) return

    Do  i = 1, iarg

      call get_command_argument(i, arg)

      if (len_trim(arg) == 0) cycle

      ieq = INDEX(arg,'=')
      var =trim(adjustl(arg(1:ieq-1)))
      value = adjustl(arg(ieq+1:))

      select case (trim(var))
        case ('atom'   ); Read(value,*) CL_atom
        case ('z'      ); Read(value,*) CL_z
        case ('h'      ); Read(value,*) CL_h
        case ('ns'     ); Read(value,*) CL_ns
        case ('ks'     ); Read(value,*) CL_ks
        case ('scf_tol'); Read(value,*) CL_scf_tol
        case ('orb_tol'); Read(value,*) CL_orb_tol
        case ('end_tol'); Read(value,*) CL_end_tol 
        case ('acc_max'); Read(value,*) CL_acc_max
        case ('max_it1'); Read(value,*) CL_max_it1
        case ('max_it2'); Read(value,*) CL_max_it2
        case ('varied1'); Read(value,*) CL_varied1
        case ('varied2'); Read(value,*) CL_varied2
      end select
    END DO

  END subroutine read_CL_param
    

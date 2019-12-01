!======================================================================
  PROGRAM Spline_MCHF
!======================================================================
!
!  This program computes the radial functions for simple Hartree-Fock
!  cases.
!
!   SUBROUTINE called:
!       open_files
!       read_mchf_param
!       get_case
!       define_spline
!       get_estimates
!       refine_grid
!       write_spline_param
!       solve_MCHF_equations
!       write_mchf_param
!
!   Written by:  by Charlotte  Froese Fischer
!   Date:        April, 2010 
!                                                                
!----------------------------------------------------------------------
!
    USE mchf_inout
    USE mchf_param
    IMPLICIT NONE
   

    ! read Command-Line arguments
    CALL read_CL_param
    
    CALL open_files

    ! .. read parameters (optional) or set defaults
    CALL read_mchf_param

    ! .. get data about the problem to be solved and how the spline
    ! .. calculation is to be performed, including the grid points
    CALL get_case

    ! .. Initialize the basic B-spline environment and
    ! .. evaluate spline matrix elements for typical operators in
    ! .. in non-relativistic atomic structure theory using the
    ! .. spline-galerkin method
    CALL define_spline

    ! Compute for good and excellent accuracy
    Do acc = 1,acc_max
      If (acc==1) then
         Call get_estimates
      Else
         Call refine_grid
      End if
      Call write_spline_param
      Call solve_MCHF_equations

    End DO
    Call write_mchf_param        ! for rerunning the highest level
    
  END PROGRAM Spline_MCHF

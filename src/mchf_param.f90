!====================================================================
    MODULE mchf_param
!====================================================================
!
!   Parameters for solving the HF equations
!     Tests for convergence:
!     . scf_tol :  relative change in the total energy
!     . orb_tol :  maximum change in an orbital relative to its
!                  maximum value
!     . cfg_tol :  maximum change in mixing coefficient
!     . end_tol :  values in the tail all less than this value in
!                  magnitude are set to zero.
!     Other parameters:
!     . acc_max :  Levels of accuracy (1 for redo)
!     . etl_sum :  Weighted sum of total energy (for redo)
!
!--------------------------------------------------------------------
       IMPLICIT NONE

       CHARACTER(LEN=4) :: varied1 = 'all', varied2 = 'all'
       CHARACTER(LEN=80):: vstring = ' '
       REAL(KIND=8)     :: scf_tol = 1.d-11, &
                           orb_tol = 1.d-05, &
                           cfg_tol = 1.d-05, &
                           end_tol = 1.d-05, &
                           etl_sum = 0.d0
       INTEGER          :: acc, acc_max = 2
       !INTEGER          :: acc, acc_max=2, max_it1, max_it2

    END MODULE mchf_param

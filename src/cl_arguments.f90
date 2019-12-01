!====================================================================
    MODULE CL_arguments
!====================================================================
!
!   Parameters defining the problem   (not yet implemented)
!     Tests for convergence:
!     . CL_atom      : name of the case.
!         input files: cfg.inp and optionally bsw.inp
!     . CL_Z         : atomic number
!        output files: name.ZZ.bsw, name.ZZ.Term.l, name.ZZ.log, name.dat
!
!   Parameters defining the iteration  (first level)
!     . CL_scf_tol :  maximum relative change in the energy defining
!                     convergence
!     . CL_orb_tol :  maximum change in an orbital relative to its
!                     maximum value
!     . CL_cfg_tol :  maximum change in mixing coefficient
!     . CL_end_tol :  values in the tail all less than this value in
!                     magnitude are set to zero.
!     . CL_acc_max :  Levels of accuracy (1 for redo)
!     . CL_max_it1 :  maximum number of iterations for SCF phase
!     . CL_max_it2 :  maximum number of iterations for ALL_NR phase
!     . CL_varied1 :  orbitals to be varied in the SCF phase
!     . CL_varied2 :  orbitals to be varied in the ALL_NR phase
!
!    Paramerters defining the spline grid
!     . CL_h       :  grid step-size parameter
!     . CL_ns      :  number of B-spline basis functions
!     . CL_ks      :  order of the B-splines
!
!--------------------------------------------------------------------
       IMPLICIT NONE

       CHARACTER(LEN=2) :: CL_atom = '-1'

       CHARACTER(LEN=4) :: CL_varied1 = '-1', CL_varied2 = '-1'
       REAL(KIND=8)     :: CL_scf_tol, CL_orb_tol, CL_cfg_tol, CL_end_tol, &
                           CL_Z
       INTEGER          :: CL_acc_max, CL_max_it1, CL_max_it2, &
                           CL_h, CL_ns, CL_ks

    END MODULE CL_arguments

!    SUBROUTINE DVDSON(A, INDROW, INDCOL, IUPPER, NZER, TM, TP, N, LIM, DIAG, &
!        ILOW, IHIGH, ISELEC, NIV, MBLOCK, CRITE, CRITC, CRITR, ORTHO, MAXITER&
!        , WORK, IWRSZ, IWORK, IIWSZ, HIEND, NLOOPS, NMV, IERR)
!  Stored in  (initialize dby call alloc_diagonalize)
! block_param:  
!    nze, lim, iiwsz, iworksz (parameters)
!    tm, tp, diag, en, wt, iwork(iiwsz)
!    hmx -- sparse form of the matrix (nze)
!    eigvec - (nume+1)*ncfg 
!    arrays to be allocted by the routine that computes hmx 
!====================================================================
    MODULE dvdson_param 
!====================================================================
!
!   Parameters for the dvdson code
!     Tests for convergence:
!     Size of work arrrays
!     Work arrays
!--------------------------------------------------------------------
    USE block_param
    IMPLICIT NONE
   
    REAL(KIND=8) :: crite =1.d-15, &
                    critc =1.d-09, &
                    critr =1.d-09, &
                    trhold=1.d-00, &
                    ortho =1.d-00
    INTEGER      :: ilow = -1, ihigh=-1
    INTEGER      :: iupper, iselec, niv, mblock, maxiter, nloops, nmv
    LOGICAL      :: hiend = .false.
   
    REAL(KIND=8), DIMENSION(:), allocatable, target:: hmx, eigvec
    END MODULE dvdson_param 

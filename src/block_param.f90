!====================================================================
    MODULE block_param
!====================================================================
!   Information about the various blocks that depend on paramters
!   nblock -- the number of blocks (each a given term and parity)
!   meig  -- the maximum number of eigenvalues in a block
!   nume  -- maximum index of eigenstates within a block (<= meig)
!   niv_bl -- number of eigenstates within a block  (<= nume)
!   leigen-- logical variable to define whethe a given eigenvalue is
!             to be considered
!   eigst_weight -- weight of an eigen state
!   en    --  list of total energies for states
!   isom_shft    -- isotope shift parameter for an eigen state
!
!--------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER :: meig, nblock, ncfg, nze

    INTEGER, DIMENSION(:), allocatable :: nze_bl, ncfg_bl, cf_tot, &
                           nze_max, nume, niv_bl
    REAL(KIND=8), DIMENSION(:,:), allocatable :: eigst_weight, isom_shift
    LOGICAL,      DIMENSION(:,:), allocatable :: leigen
    CHARACTER(LEN=3), DIMENSION(:), allocatable :: term_bl

    REAL(KIND=8), DIMENSION(:), allocatable, target :: wt, en
    INTEGER, DIMENSION(:), pointer      :: jptr

    CONTAINS
!====================================================================
       SUBROUTINE alloc_block
!====================================================================
!      Allocate arrays that depend on the number of blocks
!--------------------------------------------------------------------
       IMPLICIT NONE

       if (allocated(nze_bl)) deallocate(nze_bl, ncfg_bl, cf_tot,  &
                                         nze_max, nume, term_bl)
       allocate(nze_bl(nblock), ncfg_bl(nblock), cf_tot(nblock),   &
                nze_max(nblock), nume(nblock), niv_bl(nblock),     &
                term_bl(nblock))

       END SUBROUTINE alloc_block

!====================================================================
       SUBROUTINE alloc_eigenstates
!====================================================================
!      Allocate arrays that depend on meig and the number of blocks
!--------------------------------------------------------------------
       IMPLICIT NONE

       allocate(eigst_weight(meig,nblock), isom_shift(meig, nblock), &
                leigen(meig,nblock))
       eigst_weight = 0.d0; isom_shift=0.d0

       END SUBROUTINE alloc_eigenstates

!====================================================================
       SUBROUTINE alloc_wt
!====================================================================
!      Allocate arrays for the matrix diagonalization using
!      Dvdson
!--------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER :: i, j, k

      i = sum(nume*ncfg_bl)
      j = sum(ncfg_bl)
      k = sum(niv_bl)
      if (i.eq.1) i=2*nblock
      allocate(wt(1:i), jptr(1:j), en(1:k))

       END SUBROUTINE alloc_wt

    END MODULE block_param


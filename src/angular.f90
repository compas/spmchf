!====================================================================
    MODULE angular_data
!====================================================================
!   Arrays associated with angular data
!--------------------------------------------------------------------
       USE memory_use
       USE block_param

       IMPLICIT NONE

       INTEGER :: idim, ncodim  ! parameters read from cfg.h
       INTEGER, PARAMETER    :: maxclst = 250000000  ! maximum size
       INTEGER, DIMENSION(6)                      :: intptr
       INTEGER, DIMENSION(:), ALLOCATABLE   :: kval, jh, nijptr, inptr
       INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: ih, ico
       LOGICAL, DIMENSION(:), ALLOCATABLE         :: lused
       REAL(KIND=8), DIMENSION(:), ALLOCATABLE    :: value, coeff

    CONTAINS
!----------------------------------------------------------------------
!     A L C S T S
!----------------------------------------------------------------------
!    This routine allocates arrays associated with states.  For the
!    yint.lst arrays, memory needs to be allcoated for all blocks,
!    but c.lst arrays are read in groups of size lsdim.

!    Assumes that all arrays can be in memory.

       SUBROUTINE ALCSTS
          USE mchf_inout
          IMPLICIT NONE

          INTEGER :: nze_tot, ncfg_tot, nnn, nze_max_bl, nze_max_col, ierr

          nze_tot = sum(nze_bl(1:nblock))
          ncfg_tot = sum(ncfg_bl(1:nblock))
          nnn = sum(cf_tot(1:nblock)); 
          !  nze_max_bl is maximum nze of a block
          nze_max_bl = maxval(nze_bl(1:nblock))
          !  nze_max_col is maximum nze of a cloumn of the matrix
          nze_max_col = maxval(nze_max(1:nblock))

          ! allocate for all in memory
          allocate (kval(idim), value(idim), lused(idim))
          allocate (ico(nze_tot), ih(nze_tot), STAT=ierr)
          if (ierr /= 0) then
             write (UNIT=err, FMT='(A/A,I12)') 'ALCSTS: Not enough memory', &
                'NZE_TOT =', nze_tot
             STOP
          end if
          allocate (coeff(nnn), inptr(nnn), STAT=ierr)
          if (ierr /= 0) then
             write (UNIT=err, FMT='(A)') 'ALCSTS: Not enough memory', &
                'CF_TOT =', nnn
             STOP
          end if
          clst_memory = .true.

       END SUBROUTINE ALCSTS
    End MODULE angular_data

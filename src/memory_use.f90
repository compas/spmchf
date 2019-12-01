!====================================================================
    MODULE memory_use
!====================================================================
!
!   Logical parameters controlling memory/disk usage:
!     Default value:: .true.
!
!--------------------------------------------------------------------
    IMPLICIT NONE

    ! .. all blocks in memory
    Logical :: hmx_memory = .true. , &  ! Matrix in memory
               ico_memory = .true. , &  ! ico in memory
                ih_memory = .true. , &  ! ih.nb.lst  in memory
              clst_memory = .true.       ! c.lst in memory
    ! .. largest  block in memory
    Logical :: diag_hmx_memory = .true., &
               diag_ih_memory  = .true., &
               diag_ico_memory = .true.

    INTEGER :: idisk

    END MODULE memory_use

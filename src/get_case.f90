!======================================================================
      SUBROUTINE get_case
!======================================================================
!
!     This routine obtains information about the problem to be solved
!     and how the spline methods are to be applied

!----------------------------------------------------------------------
      Use av_energy
      Use mchf_atomic_state
      Use orbitals
      
      LOGICAL :: done
      INTEGER :: iselec

      ! .. Header information
      WRITE(6,9)
9     FORMAT(///////22X,'===========================================',&
                   /22X,'       S P L I N E  M C H F  : 2010',&
                   /22X,'==========================================='/)
 
      CALL init

      ! .. get the atomic state problem
      CALL get_atom

      CALL intgrl

      ! .. define the spline step and basis size
      CALL get_spline_param(z)

      ! .. now we know both the nwf and ns
      CALL allocate_orbital_arrays
    END SUBROUTINE get_case

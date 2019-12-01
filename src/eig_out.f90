!=======================================================================
      subroutine eig_out(iblock,unit_n)
!=======================================================================

!     This routine writes first the eigenvectors, and then
!   the eigenvalues to Term.l 
!-----------------------------------------------------------------------
      USE mchf_inout
      USE block_param
      USE mchf_atomic_state
      IMPLICIT none
      INTEGER, INTENT(IN) :: iblock, unit_n
      
      integer q(8)
      CHARACTER (LEN=3)         :: elc(8), couple(15)
      character (len = 64)	:: config_tmp, string
      character (len = 66), allocatable, dimension(:) :: config0
      character (len = 66)      :: config(MEIG)
      integer, dimension(1)	:: mx
      integer, dimension (meig) :: max_comp
      integer		 	:: max_read, neig, i,j, ics, ice, &
                                   ic1, nocc, jj, nnel, lstring

      ncfg = ncfg_bl(iblock) 
      neig = niv_bl(iblock)

      !write(60,'(I5,3X,A3)') neig, term_bl(iblock)

!     ... read 2 blank lines from file 'cfg.inp'
      Read (iuc,'(A64)') string
      Read (iuc,'(A64)') string

!<<< find the configurations with max comp of the eigen vector
      ics = 1
      ice = ncfg 
!     ... find the locations of extrema for all eigenvalues
      do ic1 = 1, nume(iblock)
	mx = MaxLoc(abs(wt(ics:ice)))
        max_comp(ic1) = mx(1)
	ics = ics + ncfg
        ice = ice + ncfg 
      end do

!     ... find the maximum of max_comp to read max_read lines  
      max_read = MaxVal(max_comp)
!
!     allocate enough space to read the cfg
      allocate(config0(max_read+1));

       !print*, max_read , '::::::::'

!     ... read all lines.le.max_read of file 'cfg.inp'
      do i = 1, max_read 
        READ(iuc,'(8(1X,A3,1X,I2,1X))') (ELC(j),Q(j),j=1,8)
        READ(iuc,'(15(1X,A3))') (COUPLE(J),J=1,15)
        
!     ... find NOCC
        nocc = 0
        do while(ELC(nocc+1).ne. '   ')
           NOCC = NOCC + 1
        end do

!     ... pack config 
        call pack(nocc,elc,q,couple,config_tmp)
!        lstring = Len_Trim(config_tmp)
        config0(i) = trim(config_tmp)
      end do
!>>>
!     ... order each configurations to correspond to max comp of the eig vector
      do i = 1, niv_bl(iblock)
        config(i) = config0(max_comp(i))
      end do
      
      deallocate(config0);
      jj = 0
      nnel = 0
      write(unit_n,'(2X,A6,A,F5.1,A,I3,A,I6)' ) ATOM,'  Z = ',Z,&
             '  NEL = ', nnel, '   NCFG = ',NCFG
      WRITE (unit_n, '(//A8,I4,2X,A8,I4)' ) '  2*J = ',JJ, &
           'NUMBER =', neig;
!     ... write $TERM.l file for each block
      do i = 1, nume(iblock);
         if (leigen(i,iblock)) then;
           config_tmp = config(i) 
           lstring = Len_Trim(config_tmp)
           write (unit_n,'(4x,A8,f16.9)') 'Ssms = ',&
                 isom_shift(i,iblock);
           WRITE (unit_n,'(i6,f16.9,2x,A)') &
                  max_comp(i), en(i), config_tmp(1:lstring) 
           WRITE (unit_n,'(7F11.8)') (wt((i-1)*ncfg+j), j=1,ncfg)
!          WRITE(iscw,'(/1X,F19.12,2X,A/(1X,7F11.8))' ) EIGVAL
         end if
      end do
!#####
      do while (string(1:1).ne.'*')
         Read (iuc,'(A64)') string
      end do

      return
      end

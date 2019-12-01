!----------------------------------------------------------------------
!     I N T G R L
!----------------------------------------------------------------------
!     This routine reads the yint.lst data storing the concatenated
!     arrays for the different blocks

      SUBROUTINE INTGRL
      USE  mchf_inout
      USE mchf_atomic_state, ONLY: nwf, nclosd
      USE memory_use
      USE  block_param
      USE  angular_data
      IMPLICIT  NONE

      INTEGER :: maxorb, ib1, ib2, ncl, mwf, nbb, lsd, &
                 count, nij, n, iblock, ierr, i, ncol, &
                 lasti, int, icount, itype, nint, ibe, &
                 n_left, n_read, j
      CHARACTER(LEN=72) :: string
      CHARACTER(LEN=2)  :: ih_file

      rewind(iuy)
      read(iuy) ncl, mwf, nbb, lsd
      maxorb = nwf-nclosd
      if (ncl.ne.nclosd.or.mwf.ne.maxorb.or.nbb.ne.nblock) then
         write(err,*) 'Yint.lst data not consistent with cfg.h data'
         write(err,*) 'Yint.lst :', ncl, mwf, nbb
         write(err,*) 'cfg.h    :', nclosd, nwf, nblock
         stop
      end if
      read(iuy) string

      do while (maxorb > 0) 
        read(iuy) string
        maxorb = maxorb - 24
      end do

      ! Allocate memory for c.lst and pointer data
      call alcsts
      jptr = 0; kval = 0
       ib1 = 0; ib2 = 0

      !  read for each block
      do iblock = 1, nblock
        if (hmx_memory) then
          write(ih_file,'(I2.2)') iblock 
          open(iih,file='ih.'//ih_file//'.lst',status='old', &
               form='unformatted');
          nij = 0;
          n = lsd;
          do while(nij < nze_bl(iblock));
            read(iih,iostat=ierr) n, (ih(i+nij),i=ib1+1,ib1+n);
            if (ico_memory) then
              read(iic,iostat=ierr) n, (ico(i+nij),i=ib1+1,ib1+n);
            end if;
            nij = nij + n;
          end do;
        end if;
      read(iuy,iostat=ierr) ncol,(jptr(i),i=ib2+1,(ib2+ncol))
      ib1 = ib1 + nij
      ib2 = ib2 + ncol
      close (iih);
      end do
      
        nij = 0
 
! ***** READ  THE LIST OF INTEGRALS
 
        lasti = 0
        
!          ...F, G, R, or L integrals....
      DO INT = 1,4
        icount = 1
        read(iuy) itype, nint

!       .. ipackn of nonh is called kval here 
!       .. integral pack number and logical variable indicating usage
        read(iuy) (kval(i),i=lasti+1,nint), &
                  (lused(i),i=lasti+1,nint)
        lasti = nint
        icount = icount + 1
        intptr(int) = lasti
      end do

!      .. adjust for the fact that overlaps have been removed

      intptr(6) = intptr(4)
      intptr(5) = intptr(3)
      intptr(4) = intptr(2)
      intptr(3) = intptr(2)
       
      if (clst_memory) then
        ibe = 0;
        do iblock = 1, nblock
          print *, 'iblock, cf_tot, nocdim', iblock, cf_tot(iblock), ncodim
          n_left = cf_tot(iblock)
          if (cf_tot(iblock) > ncodim) then
            n_read = ncodim; 
          else
            n_read = cf_tot(iblock);
          end if
          print *, 'n_left, n_read', n_left, n_read

          do while (n_left > 0)
            if (n_left<ncodim) n_read = n_left;
            read(icl) n, (coeff(j),j = ibe+1,ibe+n_read), &
                          (inptr(j), j = ibe+1,ibe+n_read)
            ibe = ibe + n_read;
            n_left = n_left - n_read
          end do
       !print *, 'Intgrl: inptr, coeff'
       !print '(8i10)',(inptr(i),i=1,n_read)
       !print '(4F12.8)', (coeff(i),i=1,n_read) 
       end do
      end if
      end SUBROUTINE INTGRL


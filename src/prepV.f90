!     ==================================================================
      SUBROUTINE prepV(last)
!     ==================================================================

         USE angular_data
         USE mchf_atomic_state
         USE mchf_inout
         USE memory_use
         USE mchf_param
         USE block_param
         USE dvdson_param
         IMPLICIT NONE
         Logical, INTENT(IN) :: last

         INTEGER :: ivs, ijp, nz, ibe, ncoef, iblock, max_col, nijcurr, i, &
                    ies, unit_n
         CHARACTER(LEN=5) :: name_dotL

         ivs = 1; ijp = 1; nz = 1
         ibe = 0; ncoef = 0; 
!<<<<<<<<<<<<<<<<<<<updatc>>>>>>>>>>>>>>>>>>>>>>>>>>>.
         do iblock = 1, nblock
            ncfg = ncfg_bl(iblock)
            if (clst_memory) then
               CALL UPDATC_memory_all(iblock, ivs, ncfg, nume(iblock), ijp, &
                                      nz, max_col, nijcurr, last)
            end if
            !ivs = ivs + ncfg*nume(iblock)
            ivs = ivs + ncfg*niv_bl(iblock)
            ijp = ijp + ncfg_bl(iblock)
            nz = nz + nze_bl(iblock)
         end do
!<<<<<<<<<<<<<<<<<<end updatc>>>>>>>>>>>>>>>>>>>>>>>>>>>>
         DO I = 1, nwf
            IF (qSUM(I) .EQ. 0.0) write (err, *) 'DIAG: SUM=', qSUM(I), ' &
               for ', EL(I)
         end do

         if (last) then
            ivs = 1
            ies = 1
            do iblock = 1, nblock
               unit_n = 40 + iblock
               name_dotL = term_bl(iblock)//'.l'; 
               open (unit_n, file=trim(name_dotL), status='unknown', form='formatted')
               call eig_out(iblock, unit_n)
               ivs = ivs + ncfg_bl(iblock)*niv_bl(iblock)
               ies = ies + niv_bl(iblock)
            end do
         endif

         !print *, 'Values:'
         !print '(3(I6,F18.15))', (i,value(i),i=1,intptr(6))

      END SUBROUTINE prepV

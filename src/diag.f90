
!     ==================================================================
      SUBROUTINE DIAG(sum_energy)
!     ==================================================================
!     neiv : maximum index of an eigenstate
!     sum_energy  -- weighted sum of the total energy
!     ------------------------------------------------------------------

         USE angular_data
         USE mchf_atomic_state
         USE mchf_inout
         USE memory_use
         USE mchf_param
         USE block_param
         USE dvdson_param

         IMPLICIT NONE
         REAL(KIND=8), INTENT(OUT) :: sum_energy
         INTEGER :: I, J, iblock, ibe, ibw, ibv, nz, nj, ijp, ivs, ies, &
                    ncoef, neiv, nijcurr, max_col, jjh, icc, n_count_tmp
         INTEGER, DIMENSION(1) ::  jloc
         REAL(KIND=8) :: etl, shift
!------------------------------------------------------------------------
!     Reminder:  eigvec is for the local block (all eigenvectors 1..nume)
!                wt is for the set of blocks but contains
!                   only the selected eigenvectors and eigenvalues
!                niv counts the selected eigenvectors globally
!------------------------------------------------------------------------

         idisk = 0; 
         etl = 0; sum_energy = 0; nz = 0; nj = 0
         ijp = 1; ivs = 1; ies = 1; 
         ibe = 0; ncoef = 0; 
         ibw = 0; ibe = 0
         do iblock = 1, nblock
!       .. set parameters to current block
            ncfg = ncfg_bl(iblock)
            nze = nze_bl(iblock)
            neiv = nume(iblock)
            allocate (hmx(nze), eigvec((ncfg + 1)*neiv))
            nijcurr = 1
            max_col = 1; 
            jjh = 1
            diag_ih_memory = .true.
            diag_ico_memory = .true.
            diag_hmx_memory = .true.

            if (clst_memory) call compute_hmx

            !print *, 'hmx:'
            !print '(3F20.16)', (hmx(i),i=1,jptr(ncfg))
            !print *, 'ih:'
            !print '(4I6)', (ih(i+nz),i=1,jptr(ncfg))
            !print *, 'jptr:', (jptr(i+nj),i=1,ncfg)
            !print *, 'iupper,nze,ncfg,lim:',iupper,nze,ncfg,lim
            !print *, 'diag :, ncfg'
            !print '(3F20.16)',(hii(i),i=1,ncfg)
            !print *, 'ilow,ihigh,iselec(1:ie)', &
            !       ilow,ihigh,iselec(1:ie)

            if (ncfg == 1) then
               wt(ibw + 1) = 1.d0
               en(ibe + 1) = hmx(1) + ec
               eigvec(1) = 1.d0
               eigvec(2) = hmx(1) + ec
               sum_energy = sum_energy + eigvec(2)*eigst_weight(1, iblock)
               !print *, 'wt, en, sum', wt(ibw+1), en(ibe+1), sum_energy, &
               !          eigst_weight(i,iblock)
               ibw = ibw + ncfg
               ibe = ibe + 1
            else

               ! full matrix
               call iniest2(2000, ncfg, neiv)
!        .. save eigenvectors and eigenvalues
!        wt -- global (selected), eigvec -- local block(all)
               ibv = 0
               do icc = 1, nume(iblock)
                  if (leigen(icc, iblock)) then
                     jloc = maxloc(abs(eigvec(ibv + 1:ibv + ncfg)))
                     if (eigvec(jloc(1)) .lt. 0) then
                        eigvec(ibv + 1:ibv + ncfg) = -eigvec(ibv + 1:ibv + ncfg)
                     end if
                     wt(ibw + 1:ibw + ncfg) = eigvec(ibv + 1:ibv + ncfg)
                     etl = eigvec(ncfg*nume(iblock) + icc) + ec
                     ibe = ibe + 1
                     en(ibe) = etl
                     sum_energy = sum_energy + etl*eigst_weight(icc, iblock)
                     ibw = ibw + ncfg
                  end if
                  ibv = ibv + ncfg
               end do
            end if
         end do
         deallocate (hmx, eigvec)

      CONTAINS
!======================================================================
         subroutine compute_hmx
!======================================================================
            IMPLICIT NONE
            INTEGER :: n_cf, nze, ii

            n_cf = cf_tot(iblock)
            nze = nze_bl(iblock)
            hmx(1:nze) = 0.0
            max_col = 1; 
            do ii = 1, n_cf; 
!       .. test for next non-zero matrix element
               n_count_tmp = ncoef + ii
               if (ii .gt. ico(nijcurr + nz)) nijcurr = nijcurr + 1
               hmx(nijcurr) = hmx(nijcurr) + &
                              coeff(n_count_tmp)*value(inptr(n_count_tmp))
               if (nijcurr .gt. jptr(jjh)) then
                  jjh = jjh + 1; 
                  max_col = max_col + 1
               end if
               ! print '(3I5,3F20.15,I4)',  &
               ! ii,nijcurr,inptr(ncoef+ii),coeff(ncoef+ii),&
               !             value(inptr(ncoef+ii)), &
               ! hmx(nijcurr)
            end do

            ncoef = ncoef + n_cf; 
         END SUBROUTINE compute_hmx
      END SUBROUTINE DIAG

!=======================================================================
      SUBROUTINE UPDATC_memory_all(iblock,ivs,ncfg1,maxev,ijp,nz,   &
                 max_col,nijcurr,last)
!=======================================================================
!
!     Evaluate all coeffients of integrals.  This involves a
!     sequential pass through all the c.lst data. Use same array
!     as "value" now called "coef"
!
! NOTE: jptrp, icop and associate aruments are not used.
!-----------------------------------------------------------------------
      USE angular_data, coef=>value
      USE mchf_atomic_state,  ONLY: qsum, clsd, l
      USE block_param
      USE dvdson_param

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: iblock, ivs, ncfg1, maxev, ijp, nz
      INTEGER, INTENT(OUT):: max_col, nijcurr
      LOGICAL, INTENT(IN) :: last

      INTEGER,     DIMENSION(:), pointer    :: jptrp, ihh, icop
      REAL(KIND=8),DIMENSION(:), pointer    :: wtt
      REAL(KIND=8),DIMENSION(:),allocatable :: coef_tmp
      INTEGER :: ncoef, iih, nij_tmp, jjh, n_count_tmp, i, im1, &
                 j, ibv, itmp_s, itmpu, im2, istart, ibegin, iend, &
                 kv, iel1, iel2, iel3, iel4
      REAL(KIND=8) :: wcoef, w,  t, w0, t0
      save ncoef

      if ((iblock.eq.1)) coef(1:intptr(6)) = 0.d0
      wtt => wt(ivs:); jptrp => jptr(ijp:)
      ihh => ih(nz:);  icop => ico(nz:)

      if (last) then
         allocate(coef_tmp(niv_bl(iblock)*idim))
         coef_tmp = 0.d0
      end if

      nijcurr = 1; iih = 1 ; nij_tmp = 1; jjh = 1  ; max_col = 1;
      n_count_tmp = 0;
      if (iblock == 1) ncoef = 0;
      !        .. add contribution from this record to coef
      !        .. make use of the fact that all coefficients are read in
      !           increasing order of the nze elements of the matrix
      do i = 1, cf_tot(iblock);
            !           ..test for next non-zero element
            n_count_tmp = ncoef+i
            if (i.gt.ico(nijcurr)) then
                nijcurr = nijcurr + 1
                !               .. have we also changed column?
                if (nijcurr.gt.jptr(max_col)) then
                   jjh = jjh + 1;
                   max_col = max_col+1
                end if
             end if
             iih = ihh(nijcurr)
             im1 = 0; ibv=0
             do j = 1,maxev
               if (leigen(j,iblock)) then
                 wcoef = eigst_weight(j,iblock)
                 W = wcoef*wtt(ibv+iih)*wtt(ibv+jjh)
                 T = W*coeff(n_count_tmp)
                 IF (IIH .NE. JJH) T = T+T
                 coef(inptr(n_count_tmp)) = coef(inptr(n_count_tmp)) + T
                 if (last) then
                   W0 = wtt(ibv+iih)*wtt(ibv+jjh);
                   T0 = W0*coeff(n_count_tmp);
                   IF (IIH .NE. JJH) T0 = T0 + T0;
                   itmp_s = (inptr(n_count_tmp))+idim*im1;
                   coef_tmp(itmp_s) = coef_tmp(itmp_s) + T0;
                   im1 = im1 + 1;
                 end if
                 ibv = ibv + ncfg1
               end if
             end do
 	 !    if (inptr(i) .gt. intptr(5)) then
         ! print '(2I4,3F15.10,I6,F15.10,i6)', &
         !iih,jjh,wtt(iih+ibv),wtt(jjh+ibv),coeff(n_count_tmp), &
         !inptr(n_count_tmp), coef(inptr(n_count_tmp)),nijcurr
  	! end if
      end do

      ncoef = ncoef + cf_tot(iblock);

      if (last) then
         im2 = 0;
         do im1 = 1,maxev
           if (leigen(im1,iblock)) then
             istart = 1+idim*im2;
             call prprty(im1,iblock, isom_shift(im1,iblock));
             im2 = im2 + 1;
           end if
         end do
         deallocate(coef_tmp);
      end if
!
!  *****  DEFINE qSUM(I)
!
         IBEGIN = INTPTR(5)+1
         IEND = INTPTR(6)

         DO I = IBEGIN,IEND
           call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
           IF (IEL1 .EQ. IEL2) then
              qSUM(IEL1) = -2*COEF(I)
              if ( abs(qsum(iel1) - 2*(2*l(iel1)+1)) <= 1.0D-5) &
                   clsd(iel1) = .true.
           end if
         end do

        !print *, 'Node:',' intptr(6)',intptr(6), 'iblock',iblock
        !print '(4f20.16)', (coef(i),i=1,intptr(6))
        print '(A,(4F20.16))', 'QSUM:', qsum
        print *, 'CLSD:', clsd


   CONTAINS

!------------------------------------------------------------------------
      SUBROUTINE PRPRTY(ieigval,iblock,SS)
!------------------------------------------------------------------------
!
      USE av_energy
      USE mchf_inout
      USE mchf_atomic_state
      USE block_param
      USE angular_data

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ieigval, iblock
      REAL(KIND=8), INTENT(OUT) :: ss
      REAL(KIND=8), EXTERNAL :: quadr, grad

      INTEGER :: i, ibegin, iend, j, k, kv, iel1, iel2, iel3, iel4
      REAL(KIND=8) :: rmean1, rmean2, R_R, sumi, sumj, rr1, gr1, c
!

      WRITE(log,'(//15X,A/)')  'Some WaveFunction Properties'
      rmean1 = 0.d0
      rmean2 = 0.d0
      SS = 0.d0
      R_R = 0.0

      Do I = 1,NCLOSD
	sumi = 4*L(i)+2
	rmean1 = rmean1 + sumi*quadr(i,i,1)
	rmean2 = rmean2 + sumi*quadr(i,i,2)
      END DO
!
!     a) From Common core interactions: 2-body  operators
      Do I = 2,NCLOSD
        sumi = 4*L(i)+2
        Do J =  1,I-1
          sumj = 4*L(j)+2
          DO K=IABS(L(I)-L(J)),L(I)+L(J),2
            IF (K.EQ.1) THEN
              C= -CB(L(I),L(J),K)
!             .. this minus is from exchange term
              SS = SS - C*SUMI*SUMJ*GRAD(I,J)**2
	      R_R = R_R + c*SUMI*SUMJ*quadr(i,j,1)**2
!	      print *, i,j, -c*sumi*sumj,-grad(i,j)**2, ss
            ENDIF
          END DO
        END DO
      END DO

!
!     b) From outer electrons
!       i) From G-integrals
      IBEGIN = INTPTR(1)+1
      IEND = INTPTR(2)
      Do I = IBEGIN, IEND
	call unpacki(2,i,kv,iel1,iel2,iel3,iel4)
	IF (kv .eq. 1) then
!         .. this minus is from exchange form
	  SS = SS - coef_tmp(i)*grad(iel1,iel2)**2
	  R_R = R_R +coef_tmp(i)*quadr(iel1,iel2,1)**2
!	  print *, iel1, iel2, coef_tmp(i),-grad(iel1,iel2)**2, ss
	END IF
      END DO

!      ii) From R-integrals
      IBEGIN = INTPTR(4)+1
      IEND = INTPTR(5)
      Do I = IBEGIN, IEND
	IF (coef_tmp(i) .ne. 0.d0 ) THEN
	  call unpacki(5,i,kv,iel1,iel2,iel3,iel4)
	  IF (kv .eq. 1) then
	    SS = SS + coef_tmp(i)*grad(iel1,iel3)*grad(iel2,iel4)
	    R_R = R_R +coef_tmp(i)*quadr(iel1,iel3,1)*quadr(iel2,iel4,1)
!	    print *, iel1, iel3, iel2, iel4, coef_tmp(i),
!    : 	    grad(iel1,iel3)*grad(iel2,iel4), ss
	  END IF
	END IF
      END DO

!     c) Contributions from outer electrons and common core

      IBEGIN = intptr(5)+1
      IEND = intptr(6)
      DO i = ibegin, iend
        call unpacki(6,i,kv,iel1,iel2,iel3,iel4)
        SUMI= -2*coef_tmp(i)
	rmean1 = rmean1 + sumi*quadr(iel1, iel2, 1)
	rmean2 = rmean2 + sumi*quadr(iel1, iel2, 2)
!	print *, i,iel1,iel2,sumi
        DO  J=1,NCLOSD
          SUMJ=4*L(J)+2
          RR1 = 0.d0
          GR1 = 0.d0
          DO K=IABS(L(iel1)-L(J)),L(iel1)+L(J),2
              IF (K.EQ.1) THEN
                C= -CB(L(IEL1),L(J),K)
                IF (IEL1 .EQ. IEL2 ) THEN
!                    .. we have a diagonal case
                     RR1=C*QUADR(IEL1,J,1)*QUADR(J,IEL2,1)
!                    .. this is from exchange form
                     GR1=-C*GRAD(IEL1,J)**2
!                    Print *, iel1,j,j,iel2,-c*sumi*sumj,
!    :                     GRAD(IEL1,J)**2
                  ELSE
!                     .. we have an off-diagonal case
                     GR1=  C*GRAD(IEL1,J)*GRAD(J,IEL2)
		     RR1=  C*Quadr(IEL1,J,1)*Quadr(J,IEL2,1)
!		     Print *, iel1,j,j,iel2,c*sumi*sumj,
!    :		     GRAD(IEL1,J)*GRAD(J,IEL2)
                  ENDIF
               ENDIF
         END DO
            SS = SS + SUMI*SUMJ*GR1
            R_R = R_R + SUMI*SUMJ*RR1
        END DO
      END DO

!     Change sign of SS to adhere to our definition
      SS = - SS
      write (log,'(A8,1X,A3)')  '    Term ', term_bl(iblock)
      WRITE(log,'(A,F15.8)') '    Mean radius             =', rmean1
      WRITE(log,'(A,F15.8)') '    Mean square radius      =', rmean2
      WRITE(log,'(A,F15.8)') '    Mean R.R parameter      =', rmean2 + 2*r_r
      WRITE(log,'(A,F15.8)') '    Isotope Shift parameter =', SS
      END SUBROUTINE prprty
      END  SUBROUTINE updatc_memory_all

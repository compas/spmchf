!=====================================================================
      subroutine update
!=====================================================================
!     Compute the radial integrals
!------------------------------------------------------------------
         !USE spline_param, ONLY: h,hmax,rmax,ns,ks,ml
         !USE spline_galerkin, ONLY: sb,bb
         USE mchf_atomic_state
         USE mchf_inout
         USE orbitals
         USE mchf_param
         USE angular_data

         IMPLICIT NONE
         INTEGER :: len, ibegin, iend, i, k1, k2, kv, iel1, iel2, iel3, iel4
         REAL(KIND=8), EXTERNAL :: hlc, rk, fk, gk, quadr
         INTEGER, DIMENSION(4) :: ir
         INTEGER :: ii, jj
         CHARACTER(LEN=4) :: case

         LOGICAL change
!
         len = intptr(6)
         value(1:len) = 0.d0
         IBEGIN = 1
         IEND = INTPTR(3)
         ! Fk, Gk
         DO I = IBEGIN, IEND
            if (lused(i)) then
               call unpacki(1, i, kv, iel1, iel2, iel3, iel4)
               IF (I .LE. INTPTR(1)) THEN
                  VALUE(I) = FK(IEL1, IEL2, KV)
               ELSE IF (I .LE. INTPTR(2)) THEN
                  VALUE(I) = GK(IEL1, IEL2, KV)
               ELSE
                  VALUE(I) = QUADR(IEL1, IEL2, 0)**KV
               END IF
            end if
         End do
!
         ! Quadr*Quadr
         IBEGIN = IEND + 1
         IEND = INTPTR(4)
         DO I = IBEGIN, IEND
            if (lused(i)) then
               call unpacki(4, i, kv, iel1, iel2, iel3, iel4)
               K1 = KVAL(I)/16
               K2 = KVAL(I) - 16*K1
               VALUE(I) = QUADR(IEL1, IEL2, 0)**K1 &
                          *QUADR(IEL3, IEL4, 0)**K2
            end if
         End do

         ! Rk
         IBEGIN = IEND + 1
         IEND = INTPTR(5)
         DO I = IBEGIN, IEND
            if (lused(i)) then
               call unpacki(5, i, kv, iel1, iel2, iel3, iel4)
               VALUE(I) = RK(IEL1, IEL2, IEL3, IEL4, KV)
            end if
         End do
!
         ! HLc
         IBEGIN = IEND + 1
         IEND = INTPTR(6)
         DO I = IBEGIN, IEND
            if (lused(i)) then
               call unpacki(6, i, kv, iel1, iel2, iel3, iel4)
               ! IF (vary(IEL1) .OR. vary(IEL2))     &
               VALUE(I) = HLC(IEL1, IEL2)
            end if
         End do
         ! print *, 'Values:'
         ! print '(3(I4,F20.15))', (i,value(i),i=1,intptr(6))

!      ... Test if any of the core functions have changed
!
         CHANGE = .FALSE.
         DO I = 1, NCLOSD
            CHANGE = CHANGE .OR. vary(I)
         End do
         IF (CHANGE .OR. (EC == 0.D0 .and. nclosd > 0)) CALL ECORE
!

      END subroutine update


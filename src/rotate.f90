!======================================================================
  SUBROUTINE rotate(i, j)
!======================================================================
!   Determine the rotation parameter for a  a stationary solution for
!   the (j,i) orbitals, i <= j
!   Rules derived from transformation of densities
!
!  Let P*i = (Pi - ePj)/(1+e^2)2
!      P*j = (Pj + ePi)/(1+e^2)2
!
!   We compute three values:
!   F(1) : sum of contributions to energy from terms that depend on i,j
!          with e = 0         (F0)
!   F(2) : sum with e=  +eps  (F+)
!   F(3) : sum with e = -eps  (F-)
!   where eps is a numerically selected value
!   From F(e) = F0  +e F'(0) + e^2 F"(0)/2  we get
!   F"(0) = [F+ - 2F0 - F-]/e^2
!   F'(0) = [F+ - F-]/(2e)
!----------------------------------------------------------------------
!
     Use spline_param
     USE spline_galerkin, ONLY: sb
     Use mchf_atomic_state
     USE orbitals, ONLY: p
     USE mchf_param
     USE av_energy
     USE angular_data, coef => value

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i, j
     INTEGER :: m, k, ip, jp, ibegin, iend, up, kv, ms
     INTEGER, DIMENSION(4) :: ir
     REAL(KIND=8), EXTERNAL :: hlc, fk, gk, rk, bhl, bvmv
     REAL(KIND=8), DIMENSION(ns) :: Pi, Pj, y
     !REAL (KIND=8), DIMENSION(ns,ks) :: d, dsum, hd
     !REAL (KIND=8), DIMENSION(ns,ns) :: dx, hfm, dxsum
     REAL(KIND=8):: eps, g, dg, ec_save, c
     LOGICAL :: found
     REAL(KIND=8), DIMENSION(3) :: F

     if (i > j) STOP 'rotate called with i < j'
     if (clsd(i) .and. clsd(j)) return

     print *, '--------------------- Rotate-------------------', i, j

     Pi = p(:, i); Pj = p(:, j)
     ec_save = ec; ec = 0.0d0

     eps = abs(rdg(i, j))
     if (eps == 0.d0) eps = 0.08
     Do ip = 1, 3
        if (ip == 2) then
           P(:, i) = (Pi - eps*Pj)/sqrt(1.+eps*eps)
           P(:, j) = (Pj + eps*Pj)/sqrt(1.+eps*eps)
        else if (ip == 3) then
           P(:, i) = (Pi + eps*Pj)/sqrt(1.+eps*eps)
           P(:, j) = (Pj - eps*Pj)/sqrt(1.+eps*eps)
        end if

        F(ip) = 0; 
        if (i <= nclosd) then

           F(ip) = -qsum(i)*bhl(i, i)/2
           do jp = 1, nclosd
              ! direct k = 0
              If (jp /= i) then
                 c = qsum(i)*qsum(jp)
              else
                 c = qsum(i)*(qsum(i) - 1)/2
              end if
              F(ip) = F(ip) + c*fk(i, jp, 0)
           end do

           ! direct k >0
           if (l(i) > 0) then
              do k = 2, 2*l(i), 2
                 c = -qsum(i)*(qsum(i) - 1)*ca(l(i), k)/2
                 F(ip) = F(ip) + c*fk(i, i, k)
              end do
           end if

           ! exchange
           do k = 0, kmax
              do jp = 1, nclosd
                 if (i == jp) cycle
                 if (k > l(i) + l(jp)) cycle
                 if (mod(abs(l(i) - l(jp)) - k, 2) /= 0) cycle
                 c = -qsum(i)*qsum(jp)*cb(l(i), l(jp), k)
                 F(ip) = F(ip) + c*gk(i, jp, k)
              end do
           end do

        end if

        IBEGIN = 1
        IEND = INTPTR(2)
        DO jp = IBEGIN, IEND
           if (lused(jp) .and. abs(coef(jp)) > orb_tol) THEN
              call unpacki(1, jp, kv, ir(1), ir(2), ir(3), ir(4))
              if (product(ir(1:2) - i) == 0 .or. product(ir(1:2) - j) == 0) then
                 iF (jp .LE. INTPTR(1)) THEN
                    F(ip) = F(ip) + coef(jp)*Fk(ir(1), ir(2), kv)
                 ELSE IF (jp .GT. INTPTR(1)) THEN
                    F(ip) = F(ip) + coef(jp)*Gk(ir(1), ir(2), kv)
                 END IF
              end if
           end if
        end do
        !
        IBEGIN = intptr(4) + 1
        IEND = INTPTR(5)
        DO jp = IBEGIN, IEND
           if (lused(jp) .and. abs(coef(jp)) > orb_tol) THEN
              call unpacki(5, jp, kv, ir(1), ir(2), ir(3), ir(4))
              if (product(ir(1:4) - i) == 0 .or. product(ir(1:4) - j) == 0) then
                 F(ip) = F(ip) + coef(jp)*Rk(ir(1), ir(2), ir(3), ir(4), KV)
              end if
           end if
        end do
        !
        IBEGIN = IEND + 1
        IEND = INTPTR(6)
        DO jp = IBEGIN, IEND
           if (lused(jp) .and. abs(coef(jp)) > orb_tol) THEN
              call unpacki(6, jp, kv, ir(1), ir(2), ir(3), ir(4))
              if (product(ir(1:2) - i) == 0 .or. product(ir(1:2) - j) == 0) then
                 F(ip) = F(ip) + coef(jp)*hlc(ir(1), ir(2))
              end if
           end if
        end do
     end do
     ! F(1) == F0, F(2) = F+, F(3) = F-

     Print *
     Print '(A,F12.8,A)', '         -h                0               h  (h=', &
        eps, ' )'
     Print *, '      --------------------------------------------'
     Print '(A,3F15.8)', 'F(e) :', F(3), F(1), F(2)
     dg = (F(2) - 2*F(1) + F(3))/(2*eps*eps)
     g = (F(2) - F(3))/(2*eps)
     eps = -g/(2*dg)
     Print '(A,3F15.8)', 'Coeff:', F(1), g, dg
     Print '(A,30X,F15.8)', 'Min  :', eps
     Print *, '      --------------------------------------------'

     rdg(i, j) = eps

     p(:, i) = (Pi - eps*Pj)/sqrt(1 + eps*eps)
     p(:, j) = (Pj + eps*Pi)/sqrt(1 + eps*eps)

     Print '(I5,2F16.9)', 1, qsum(1), -bhl(1, 1)/2.
     Print '(I5,2F16.9)', 2, qsum(2), -bhl(2, 2)/2.
     Print '(A, F16.9)', 'F0(1s,2s) ', fk(1, 2, 0)
     Print '(A, F16.9)', 'G0(1s,2s) ', gk(1, 2, 0)

     ! restore ec
     ec = ec_save
  end SUBROUTINE rotate


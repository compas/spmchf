!==================================================================
       SUBROUTINE scf_nr(i,north, jorth, hfm, bb, hx, rhs, md, v, eii)
!==================================================================
!    Improve the estimate v by the Newton_Raphson method
!    subject to orthogonality where eij and eji are computed by 
!    the lagrange routine to ensure the same values.
!------------------------------------------------------------------
        USE angular_data, coef=>value
        USE av_energy
        USE spline_galerkin, ONLY: sb
        Use spline_param
        USe mchf_atomic_state
        USE mchf_inout
        Use orbitals
        IMPLICIT NONE      
        
        INTEGER, INTENT(in) :: i, north, md
        INTEGER, DIMENSION(:), INTENT(in) :: jorth
        REAL(KIND=8), DIMENSION(ns,ns), INTENT(in) :: hfm,  bb, hx
        REAL(KIND=8), DIMENSION(ns), INTENT(inout)    ::  rhs
        REAL(KIND=8), DIMENSION(ns), INTENT(inout) :: v
        REAL(KIND=8), INTENT(out) :: eii

        REAL(KIND=8),  DIMENSION(md,md) :: aa, aaa
        REAL(KIND=8),  DIMENSION(md) :: res, rh
        !REAL(KIND=8), DIMENSION(3*md) :: work
        INTEGER(4),  DIMENSION(md) :: ipiv,  ms

        REAL(KIND=8), DIMENSION(ns) :: x, y
        REAL(KIND=8), DIMENSION(ns,ns) :: dx, dxsum
        REAL(KIND=8) :: eij, eeii,  cc 
        
        INTEGER :: mm, k, info, j, jp, ibegin, iend, kv, i1,i2,i3, i4
        
        !  .. use the same sign conventions as in scall_nr


        aa = 0.d0
       
        rhs(1) = 0.d0; rhs(ns-1:ns) = 0.d0
        x = matmul(hfm,v) + rhs;  x(1) = 0.d0; x(ns-1:ns) = 0.d0

        eeii= dot_product(x,v)
        call bxv(ks,ns,sb,v,y);  y(1) = 0.d0; y(ns-1:ns) = 0.d0
        res(1:ns) = x-eeii*y

        ! compute the residuals     
        aa(1:ns,1:ns) = (hfm -eeii*bb)
        call apply_bc(aa(1:ns,1:ns),'n')

        mm = ns+1
        ! Add normalization constraint for orbital i
        call bxv(ks,ns,sb,v,y); y(1) = 0.d0; y(ns-1:ns) = 0.d0

        aa(1:ns,mm) = -y; aa(mm,1:ns) =-y
        res(mm) = 0.d0

        !Add orthogonality constraints with orbital j
        do jp = 1,north
          mm = mm+1
          j = jorth(jp)

          eij = dot_product(p(:,j),x(:))
          ! check Lagrange multipliers
   
          eij = e(i,j)
          call bxv(ks,ns,sb,p(:,j),y)
          y(1) = 0.d0; y(ns-1:ns) = 0.d0

          aa(1:ns,mm) = -y; aa(mm,1:ns) = -y
          res(1:ns) = res(1:ns) -eij*y
          res(mm) = 0.d0
        end do

        If ( mm /= md) then
            Print *, "mm /= md"
            STOP
        end if
   
        ! Add Haa contributions
         call apply_bc(hx,'n')
         aa(1:ns,1:ns) = aa(1:ns,1:ns) + hx

       !  CALL DSYTRF('U',md, aa, md, ipiv, work, 3*md, info)
          CALL DGETRF(md, md, aa, md, ipiv, info)
         If (info /= 0) then
           WRITE(UNIT=err,'(A)') ' Error in factorization: DGETRF in scf_nr'
           STOP
         End if
         !CALL DSYTRS('U',md, 1, aa, md, ipiv, res, md, info)
         CALL DGETRS('N', md, 1, aa, md, ipiv, res, md, info)
         If (info /= 0) then
           WRITE(UNIT=err,'(A)') ' Error in solve routine : DGETRS in scf_nr'
           STOP
         End if
         
         v = v-res(1:ns)
         call normalize(v)
         eii = (eeii -res(ns+1))

     END SUBROUTINE  scf_nr

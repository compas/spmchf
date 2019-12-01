    !==================================================================
      SUBROUTINE apply_orthogonality(hfm, rhs, v)
    !==================================================================
    !   Redefine the matrix to project out the vector v so that
    !   eigenvalues of  (A' - lamda S)x = 0 are orthogonal to v.
    !------------------------------------------------------------------
        USE spline_param, ONLY: ks, ns
        USE spline_galerkin, ONLY: sb

        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(ns,ns), INTENT(INOUT) :: hfm
        REAL(KIND=8), DIMENSION(ns), INTENT(INOUT) :: rhs
        REAL(KIND=8), DIMENSION(ns), INTENT(IN) :: v
        REAL(KIND=8), DIMENSION(ns) :: z
        REAL(KIND=8), DIMENSION(ns,ns) :: cc
        INTEGER :: i, j, m

        cc = 0.d0
        DO i=1,ns
          cc(i,i) = 1.d0
        END DO

        call bxv(ks,ns,sb,v,z)
   
        DO i = 2,ns-2
          DO j = 2,ns-2
            cc(i,j) = cc(i,j) -v(i)*z(j)
          END DO
        END DO

        hfm = MATMUL(hfm,cc)
        cc = transpose(cc)

        hfm = MATMUL(cc,hfm)
        rhs = MATMUL(cc,rhs)
   

       END SUBROUTINE apply_orthogonality

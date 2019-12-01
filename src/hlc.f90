!======================================================================
    REAL(KIND=8) FUNCTION hlc(i, j)
!======================================================================
!
!   Computes hl(i,j) modified by the interactions with the closed shell
!      where H= -(1/2) [d^2/dr^2 +2Z/r -L(L+1)/r^2]
!
!---------------------------------------------------------------------
!
       USE mchf_atomic_state
       USE av_energy
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: i, j

       ! .. local variables
       INTEGER :: ip, k
       REAL(KIND=8) :: sumip, tm
       REAL(KIND=8), EXTERNAL :: rk, bhl
       REAL(KIND=8), INTRINSIC :: ABS

       hlc = bhl(i, j)
       !print *, 'hlc:', i,j
       DO ip = 1, nclosd
          sumip = 4*l(ip) + 2
          tm = rk(i, ip, j, ip, 0)
          ! print *, 'rk(i,ip,j,jp,0)', tm
          DO k = ABS(l(i) - l(ip)), l(i) + l(ip), 2
             tm = tm - cb(l(i), l(ip), k)*rk(i, ip, ip, j, k)
             !   print *, 'cb(l(i),l(ip),k),rk',cb(l(i),l(ip),k),rk(i,ip,ip,j,k),tm
          END DO
          hlc = hlc - 2.d0*sumip*tm
          ! print *, 'hlc',  hlc
       END DO

    END FUNCTION hlc

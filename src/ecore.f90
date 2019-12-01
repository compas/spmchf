!======================================================================
      SUBROUTINE ECORE
! =====================================================================
!     COMPUTES THE ENERGY OF THE COMMON CLOSED SHELLS
!
!----------------------------------------------------------------------
      use mchf_atomic_state
      use av_energy
      IMPLICIT NONE
      REAL(KIND=8), EXTERNAL :: rk,fk, gk, bhl

      INTEGER:: i, j, k
      REAL(KIND=8) :: TI, TIJ, sumi, sumj

      EC = 0.
      DO I = 1,NCLOSD
         SUMI = 4*L(I)+2
         TI   = FK(I,I,0)
         !Print *, ' Fk:', i, i, ti
         DO K = 2,2*L(I),2
            TI = TI - CA(L(I),K)*FK(I,I,K)
         !   print *, ' Fk:', i,i,k,ca(l(i),k), FK(I,I,0), Ti
         END DO
         EC = EC + SUMI*((SUMI-1)*TI - BHL(I,I))/2.0
         DO J = 1,I-1
            SUMJ = 4*L(J)+2
            TIJ = FK(I,J,0)
         !   print *, ' Fk:', i,j, FK(I,J,0), Tij
            DO K=IABS(L(I)-L(J)),L(I)+L(J),2
               TIJ = TIJ -CB(L(I),L(J),K)*GK(I,J,K)
         !   print *, ' Gk:', i,j,k,cb(l(i),l(j),k), GK(I,I,K), Ti
            End do
            EC = EC + SUMI*SUMJ*TIJ
         END DO
      END DO
      END SUBROUTINE ECORE

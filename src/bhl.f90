!======================================================================
      Function BHL(i, j)
!======================================================================
!
!       <P_i| L |P_j> with inclusion of rel.shift if rel = .true.
!
!----------------------------------------------------------------------

         USE spline_param
         USE mchf_atomic_state
         USE orbitals
         USE spline_hl

         IMPLICIT NONE

         INTEGER, INTENT(in) :: i, j
         REAL(KIND=8) :: BHL
         REAL(KIND=8), EXTERNAL :: BVMV

         if (iabs(L(i) - L(j)) .NE. 0) Stop ' HL:  LI <> LJ'

         Call HLM(l(i))
!      print '(A,2I5,2F10.3)', 'i,j,p(1,i),p(1,j)',i,j, p(1,i), p(1,j)
!      print *, 'p(1,i)', P(1:10,i)
!      print *, 'p(1,j)', P(1:10,j)

         BHL = BVMV(ns, ks, hl, 's', p(1, i), p(1, j))
!     print *, 'l(i), bhl', l(i), bhl
      END FUNCTION BHL


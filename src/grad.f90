!=======================================================================
    Real(8) FUNCTION grad(i,j)
!=======================================================================
!
!   <p(i)| d + [l(j)(l(j)+1)-l(i)*(l(i)+1)]/2r |p(j)>
!
!-----------------------------------------------------------------------

    USE mchf_atomic_state, ONLY: L
    USE spline_param
    USE spline_galerkin
    USE orbitals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,j

    ! .. local variables

    INTEGER :: ll,ii,jj,jp

      if ( iabs(L(I) - L(J)) .ne. 1)  Stop ' BGRAD: L(i) - L(j) <> 1'

      LL = (L(j)*(L(j)+1)-L(i)*(L(i)+1))/2

      ! .. form the sums

      grad = ll * SUM(p(:,i)*p(:,j)*rm1(:,ks))

      do jp=1,ks-1;  do ii=ks+1-jp,ns;  jj=jp+ii-ks
         grad =  grad &
           +     db1(ii,jp) * (p(jj,j)*p(ii,i) - p(jj,i)*p(ii,j)) &
           +  ll*rm1(ii,jp) * (p(jj,j)*p(ii,i) + p(jj,i)*p(ii,j))
      end do;  end do


    END FUNCTION grad



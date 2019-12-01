!=====================================================================
      subroutine lagrange
!=====================================================================
!     Compute lagrange multipliers for orbital pairs with e(i,j > 1.d-12
!
!     Each equation is  (H - eii B )P_i - eij B P_j - RHS =0
!     If we divide each equation by qsum_i, the leading term in - (hlm)P_i
!     and cancels when Lagrange multipliers are computed when qsum_i /= qsum_j
!
!     Alternatively, we could use the "sum" rule but compute lagrange multipliers
!     at the beginning of very iteration sot that the same lagrange multiplier is
!     use for each orbital pair.
!
!     ie: compute vectors HP_i and RHS_i  for orbital pairs
!     Hv(i) = H P_i -RHS_i,
!-----------------------------------------------------------------------
      USE spline_param
      USE mchf_atomic_state
      USE orbitals
      USE spline_hl
      USE av_energy
      USE angular_data, coef => value
      USE mchf_inout

      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(ns,nwf) :: hfv

      REAL(KIND=8), EXTERNAL :: fk, gk, bhl, rk
      REAL(KIND=8), DIMENSION(ns,ks) :: hd, d, dsum
      REAL(KIND=8), DIMENSION(ns) :: y
      REAL(KIND=8) :: c, cc, e11=0.d0, e22=0.d0
      INTEGER :: i, ii, j, jp, jv, jbegin, k,kv, i1,i2,i3,i4,itmp, ms

     print *, '------------- Entering lagrange ----------------'
     hfv = 0.d0; dsum =0.d0; hd=0;
     e11 = 0.d0; e22 = 0.d0

     call Icore

     call from_I_integrals

     call from_F_integrals

     call from_G_integrals

     call from_R_integrals

     print *, 'Calling compute_lagrange'
     call compute_lagrange

   CONTAINS
!=====================================================================
      subroutine Icore
!=====================================================================
!     contributions to hfv from closed shells
!---------------------------------------------------------------------
      IMPLICIT NONE

      do i = 1, nclosd
         call hlm(l(i))
         hd =  -qsum(i)*hl/2
         call bxv(ks, ns, hd,p(:,i), y)
         hfv(:,i) = hfv(:,i) +y
      end do

      do j = 1, nclosd
        ! direct k = 0
        hd = 0.d0
        call density(ns, ks, d, p(:,j), p(:,j), 's', ms)
        call add_rkm_d(0, d, hd)
        do i = 1,nclosd
          If (j /= i) then
            c = qsum(i)*qsum(j)
          else
            c = qsum(i)*(qsum(i)-1)
          end if
          call bxv(ks,ns,hd,p(:,i), y)
          hfv(:,i) = hfv(:,i) + c*y
        end do
      end do

      ! direct k >0
      hd = 0.d0
      do i = 1,nclosd
        if (l(i) > 0) then
          do k = 2, 2*l(i),2
            c = -qsum(i)*(qsum(i)-1)*ca(l(i),k)
            call density(ns, ks, d, p(:,i), p(:,i), 's', ms)
            call add_rkm_d(k,c*d,hd)
          end do
          call bxv(ks,ns,hd,p(:,i), y)
          hfv(:,i) = hfv(:,i) + y
        end if
      end do

      ! exchange
      do i = 1,nclosd
        do j = i+1,nclosd
          hd = 0.d0
          call density(ns,ks,d,p(:,j), p(:,i), 's', ms)
          do k = abs(l(i)-l(j)), l(i)+l(j), 2
            c = -qsum(i)*qsum(j)*cb(l(i), l(j),k)
            call add_rkm_d(k, d, hd)
          end do
          call bxv(ks,ns,hd,p(:,j), y)
          hfv(:,i) = hfv(:,i) + c*y
          call bxv(ks,ns,hd,p(:,i), y)
          hfv(:,j) = hfv(:,j) + c*y
        end do
     end do
     end subroutine Icore
!=====================================================================
      subroutine from_I_integrals
!=====================================================================
!     contributions to hfv (from diagonal ) I(i,j) integrals
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ip
      REAL(KIND=8) :: cc


      do jp = intptr(5)+1, intptr(6)
        if (.not. lused(jp)) cycle
        call unpacki(6,jp, kv, i1, i2, i3, i4)
        ! contributions ot hfv, i=1,nclose
        if ( nclosd > 0 ) call Icoreval

        ! contribution to i1 and i2
        call Ivalcore
       end do

       end subroutine from_I_integrals
!=====================================================================
      subroutine Icoreval
!=====================================================================
!     contributions to hfv for orbitals in closed shells
!     related to outer I-integrals (i1,i2)  (i < nclosd)
!---------------------------------------------------------------------
      IMPLICIT NONE

      ! direct k = 0   Rk(.,i1;.,i2,0)
      ! effect on core
      hd = 0.d0
      call density(ns, ks, d, p(:,i1), p(:,i2), 's', ms)
      call add_rkm_d(0, d, hd)
      call bxv(ks,ns,hd,p(:,1), y)
      do i = 1,nclosd
         c = -2*coef(jp)*qsum(i)
         call bxv(ks,ns,hd,p(:,i), y)
         hfv(:,i) = hfv(:,i) + c*y
      end do

      ! exchange        Rk(i,i1;i2,i)

      do i = 1, nclosd
         call density(ns,ks,d,p(:,i1), p(:,i), 's', ms)
         do k = abs(L(i)-L(i1)), L(i)+L(i1), 2
            hd = 0.0d0
            c= 2*coef(jp)*qsum(i)*cb(l(i), l(i1),k)
            call add_rkm_d(k, d, hd)
            call bxv(ks,ns,hd,p(:,i2), y)
            hfv(:,i) = hfv(:,i) + c*y
         end do
      end do
     end subroutine Icoreval
!=====================================================================
     subroutine Ivalcore
!=====================================================================
!     contributions to hfv for orbital i outside closed shells from the
!     core related to outer I-integrals (ii,i2) where i1 = i2 =i > nclosd
!---------------------------------------------------------------------
     IMPLICIT NONE
     REAL(KIND=8)        :: fk, gk

     call hlm(l(i1))
     c = coef(jp)
     if (i1 /= i2) c = c/2
     call bxv(ks,ns, hl, p(:,i2), y)
     hfv(:,i1) = hfv(:,i1) +c*y
     if ( i1==3) then
        print *, 'I(2p)', c, dot_product(hfv(:,i1), p(:, i1)),-c*bhl(i1,i1)/2
     end if
     if (i1 /=i2) then
       call bxv(ks,ns, hl, p(:,i1), y)
       hfv(:,i2) = hfv(:,i2) +c*y
     end if

     if (nclosd > 0) then

        dsum = 0.d0
        do j = 1, nclosd
           ! direct k = 0   Rk(i1,j;i2,j,0)
           call density(ns, ks, d, p(:,j), p(:,j), 's', ms)
           dsum = dsum + d*qsum(j)
        end do
        hd = 0.0d0
        call add_rkm_d(0, dsum, hd)
        c = -coef(jp)
        if (i1 == i2) c = 2*c
        call bxv(ks,ns,hd,p(:,i2), y)
        hfv(:,i1) = hfv(:,i1) + c*y
     if ( i1==3) then
        print *, 'Fk', c,  dot_product(hfv(:,i1), p(:, i1)), &
               c*bhl(i1,i1)/2+4*fk(1,3,0)
     end if

        if (i1 /= i2) then
          call bxv(ks,ns,hd,p(:,i1), y)
          hfv(:,i2) = hfv(:,i2) + c*y
        end if

      end if

      ! exchange        Rk(i1,j;j,i2)
     do ii = 1,2
     do j = 1,nclosd
        call density(ns,ks,d,p(:,j), p(:,i2), 's', ms)
        do k = abs(l(j)-l(i1)), l(j)+l(i1), 2
          hd = 0.d0
          c=coef(jp)*qsum(j)*cb(l(i),l(j),k)
          if (i1 == i2 ) c = 2*c
          call add_rkm_d(k, d, hd)
          call bxv(ks,ns,hd,p(:,j), y)
          hfv(:,i1) = hfv(:,i1) + c*y
        end do
     end do

     if (i1 /= i2) then
        ! interchange i1 and i2
        i3 = i1; i1=i2; i2=i3
     else
        exit
     end if

     end do

     end subroutine Ivalcore
!=====================================================================
      subroutine from_F_integrals
!=====================================================================
!     contributions to hfv F-integrals that are stored in order of
!     increasing k
!---------------------------------------------------------------------
      IMPLICIT NONE

       jbegin=1
       do jp = jbegin,intptr(1)
         if (.not. lused(jp)) cycle
         call unpacki(1,jp, kv, i1, i2, i3, i4)
           !   Rk(i1,i2; i1, i2)
           c = coef(jp)
           if (i1 == i2) c = 2*c
           do ii = 1,2
              hd = 0.d0
              call density(ns,ks,d,p(:,i2),p(:,i2), 's', ms)
              call add_rkm_d(kv,d, hd)
              call bxv(ks,ns,hd,p(:,i1), y)
              hfv(:,i1) = hfv(:,i1) + c*y
              if (i1 == i2) then
                 exit
              else
                i3 = i1; i1=i2; i2=i3
              end if
           end do
       end do
     end subroutine from_F_integrals
!=====================================================================
      subroutine from_G_integrals
!=====================================================================
!     contributions to hfm from G-integrals
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: im

     ! Add exchange contribution
      jbegin=intptr(1)+1
      do jp = jbegin, intptr(2)
       if (.not. lused(jp)) cycle
       call unpacki(2,jp, kv, i1, i2, i3, i4)
           !   Rk(i1,i2; i2, i1)
           c = coef(jp)
           hd = 0.d0
           call density(ns,ks,d,p(:,i2),p(:,i1), 's', ms)
           call add_rkm_d(kv,d, hd)
           call bxv(ks,ns,hd,p(:,i2), y)
           hfv(:,i1) = hfv(:,i1) + c*y
           call bxv(ks,ns,hd,p(:,i1), y)
           hfv(:,i2) = hfv(:,i2) + c*y
       end do
     end subroutine from_G_integrals
!=====================================================================
      subroutine from_R_integrals
!=====================================================================
!     contributions to hfm R-integrals
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(ns) :: v
      INTEGER :: ie1, ie2

     ! Direct contributions
     jbegin=intptr(4)+1
     do jp = jbegin, intptr(5)
       if (.not. lused(jp)) cycle
       call unpacki(5,jp, kv, i1, i2, i3, i4)
       do ii = 1,2
           hd = 0.d0
           call density(ns, ks, d, p(:,i2), p(:, i4), 's', ms)
           call add_rkm_d(kv, d, hd)
           c = coef(jp)
           if (i1 /= i3) c = c/2
           if (i1 == i2 .and. i3 == i4 ) c = 2*c
           call bxv(ks, ns, hd, p(:,i3), y)
           hfv(:, i1) = hfv(:,i1) + c*y
           if (i1 /= i3) then
             call bxv(ks, ns, hd, p(:,i1), y)
             hfv(:, i3) = hfv(:,i3) + c*y
           end if
           if (i1 ==i2 .and. i3 ==i4) then
              exit
           else
              ! interchange
             itmp = i1; i1 = i2; i2=itmp
             itmp = i3; i3 = i4; i4=itmp
           end if
        end do
     end do


     end subroutine from_R_integrals

!=====================================================================
      subroutine compute_lagrange
!=====================================================================
!     compute e(i,j) matrix
!---------------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(nwf,nwf) :: eij
      REAL(KIND=8) :: eav, e1a, e1b, e1c, e2a, e2b, e11, e22

      eij = 0.d0
      do i = 1,nwf
        do j = 1, nwf
         if ( i /= j .and. (clsd(i) .and. clsd(j)) ) cycle
          if (l(i) == l(j)) eij(i,j) = dot_product(hfv(:,i), p(:,j))
        end do
      end do
         write(50,*) 'Energy matrix before symmetrization: SCF'
      do i = 1,nwf
         write(50, '(10F15.8)') eij(i,:)
      end do
      !write(scr,FMT='(/7X,A)') 'Energy matrix before symmetrization'
      !do i = 1,nwf
      !   write(scr, '(7X,A,I4,(8F12.5))') 'Row i=', i, eij(i,:)
      !end do
      !write(scr,*)

      ! Symmetrize the energy matrix and deal with closed shells
      Do i = 1,nwf
         Do j = i, nwf
            if (clsd(i) .and. clsd(j)) cycle
               if (l(i) == l(j)) then
                 eav = (eij(i,j)+eij(j,i))/2
                 eij(i,j) = eav
                 eij(j,i) = eav
                 ! save these values
                 e(i,j) = eav
                 e(j,i) = eav
            end if
          end do
          e(i,i) = eij(i,i)
      end do

     ! print *, 'Energy matrix After symmetriztion and closed shell consideration'
     ! do i = 1,nwf
     !   Print '(A,I4,(8F15.8))', 'Row i=', i, eij(i,:)
     !end do

    end subroutine compute_lagrange

   END SUBROUTINE lagrange

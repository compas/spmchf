!=====================================================================
      subroutine scf_matrix(i,hfm, hx, rhs)
!=====================================================================
!
!     1) Set up the hf_matrix for orbital i.  Banded matrices will be
!       summed in the array hd, converted to full hfm, and exchange
!       contributions added to hfm
!     2) hx is the matrix of extra terms for hf_nr
!     3) rhs is the rhs for orbitals that appear only once in an integral
!
!      Intgrl   Contribution to hfm      hx	Contribution to RHS
!      F(i,i)     2[ii]                 4(ii)
!      F(i,j)      [jj]
!      G(i,j)      (jj)
!      I(i.i)      I(., .)
!      I(i,j)                                   (1/2) I(.,.)j   + core ?
!      R(i,i;i,i)  [ii]                 4(ii)
!      R(i,i;i,j)  (1/2)([ij]+{i,j})  same as hfm
!      R(i,i;a,b)  (1/2){ab}          same as hfm
!      R(i,a;i,b)  [ab]
!      R(i,a;b,c)                                 (1/2)R(.,a;.,c)b
!
! Example:
!<  2 |H|  1 > = < 2s( 1) 2p( 2) 3s( 1) |H|  2s( 2) 2p( 2) >
!
!
!               COEFFICIENTS OF VARIOUS INTEGRALS
!
!   E-TOTAL(ET)  E-AVERAGE(EAV)  (ET - EAV)    INTEGRALS WITH OVERLAPS
!
!     1.414214      0.000000      1.414214    I ( 3s, 2s)
!     0.707107      0.000000      0.707107    R 0( 2s, 3s; 2s, 2s)
!     0.707107      0.000000      0.707107    R 0( 2s, 3s; 2s, 2s)
!    -0.471405      0.000000     -0.471405    R 1( 2p, 3s; 2s, 2p)
!     2.828427      0.000000      2.828427    R 0( 2p, 3s; 2p, 2s)
!     2.828427      0.000000      2.828427    R 0( 1s, 3s; 1s, 2s)
!    -1.414214      0.000000     -1.414214    R 0( 1s, 3s; 2s, 1s)
!
!-----------------------------------------------------------------------
      USE spline_param
      USE mchf_atomic_state
      USE orbitals
      USE spline_hl
      USE av_energy
      USE angular_data, coef => value

      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(ns,ns),INTENT(out) :: hfm, hx
      REAL(KIND=8), DIMENSION(ns),INTENT(out)    :: rhs
      INTEGER, INTENT(IN) :: i

      REAL(KIND=8), EXTERNAL :: fk, gk, bhl, rk
      REAL(KIND=8), DIMENSION(ns,ks) :: hd, d, dsum, ihd
      REAL(KIND=8), DIMENSION(ns,ns) :: dx,dxsum, ihfm
      REAL(KIND=8), DIMENSION(ns) :: y
      REAL(KIND=8) :: c, cc, e11, rhs1, eii
      INTEGER :: j, jp, jv, jbegin, k, kv, ie, i1,i2,i3,i4,itmp, ms
      LOGICAL :: found, foundx

     hfm = 0.d0; rhs=0.d0; dsum =0.d0; dxsum=0.d0; hd=0; hx = 0
     ihd=0; ihfm=0

      e11=0; rhs1=0

      call hlm(l(i))

      if (i <= nclosd) call Icore(i)

     call from_I_integrals
     call bxv(ks,ns,hd, p(:,i),y)
     print *, 'From I_integrals', dot_product(p(:,i),y)
     y = matmul(hfm, p(:,i))
     print *, 'From I_integrals: exchage', dot_product(p(:,i),y)


     call from_F_integrals
     call bxv(ks,ns,hd, p(:,i),y)
     print *, 'From F_integrals', dot_product(p(:,i),y)

     call from_G_integrals
     y = matmul(hfm, p(:,i))
     print *, 'From G_integrals', dot_product(p(:,i),y)

     call from_R_integrals
     call bxv(ks,ns,hd, p(:,i),y)
     print *, 'From R_integrals:hd', dot_product(p(:,i),y)
     y = matmul(hfm, p(:,i))
     print *, 'From R_integrals:hx ', dot_product(p(:,i),y)


     Call add_band_to_full(hd,hfm)
     y = matmul(hfm, p(:,i))
     print *, ' From hfm', dot_product(p(:,i), y)

!
!
!     call from_I_integrals
!
!     call from_F_integrals
!
!     call from_G_integrals
!
!     call from_R_integrals
!
!     Call add_band_to_full(hd,hfm)

   CONTAINS
!=====================================================================
      subroutine Icore(i)
!=====================================================================
!     contributions to hfm for orbital i in the closed shells
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: i

      hd =  hd -qsum(i)*hl/2

      ! direct  (k=0)
      dsum = 0.d0
      do j = 1, nclosd
        ! direct k = 0
        call density(ns, ks, d, p(:,j), p(:,j), 's', ms)
        If (j /= i) then
          c = qsum(i)*qsum(j)
        else
          c = qsum(i)*(qsum(i)-1)
        end if
        dsum = dsum +c*d
        if (j ==i) then   !  compute the contribution to hx
          c = 2*c
          call density(ns,ns,dx, p(:,j), p(:,j), 'x', ms)
          call add_rkm_x(0,c*dx, hx)
        end if
      end do
      call add_rkm_d(0, dsum, hd)

      ! direct k >0
      if (l(i) > 0) then
        call density(ns, ks, d,  p(:,i), p(:,i), 's', ms)
        call density(ns, ns, dx, p(:,i), p(:,i), 'x', ms)
        do k = 2, 2*l(i),2
          c = -qsum(i)*(qsum(i)-1)*ca(l(i),k)
          call add_rkm_d(k,c*d,hd)
          call add_rkm_x(k,2*c*dx, hx)
        end do
      end if

      ! exchange
      do k = 0, kmax
        dxsum = 0.d0
        found = .false.
        do j = 1,nclosd
          if (i ==j) cycle
          if ( k > l(i)+l(j) .or. k < abs(l(i)-l(j)) ) cycle
          if ( mod( abs(l(i)-l(j))-k,2) /= 0) cycle
          call density(ns,ns,dx,p(:,j), p(:,j), 'x', ms)
          c = -qsum(i)*qsum(j)*cb(l(i), l(j),k)
          dxsum = dxsum +c*dx
          found = .true.
        end do
        if (found) call add_rkm_x(k, dxsum, hfm)
     end do
     end subroutine Icore
!=====================================================================
      subroutine from_I_integrals
!=====================================================================
!     contributions to hfm (from diagonal ) and rhs (from non-diagonl)
!     I(i,j) integrals
!---------------------------------------------------------------------
! <  1 |H|  1 > = < 3s( 2) 3p( 2) |H|  3s( 2) 3p( 2) >
!
!                COEFFICIENTS OF VARIOUS INTEGRALS
!
!    E-TOTAL(ET)  E-AVERAGE(EAV)  (ET - EAV)    INTEGRALS WITH OVERLAPS
!
!      2.000000      2.000000      0.000000    I ( 3s, 3s)
!      2.000000      2.000000      0.000000    I ( 3p, 3p)
!      2.000000      2.000000      0.000000    I ( 1s, 1s)
!      2.000000      2.000000      0.000000    I ( 2s, 2s)
!      6.000000      6.000000      0.000000    I ( 2p, 2p)
!      1.000000      1.000000      0.000000    F 0( 3s, 3s)
!      4.000000      4.000000      0.000000    F 0( 3s, 3p)
!     -0.666667     -0.666667      0.000000    G 1( 3s, 3p)
!      1.000000      1.000000      0.000000    F 0( 3p, 3p)
!      0.400000     -0.080000      0.480000    F 2( 3p, 3p)
!      1.000000      1.000000      0.000000    F 0( 1s, 1s)
!      4.000000      4.000000      0.000000    F 0( 1s, 3s)
!     -2.000000     -2.000000      0.000000    G 0( 1s, 3s)
!      4.000000      4.000000      0.000000    F 0( 1s, 3p)
!     -0.666667     -0.666667      0.000000    G 1( 1s, 3p)
!      4.000000      4.000000      0.000000    F 0( 1s, 2s)
!     -2.000000     -2.000000      0.000000    G 0( 1s, 2s)
!      1.000000      1.000000      0.000000    F 0( 2s, 2s)
!      4.000000      4.000000      0.000000    F 0( 2s, 3s)
!     -2.000000     -2.000000      0.000000    G 0( 2s, 3s)
!      4.000000      4.000000      0.000000    F 0( 2s, 3p)
!     -0.666667     -0.666667      0.000000    G 1( 2s, 3p)
!     12.000000     12.000000      0.000000    F 0( 1s, 2p)
!     -2.000000     -2.000000      0.000000    G 1( 1s, 2p)
!     12.000000     12.000000      0.000000    F 0( 2s, 2p)
!     -2.000000     -2.000000      0.000000    G 1( 2s, 2p)
!     15.000000     15.000000      0.000000    F 0( 2p, 2p)
!     -1.200000     -1.200000      0.000000    F 2( 2p, 2p)
!     12.000000     12.000000      0.000000    F 0( 2p, 3s)
!     -2.000000     -2.000000      0.000000    G 1( 2p, 3s)
!     12.000000     12.000000      0.000000    F 0( 2p, 3p)
!     -2.000000     -2.000000      0.000000    G 0( 2p, 3p)
!     -0.800000     -0.800000      0.000000    G 2( 2p, 3p)
!
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ip
      REAL(KIND=8) :: cc


      do jp = intptr(5)+1, intptr(6)
        if (.not. lused(jp)) cycle
        call unpacki(6,jp, kv, i1, i2, i3, i4)

        if ( i <= nclosd) then
            call Icoreval(i)
         else if (i1 == i .and.  i2 == i) then
            call Ivalcore(i)
         else if (i1 ==i .or. i2 ==i) then
            call Irhs(i)
         end if
       end do

       end subroutine from_I_integrals
!=====================================================================
      subroutine Icoreval(i)
!=====================================================================
!     contributions to hfm for orbital i in the closed shells
!     related to outer I-integrals (i1,i2)  (i < nclosd)
!---------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i

      ! direct k = 0   Rk(.,i1;.,i2,0)
      if (nclosd > 0) then
          call density(ns, ks, d, p(:,i1), p(:,i2), 's', ms)
          c = -2*coef(jp)*qsum(i)
          call add_rkm_d(0, c*d, hd)
      end if

      ! exchange        Rk(.,.;i2,i1)
      do k = abs(L(i)-L(i1)), L(i)+L(i1), 2
         call density(ns,ns,dx,p(:,i2), p(:,i1), 'x', ms)
         c= 2*coef(jp)*qsum(i)*cb(l(i), l(i1),k)
         if (i1 /= i2) dx = (dx+transpose(dx))/2
         call add_rkm_x(k, c*dx, hfm)
      end do
     end subroutine Icoreval
!=====================================================================
     subroutine Ivalcore(i)
!=====================================================================
!     contributions to hfm for orbital i outside closed shells from the
!     core related to outer I-integrals (ii,i2) where i1 = i2 =i > nclosd
!---------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i
     REAL(KIND=8)        :: fk, gk

     hd = hd + coef(jp)*hl
     call bxv(ks,ns,hd, p(:,i),y)

     dsum = 0.d0
     do j = 1, nclosd
        ! direct k = 0   Rk(.,j;.,j,0)
        call density(ns, ks, d, p(:,j), p(:,j), 's', ms)
        c = -2*coef(jp)*qsum(j)
        dsum = dsum +c*d
      end do
      if (nclosd > 0) call add_rkm_d(0, dsum, hd)

      ! exchange        Rk(.,.;i2,i1)
      do k = 0, kmax
        dxsum = 0.d0
        found = .false.
        do j = 1,nclosd
          if (k > l(i)+l(j)) cycle
          if ( mod( abs(l(i)-l(j))-k,2) /= 0) cycle
          call density(ns,ns,dx,p(:,j), p(:,j), 'x', ms)
          c=2*coef(jp)*qsum(j)*cb(l(i),l(j),k)
          dxsum = dxsum + c*dx
          found = .true.
        end do
        if (found) call add_rkm_x(k, dxsum, hfm)
     end do

     end subroutine Ivalcore
!=====================================================================
     subroutine Irhs(i)
!=====================================================================
!     contributions to rhs for orbital i outside closed shells from the
!     core related to outer I-integrals (ii,i2) where i1 /= i2 .and.
!     i > nclosd
!---------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: i
     REAL(KIND=8), DIMENSION(ns) :: v

     if (i1 /= i) then
        i3 = i1; i1 = i2; i2 = i3   ! interchange
     end if
     c = coef(jp)/2
     ihd = ihd + c*hl

     c = -2*c
     dsum = 0.d0
     do j = 1, nclosd
        ! direct k = 0   Rk(.,j;.,j,0)
        call density(ns, ks, d, p(:,j), p(:,j), 's', ms)
        dsum = dsum +c*qsum(j)*d
     end do

     if (nclosd > 0) call add_rkm_d(0, dsum, ihd)
     call bxv(ks,ns,ihd,p(:,i2),v)
     rhs = rhs + v

     ! exchange        Rk(.,.;j,j)
     do k = 0, kmax
        dxsum = 0.d0
        foundx = .false.
        do j = 1,nclosd
          if ( mod( abs(l(i)-l(j))-k,2) /= 0) cycle
          call density(ns,ns,dx,p(:,j), p(:,j), 'x', ms)
          dxsum = dxsum - c*qsum(j)*cb(l(i), l(j),k)*dx
          foundx = .true.
        end do
        if (foundx) call add_rkm_x(k, dxsum, ihfm)
     end do
     If (foundx) rhs = rhs + matmul(ihfm,p(:,i2))

     end subroutine Irhs
!=====================================================================
      subroutine from_F_integrals
!=====================================================================
!     contributions to hfm F-integrals that are stored in order of
!     increasing k
!---------------------------------------------------------------------
      IMPLICIT NONE

       k = 0; jbegin=1; found = .false.; dsum = 0.d0
       do jp = jbegin,intptr(1)
         if (.not. lused(jp)) cycle
         call unpacki(1,jp, kv, i1, i2, i3, i4)
         if (kv >k) then
            if (found) call add_rkm_d(k, dsum, hd)
            dsum = 0
            found = .false.
            k = kv
         end if
         ie = 0
         If (i1 == i) then
            ie = i2
         else if (i2 == i) then
            ie = i1
         end if
         if (ie /= 0) then
           c = coef(jp)
           if (ie == i) c = 2*c
           call density(ns,ks,d,p(:,ie),p(:,ie), 's', ms)
           dsum = dsum +c*d
           found = .true.
           if (ie ==i) then
              call density(ns,ns,dx,p(:,ie), p(:,ie), 'x', ms)
              call add_rkm_x(k,2*c*dx, hx)
           end if
         end if
       end do
       if (found) call add_rkm_d(k, dsum, hd)
       call bxv(ks,ns,hd,p(:,i),y)

     end subroutine from_F_integrals
!=====================================================================
      subroutine from_G_integrals
!=====================================================================
!     contributions to hfm from G-integrals
!---------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: im

     ! Add exchange contribution
      dxsum = 0.d0
      k = 0; found = .false.
      jbegin=intptr(1)+1
      do jp = jbegin, intptr(2)
       if (.not. lused(jp)) cycle
       call unpacki(2,jp, kv, i1, i2, i3, i4)

       if (kv >k) then
          if (found) call add_rkm_x(k, dxsum, hfm)
          found = .false.
          dxsum = 0.d0
          k = kv
       end if

       ie = 0
       if (i1 == i) then
          ie = i2
       else if (i2 == i) then
          ie = i1
       end if

       if (ie /= 0) then
         c = coef(jp)
         call density(ns,ns,dx,p(:,ie),p(:,ie), 'x', ms)
         dxsum = dxsum +c*dx
         found = .true.
       end if
     end do

     if (found) call add_rkm_x(k,dxsum,hfm)

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
     k =0
     dxsum = 0.d0; dsum  = 0.d0
     found = .false. ; foundx = .false.
     jbegin=intptr(4)+1
     do jp = jbegin, intptr(5)
       if (.not. lused(jp)) cycle
       call unpacki(5,jp, kv, i1, i2, i3, i4)
        Print '(A,5I4,F9.5,I4)', 'R: jp, i1, i2, i3, i4', &
                               jp, i1, i2, i3, i4, coef(jp), kv
       !if (kv >k) then
       !   if (found) call add_rkm_d(k,dsum,hd)
       !   if (foundx) call add_rkm_x(k,dxsum,hfm)
       !   k =kv
       !  dsum=0.d0; dxsum=0.d0; found=.false.; foundx=.false.
       !end if
       c = coef(jp)

       ! re-order in a standard order
       ie1 = 0; ie2=0
       if ((i1 == i) .or. (i3 == i)) then
          ie1 = 1
          if (i1 /= i) then
             itmp = i1; i1=i3; i3 = itmp
           end if
           if ((i1 == i) .and. (i3 == i)) ie1 = 2
       end if
       if ((i2 == i) .or. (i4 == i) ) then
          ie2 = 1
          if (i2 /= i) then
             itmp = i2; i2=i4; i4 = itmp
           end if
           if ((i2 ==i) .and. (i4 == i) ) ie2 = 2
        end if
       if (ie1+ie2 == 0) cycle

       if ( ie2 > ie1 ) then
           itmp = i1; i1 = i2; i2 = itmp
           itmp = i3; i3 = i4; i4 = itmp
           itmp = ie1; ie1=ie2; ie2=itmp
       end if

       select case (ie1+ie2)                  ! Rk(i,i2;i3,i4)
         case (1)
           call density(ns,ks,d,p(:,i2),p(:,i4),'s',ms)
           ihd = 0
           call add_rkm_d(k, d, ihd)
           call bxv(ks,ns,ihd,p(:,i3),v)
           rhs = rhs + (c/2)*v

          case (2)             ! Rk(i,i2;i.i4) or Rk(i,i;i3,i4)
            if ( i3==i ) then
               call density(ns,ks,d,p(:,i2),p(:,i4),'s',ms)
               call add_rkm_d(kv, c*d, hd)
            else if ( i2 ==i )  then
               print *, 'case(2), exchange', i3,i4
               call density(ns,ns,dx,p(:,i3),p(:,i4),'x',ms)
               if (i3 /= i4) dx = (dx + transpose(dx))/2
               call add_rkm_x(kv, c*dx, hfm)
               v=matmul(hfm, p(:,i))
               print *, 'R_int', dot_product(v,p(:,i))
            end if

          case (3)             ! Rk(i,i;i,i4)
            call density(ns,ks,d,p(:,i2),p(:,i4),'s',ms)
            call add_rkm_d(kv, c*d, hd)

            call density(ns,ks,d,p(:,i1),p(:,i3),'s',ms)
            ihd = 0.d0
            call add_rkm_d(kv, d, ihd)
            call bxv(ks, ns, ihd, p(:,i4), v)
            rhs = rhs +(c/2)*v
            print *, 'R_rhs:case(3) ', dot_product(rhs, p(:,i)), &
                      (c/2)*rk(i,i,i,i4,kv)

            call density(ns,ns,dx,p(:,i3),p(:,i4),'x',ms)
            dx = dx + transpose(dx)
            call add_rkm_x(kv, c*dx, hx)



          case (4)                 !Rk(i,i;i,i)
            call density(ns,ks,d,p(:,i2),p(:,i4),'s',ms)
            call add_rkm_d(kv, 2*c*d, hd)
            call add_rkm_x(kv, 4*c*dx, hx)
        end select
     end do

     end subroutine from_R_integrals

   END SUBROUTINE scf_matrix

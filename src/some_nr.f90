!==================================================================
       SUBROUTINE some_nr(n1,n2)
!==================================================================
!    Improve the estimates v by the Newton_Raphson method
!    subject to orthogonality where:
!      . orbitals 1:nf are fixed orbitals (not varied)
!      . orbitals nf+1:nv  are varied,
!      . orbitals nv+1:nwf are fixed
!    Matrices:
!      . hd -- second order variations (banded matrix)
!         . hda(ns,ks,1:nf,nf+1:nv)  (interaction of first fixed
!                                     with varied)
!         . hdb(ns,ks,nf+1:nv,nf+1:nwf) (interaction of varied with
!                              varied and remaining fixed orbitals)
!      . hx -- second order variations (full matrix)
!         . aa (nf+1:nv,  nf+1,nv)
!    Equation of orbital v
!      . hm(ns,ns,nf+1:nv,1:nwf) = transpose(hda(ns,ns,1:nf,nf+1,nv)
!                                  + hdb(ns,ns,nf+1:nwf)
!------------------------------------------------------------------
      USE angular_data, coef=>value
      USE av_energy
      USE spline_hl
      USE spline_param
      USE spline_galerkin
      USE mchf_atomic_state
      USE mchf_inout, ONLY: err
      USE mchf_param, ONLY: end_tol
      USE orbitals

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n1, n2

      REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: hda, hdb, hx, hm
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE      :: rhs
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE       :: res, work_nr
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE      :: aa
      INTEGER, DIMENSION(:), ALLOCATABLE            :: ipiv
      INTEGER :: iwork_nr


      REAL(KIND=8), EXTERNAL               :: a, b, bvmv, bhl, fk, gk, rk
      REAL(KIND=8), DIMENSION(ns)          :: x, y, v
      REAL(KIND=8), DIMENSION(ns,ns)       :: dx, dxsum
      REAL(KIND=8), DIMENSION(ns,ks)       :: d, dsum
      REAL(KIND=8), DIMENSION(n1+1:n2,nwf) :: eij
      REAL(KIND=8)                         ::  c, cc, ccc, e11, e22

      INTEGER :: mm, k, ms, info, i, j, ip, jp, nns, nd, nnd, north, &
                 i1, i2, i3, i4, kv, ibegin, iend, itmp, ic, &
                 ib, ie, jb, je, ii, jj

      allocate(hda(ns,ks,1:n1,n1+1:n2))                ! hda
      allocate(hdb(ns,ks,n1+1:n2,n1+1:nwf))             ! hdb
      allocate(hx(ns, ns, n1+1+n2, n1+1:n2))           ! hx
      allocate(hm(ns,ns,n1+1:n2, n1+1:ns))             ! hm
      allocate(rhs(ns, n1+1:n2))                      ! rhs
      hda = 0; hdb = 0; hx = 0; hm=0; rhs=0

!     determine orthogonality constraints between orbitals varied
      north = 0
      Do i = n1+1, n2
        Do j = i+1, n2
           if (clsd(i) .and. clsd(j)) cycle
           if ( l(i) == l(j) .and. abs(e(i,j)) >1.d-12) then
              north = north +1
           end if
        End do
      End do

      nns = (ns)*(n2-n1)
      nnd = nns + (n2-n1)
      nd = nnd + north
      iwork_nr = 3*nd
      allocate (aa(nd,nd), res(nd), work_nr(iwork_nr), ipiv(nd))
      aa=0; res=0; work_nr=0; ipiv=0
      aa=0; res=0; ipiv=0;

      ! create the gradient (hd)  and Hessian (hm+hx) matrix

      ! contributions from the closed shells
      call closed_shells


      ! contributions from I-integrals
      call I_integrals

      ! contributions for F-integrals
      call F_integrals

      ! contributions for G-integrals
      call G_integrals

      ! contributions from R-integrals
      call R_integrals

      ! compute lagrange multipliers, residuals and modify orbital equations
      Call hm_residuals

      ! Add normalization constraints
      Call hm_norm

      ! Add orthogonality constraints
      Call hm_orthog

      ! check the  matrix
      !call hm_check

       ! solve the Newton Raphson hm equations
      Call hm_solve

      ! update and test the solutions
      Call hm_update

      deallocate(hda, hdb, hx, hm, rhs)
      deallocate(aa, res, work_nr, ipiv)
   CONTAINS

!==================================================================
      SUBROUTINE closed_shells
!==================================================================
!     Contributions to hd and hdm from closed shells (upper only)
!------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: kp, k

      ! upper matrix only, all direct
      do i = 1, nclosd
        do j = i,nclosd
          if (i ==j) then
             if (i <=n1 .or. i > n2) cycle
             call hlm(l(i))
             hdb(:,:,i,i) = hdb(:,:,i,i) -qsum(i)*hl/2
             call density(ns,ks,d,p(:,i), p(:, i), 's', ms)
             call density(ns,ns,dx,p(:,i), p(:, i), 'x', ms)
             c = qsum(i)*(qsum(i)-1)
             ! direct
             do k = 0, 2*l(i), 2
               if (k == 0) then
                  cc = c
               else
                  cc = -c*ca(l(i),k)
               end if
               call add_rkm_d(k,cc*d,hdb(:,:,i,i))
               call add_rkm_x(k,2*cc*dx,hx(:,:,i,i))  ! i=j

             end do
           else
             ! direct    F0(i,j)  i <j, k=0 of j < ib < i
             k=0
             c = qsum(i)*qsum(j)
             call add_cont('d',c, k, j, j, i, i, 'Closed')
             call add_cont('d',c, k, i, i, j, j, 'Closed')
             call add_cont('x',2*c,k,i, j, i, j, 'Closed')

             !  exchange Gk(i,j)
             do k = abs(l(i)-l(j)), l(i)+l(j),2
                c = -qsum(i)*qsum(j)*cb(l(i), l(j), k)

                call add_cont('d', c, k, i, j, i, j, 'Closed')
                call add_cont('x', c, k, j, j, i, i, 'Closed')
                call add_cont('x', c, k, j, i, i, j, 'Closed')
                call add_cont('x', c, k, i, i, j, j, 'Closed')

             end do
           end if
        end do
      end do

      end subroutine closed_shells

!==================================================================
      SUBROUTINE I_integrals
!==================================================================
!     Contributions to hd and hx from I_integrals  (upper matrix)
!                                                   (i1 < i2)
!------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: kv, k

      ibegin = intptr(5)+1
      iend   = intptr(6)
      do ip = ibegin, iend
        if (.not. lused(ip)) cycle
        call unpacki(6,ip, kv, i1, i2, i3, i4)
        Print '(A,5I4,F9.5)', 'I: ip, i1, i2, i3, i4', &
                               ip, i1, i2, i3, i4, coef(ip)

        call hlm(l(i1))
        c = -2*coef(ip)
        cc = c
        if (i1 /= i2) cc = cc/2
        if (i1 > n1 .and. i1 <= n2 ) then
           hdb(:,:,i1,i2) = hdb(:,:,i1,i2) -cc*hl/2
        else if (i1 <= n1 .and. i2 <= n2 ) then
           hda(:,:,i1,i2) = hda(:,:,i1,i2) -cc*hl/2
        end if

!       !  Direct contributions to (i1,i2) from core
         dsum = 0.d0
         Do jp = 1, nclosd     ! Rk(jp,i1;jp,i2) == Rk(i1,jp;i2,jp)
           ! direct
           call density(ns,ks,d,p(:,jp), p(:,jp), 's', ms)
           dsum = dsum + cc*qsum(jp)*d
        end do
        if (i1 > n1 .and. i1 <=n2 ) then
          call add_rkm_d(0,dsum,hdb(:, :, i1,i2))        ! hd(i1,i2)
        else if (i1 <= n1 .and. i2<=n2 ) then
          call add_rkm_d(0, dsum, hda(:,:, i1,i2))
        end if

         ! Direct contribution to core from (i1,i2)
        Do jp = 1, nclosd     ! Rk(jp,i1;jp,i2) == Rk(i1,jp;i2,jp)

          cc = qsum(jp)*c
          call add_cont('d', cc,  0, i1, i2, jp, jp, 'Closed_I')   !hd(jp, jp)
          if (i1 == i2) then
             call add_cont('x',2*cc, 0, jp, i2, jp, i1, 'Closed_I')   !hx(jp,i1)
          else
           call add_cont('x', cc, 0, jp, i2, jp, i1,'Closed_I')   !hx(jp,i1)
           call add_cont('x', cc, 0, jp, i1, jp, i2,'Closed_I')   !hx(jp,i2)
          end if
        end do

      !  exchange contributions  from  Rk(jp, jp; i1, i2)
        do jp = 1,nclosd                             !   Rk(i1,jp;jp,i2)
          do k = abs(l(i1)-l(jp)), l(i1)+l(jp), 2

            cc = 2*cb(l(i1), l(jp), k)*qsum(jp)*coef(ip)

            if (i1 == i2) then
              call add_cont('d', cc, k, jp, i1, jp, i1, 'Closed_I')
              call add_cont('x', cc, k, i1, i1, jp, jp, 'Closed_I')
              call add_cont('x', cc, k, i1, jp, jp, i1, 'Closed_I')
              call add_cont('x', cc, k, jp, jp, i1, i1, 'Closed_I')

            else
              cc = cc/2
              call add_cont('d', cc, k, jp, i2, jp, i1, 'Closed_I')
              call add_cont('d', cc, k, jp, i1, jp, i2, 'Closed_I')
              call add_cont('s', cc, k, i1, i2, jp, jp, 'Closed_I')
              call add_cont('x', cc, k, i2, jp, jp, i1, 'Closed_I')
              call add_cont('x', cc, k, i1, jp, jp, i2, 'Closed_I')
              call add_cont('x', cc, k, jp, jp, i1, i2, 'Closed_I')

            end if

          end do
        end do
      end do

      end subroutine I_integrals

!==================================================================
      SUBROUTINE F_integrals
!==================================================================
!     Contributions to hd and hdm from F_integrals  (upper matrix)
!
!------------------------------------------------------------------
      IMPLICIT NONE

      ibegin = 1
      iend = intptr(1)
      dsum = 0
      do ip = ibegin, iend
        if (.not. lused(ip)) cycle
        call unpacki(1,ip, kv, i1, i2, i3, i4)
        Print '(A,5I4,F9.5)', 'F: ip, i1, i2, i3, i4', &
                               ip, i1, i2, i3, i4, coef(ip)
        c = coef(ip)

        if (i1 == i2 )then
          call add_cont('d', 2*c, kv, i1, i1, i1, i1, 'F_int') ! hd(i1,i1)
          call add_cont('x', 4*c, kv, i1, i1, i1, i1, 'F_int') !hx(i1,i1)
        else
          call add_cont('d',   c, kv, i2, i2, i1, i1, 'F_int') ! hd(i1,i1)
          call add_cont('d',   c, kv, i1, i1, i2, i2, 'F_int') !hd(i2,i2)
          call add_cont('x', 2*c, kv, i1, i2, i1, i2, 'F_int') !hx(i1,i2)
        end if

      end do
      end subroutine F_integrals

!==================================================================
      SUBROUTINE G_integrals
!==================================================================
!     Contributions to hd and hx from G_integrals  (upper matrix)
!
!------------------------------------------------------------------
      IMPLICIT NONE

      ibegin= intptr(1)+1
      iend  = intptr(2)
      do j = ibegin, iend
        if (.not. lused(j)) cycle
        call unpacki(1,j, kv, i1, i2, i3, i4)
        Print '(A,5I4,F9.5)', 'G: ip, i1, i2, i3, i4', &
                               j, i1, i2, i3, i4, coef(j)
        c = coef(j)

        call add_cont('d', c, kv, i1, i2, i1, i2, 'G-int')  ! hd(i1,i2)
        call add_cont('x', c, kv, i2, i2, i1, i1, 'G-int')  !hx(i1,i1)
        call add_cont('x', c, kv, i2, i1, i1, i2, 'G-int')  !hx(i1,i2)
        call add_cont('x', c, kv, i1, i1, i2, i2, 'G-int')  !hx(i2,i2)

      end do

     end subroutine G_integrals
!
!==================================================================
       SUBROUTINE R_integrals
!==================================================================
!     Contributions to hd and hx from R_integrals  (upper matrix)
!
!------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: itmp, ii,jj
      INTEGER, DIMENSION(4) :: ir
      CHARACTER(LEN=4) :: symm

      ibegin = intptr(4)+1
      iend = intptr(5)
      do j = ibegin, iend
        if (.not. lused(j)) cycle
        call unpacki(5,j, kv, ir(1), ir(2), ir(3), ir(4))
        ! rearrange, if necessary
        if (ir(4) < ir(2)) then
          itmp = ir(4); ir(4)=ir(2); ir(2)=itmp
        end if
        Print '(A,5I4,F9.5)', 'R: jp, i1, i2, i3, i4', &
                               j, ir(1:4), coef(j)

        c = coef(j)

       ! get case selector
       symm = 'abcd'
       do ii = 1,3
          do jj = ii+1,4
             if (ir(jj) == ir(ii)) symm(jj:jj) = symm(ii:ii)
          end do
       end do
       print *, 'Ir:', ir, symm

       select case ( symm )
        case('aaaa')
          call add_cont('d',2*c,kv,ir(1),ir(1), ir(1), ir(1),'R_int')
          call add_cont('x',4*c,kv,ir(1),ir(1), ir(1), ir(1),'R_int')

        case('aaad')
          call add_cont('d',c,  kv,ir(1),ir(4), ir(1), ir(1),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(1), ir(1), ir(4),'R_int')
          call add_cont('s',c,  kv,ir(1),ir(4), ir(1), ir(1),'R_int')
          call add_cont('x',c,  kv,ir(1),ir(1), ir(1), ir(4),'R_int')

        case('aaca')
          call add_cont('d',c,  kv,ir(1),ir(3), ir(1), ir(1),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(1), ir(1), ir(3),'R_int')
          call add_cont('s',c,  kv,ir(1),ir(3), ir(1), ir(1),'R_int')
          call add_cont('x',c,  kv,ir(1),ir(1), ir(1), ir(3),'R_int')

        case('aacc')
          call add_cont('d',c,  kv,ir(1),ir(3), ir(1), ir(3),'R_int')
          call add_cont('x',c,  kv,ir(3),ir(3), ir(1), ir(1),'R_int')
          call add_cont('x',c,  kv,ir(3),ir(1), ir(1), ir(3),'R_int')
          call add_cont('x',c,  kv,ir(1),ir(1), ir(3), ir(3),'R_int')

        case('aacd')
          call add_cont('d',c/2,kv,ir(1),ir(4), ir(1), ir(3),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(3), ir(1), ir(4),'R_int')
          call add_cont('s',c/2,kv,ir(3),ir(4), ir(1), ir(1),'R_int')
          call add_cont('x',c/2,kv,ir(4),ir(1), ir(1), ir(3),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(1), ir(1), ir(4),'R_int')
          call add_cont('x',c/2,kv,ir(1),ir(1), ir(3), ir(4),'R_int')

        case('abab')
          call add_cont('d',  c,kv,ir(2),ir(2), ir(1), ir(1),'R_int')
          call add_cont('d',  c,kv,ir(1),ir(1), ir(2), ir(2),'R_int')
          call add_cont('x',2*c,kv,ir(1),ir(2), ir(1), ir(2),'R_int')

        case('abad')
          call add_cont('d',c,  kv,ir(2),ir(4), ir(1), ir(1),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(1), ir(2), ir(4),'R_int')
          call add_cont('x',c,  kv,ir(1),ir(4), ir(1), ir(2),'R_int')
          call add_cont('x',c,  kv,ir(1),ir(2), ir(1), ir(4),'R_int')

        case('abbb')
          call add_cont('d',c/2,kv,ir(2),ir(2), ir(1), ir(2),'R_int')
          call add_cont('d',c  ,kv,ir(1),ir(2), ir(2), ir(2),'R_int')
          call add_cont('s',c,  kv,ir(1),ir(2), ir(2), ir(2),'R_int')

        case('abbd')
          call add_cont('d',c/2,kv,ir(2),ir(4), ir(1), ir(2),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(2), ir(2), ir(4),'R_int')
          call add_cont('x',c/2,kv,ir(2),ir(4), ir(1), ir(2),'R_int')
          call add_cont('x',c/2,kv,ir(2),ir(2), ir(1), ir(4),'R_int')
          call add_cont('s',c/2,kv,ir(1),ir(4), ir(2), ir(2),'R_int')
          call add_cont('x',c/2,kv,ir(1),ir(2), ir(2), ir(4),'R_int')

        case('abcb')
          call add_cont('d',c/2,kv,ir(2),ir(2), ir(1), ir(3),'R_int')
          call add_cont('d',c  ,kv,ir(1),ir(3), ir(2), ir(2),'R_int')
          call add_cont('x',c,  kv,ir(3),ir(2), ir(1), ir(2),'R_int')
          call add_cont('x',c,  kv,ir(2),ir(1), ir(2), ir(3),'R_int')

        case('abcc')
          call add_cont('d',c/2,kv,ir(2),ir(3), ir(1), ir(3),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(3), ir(2), ir(3),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(3), ir(1), ir(2),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(2), ir(1), ir(3),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(1), ir(2), ir(3),'R_int')
          call add_cont('s',c/2,kv,ir(1),ir(2), ir(3), ir(3),'R_int')

        case('abcd')
          call add_cont('d',c/2,kv,ir(2),ir(4), ir(1), ir(3),'R_int')
          call add_cont('d',c/2,kv,ir(1),ir(3), ir(2), ir(4),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(4), ir(1), ir(2),'R_int')
          call add_cont('x',c/2,kv,ir(3),ir(2), ir(1), ir(4),'R_int')
          call add_cont('x',c/2,kv,ir(1),ir(4), ir(2), ir(3),'R_int')
          call add_cont('x',c/2,kv,ir(1),ir(2), ir(3), ir(4),'R_int')

        case DEFAULT
          Write(err, FMT='(A)')  'Error: symmetry ', symm,' not found'
          stop
        END SELECT

      end do
      end subroutine r_integrals

!==================================================================
      SUBROUTINE hm_residuals
!==================================================================
!     Compute the orbital residuals and lagrange multipliers
!------------------------------------------------------------------
      USE mchf_inout
      IMPLICIT NONE
      INTEGER :: ib, iib, ie, m, i,j, jb, je, irow
      REAL(KIND=8)                     :: eav
      REAL(KIND=8), DIMENSION(ns)      :: z
!
      eij = 0.d0
      ib = 1; ie = ns
      Do i = n1+1,n2
        ! upper fixed
        rhs(:,i) = 0
        do j = 1, n1
          call bxv(ks,ns,hda(:,:,j,i),p(:,j),y)
          rhs(:,i) = rhs(:,i) +y
        end do
        do j = n1+1,nwf
          if (j<i) then
            call bxv(ks,ns,hdb(:,:,j,i),p(:,j),y)
          else
            call bxv(ks,ns,hdb(:,:,i,j),p(:,j),y)
          end if
          rhs(:,i) = rhs(:,i) +y
        end do
        rhs(1,i) = 0; rhs(ns-1:ns,i) =0
        res(ib:ie) = -rhs(:,i)
        Do j = 1,nwf
          if (i /= j .and. (clsd(i) .and. clsd(j))) cycle
          if (l(i) == l(j)) then
              eij(i,j) = dot_product(rhs(:,i), p(:,j))
          end if
        end do
        ib = ib+ns; ie = ie+ns
      end do

      write(50, *) 'Energy matrix before symmetrization: All_NR'
      do i = n1+1, n2
         write(50, '(20F15.8)')  eij(i,:)
      end do

      ! Symmetrize the energy matrix for orbitals that are varied
      Do i = n1+1,n2
         Do j = i+1, n2
            if (clsd(i) .and. clsd(j)) then
                eij(i,j) = 1.d-13
                eij(j,i) = 1.d-13
                cycle
            end if
            If (abs(e(i,j)) > 1.D-12) then
               eav = (eij(i,j) + eij(j,i))/2
               eij(i,j) = eav; eij(j,i) = eav
            end if
          end do
      end do

      ! Generate the full upper hm =hd + hx matrix for orbitals to
      ! be varied
      hm=0
      do i = n1+1,n2
         do j = i,n2
            call add_band_to_full(hdb(:,:,i,j), hm(:,:,i,j))
            hm(:,:,i,j) = hm(:,:,i,j) + hx(:,:,i,j)
         end do
      end do

      ! Modify hm matrix and residuals for lagrange multipliers
      ib = 1; ie = ns;
      do i = n1+1, n2
        jb = 1; je = ns
        Do j = 1,nwf
          if (i == j .or.  abs(e(i,j)) > 1.d-12) then
            if (j <= n1 .or. j >n2 ) then    ! orthogonality to fixed
               call bxv(ks,ns,sb,p(:,j),y)
               y(1)=0; y(ns-1:ns) =0
               res(ib:ie) = res(ib:ie) +eij(i,j)*y
            else if ( j >=i .and. j<=n2 .and. abs(e(i,j)) > 1.d-12 ) then
               hm(:,:,i,j) = hm(:,:,i,j)-eij(i,j)*bb
               !  an off-diagonal matrix element affects two residuals
               call bxv(ks,ns,sb,p(:,j),y)
               y(1)=0; y(ns-1:ns) =0
               res(ib:ie) = res(ib:ie) + eij(i,j)*y
               if ( i /= j) then
                  call bxv(ks,ns,sb,p(:,i),y)
                  y(1)=0; y(ns-1:ns) =0
                  res(jb:je) = res(jb:je) + eij(j,i)*y
               end if
            end if
          end if
          if (i == j) then
            call apply_bc(hm(:,:,i,j), 'n')
          else if (j > i .and. j <=n2) then
            call apply_bc(hm(:,:,i,j), 'o')
          end if
          jb = jb+ns; je = je+ns
        end do
        res(ib) = 0.d0; res(ie-1:ie) = 0.d0
        ib = ib+ns; ie = ie+ns
      end do

      END SUBROUTINE hm_residuals

    !==================================================================
      SUBROUTINE hm_norm
    !==================================================================
    !  Generate the NR equations that include normality constraints
    !------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: ib, iib, ie, jb, je
      ! generate  upper part of aa(1:nns:1:nns) = hm + hx
      ! could also use

      ! aa could be reshaped from the hm has block structure
       ib = 1; ie = ns
       do i  = n1+1, n2
         jb = ib; je = ie
         do j = i, n2
             aa(ib:ie, jb:je) =  hm (:,:,i,j)
             jb = jb+ns; je = je+ns
         end do
         ib = ib+ns;  ie = ie+ns
       end do

      ! Add the normalization constraint for orbital i
      mm = nns; ib=1; ie = ns
      do i = n1+1, n2
        mm = mm+1
        call bxv(ks,ns,sb,p(:,i),y)
        y(1) = 0.d0
        y(ns-1:ns)= 0.d0
        aa(ib:ie,mm) = -y
        aa(mm,ib:ie) =-y
        res(mm) = 0.d0   ! since all current orbitals are normalized
        ib = ib+ns; ie = ie+ns
      end do

      END SUBROUTINE hm_norm

    !==================================================================
       SUBROUTINE hm_orthog
    !==================================================================
    !   Add orthogonality constraint to assure the change is orthogonal
    !   to a given orbital
    !------------------------------------------------------------------

       !Add orthogonality constraints
       !  If (l(i) = l(j))
       ib=1; ie=ns
       do i = n1+1,n2
        jb =  ib; je = ie
        do j = i+1,n2
          jb = jb+ns; je = je+ns
          if (clsd(i) .and. clsd(j)) cycle
          if (l(i) == l(j) .and. abs(e(i,j))> 1.d-12 ) then
              mm = mm+1
              call bxv(ks,ns,sb,p(:,j),y)
              y(1) = 0.d0
              y(ns-1) = 0.d0
              y(ns)= 0.d0
              aa(ib:ie,mm) = -y
              aa(mm,ib:ie) = -y
              call bxv(ks,ns,sb,p(:,i),y)
              y(1) = 0.d0
              y(ns-1:ns) = 0.d0
              aa(jb:je,mm) = -y
              aa(mm,jb:je) = -y
              res(mm) =  sum(y*p(:,j))
           end if
         end do
         ib=ib+ns; ie=ie+ns
       end do

       If ( mm /= nd) then
         Print *, "mm /= nd : mm, nd", mm,nd
         STOP "Something wrong"
       end if
      END SUBROUTINE hm_orthog
   !==================================================================
        SUBROUTINE hm_check
    !==================================================================
    !   Check the aa matrix defining the Newton Raphson equations
    !------------------------------------------------------------------
       integer :: ib, ie, jb, je

        ! Check the hm-blocks
        mm = ns*(n2-n1)
        print *, 'mm, north, nns, nd', mm, north, nns, nd
        ib =1; ie=ns
        do i = n1+1, n2
           jb = ib; je=ie
           do j = i, nwf
              Print '(10X,2I3,F16.9)',i,j, dot_product(p(:,i), &
                  matmul(aa(ib:ie,jb:je),p(:,j)))
              if (i==j ) then
                 mm = mm+1
              end if
              write(52,'(A,6I3)') 'Block:', i,j, ib,ie, jb, je
              do ii = ib, ie
                write(52,'(8F12.6)') aa(ii,jb:je)
              end do
              jb = jb+ns; je=je+ns
           end do
         ! Check the residuals
           Print *, 'a^t,res(a)', dot_product(p(:,i),res(ib:ie))
           Print *, 'Norm(1):', dot_product(p(:,i),aa(ib:ie,mm)), mm
           ib = ib+ns; ie=ie+ns
        end do
        mm = ns*(n2-n1)
        write(52, *) 'Orthonormality'
        Do mm = ns*(n2-n1)+1, ns*(n2-n1) + nwf
          write(52,FMT='(8F12.6)') aa(:, mm)
        end do
        write(52, *) 'Orthogonality'
        Do mm = (ns+1)*(n2-n1)+1, (ns+1)*(n2-n1)+north
          write(52,FMT='(8F12.6)') aa(:, mm)
        end do
        write(52, *) 'Residuals'
        write(52,FMT='(8F12.6)') res
        write(52, *) 'First'
        write(52,FMT='(8F12.6)') aa(1,:)

        END SUBROUTINE hm_check

    !==================================================================
       SUBROUTINE hm_solve
    !==================================================================
    !   Solve and test the solution of the Newton Raphson equations
    !------------------------------------------------------------------
       INTEGER :: ib, ie, jb, je, irow


       CALL DSYTRF('U',nd, aa, nd, ipiv, work_nr, iwork_nr, info)
       If (info /= 0) then
        Print *, ' Error in factorization: DSYTRF in scfall_nr', info
        STOP
       End if

       CALL DSYTRS('U',nd, 1, aa, nd, ipiv, res, nd, info)
       If (info /= 0) then
         Write(err, *) ' Error in solve routine : DSYTRS in scfall_nr', info
        STOP
       End if

       END SUBROUTINE hm_solve

    !==================================================================
       SUBROUTINE hm_update
    !==================================================================
    !  Solve and test the solution of the Newton Raphson equations
    !------------------------------------------------------------------

       ! save the current set of orbitals
       ib = 1; ie = ns
       do i = n1+1, n2
        v = p(:,i)+res(ib:ie)
        do  j = ns,ns/2,-1
          if (abs(v(j)) > end_tol) exit
        end do
        maxr(i) = min(ns,j+1)
        v(j+1:ns) = 0.d0
        e(i,i) = (eij(i,i) +res(nns+i))
        dpm(i) = maxval( abs(p(:,i)-v(:)))/maxval(abs(p(:,i)))
        p(:,i) = v
        call normalize(p(:,i))
        ib = ib+ns; ie = ie+ns
       end do
       END SUBROUTINE hm_update

    !==================================================================
       SUBROUTINE add_cont(type, c, k, i1, i2, i3, i4, name)
    !==================================================================
    !  Add the contribution from an integral to hd or hx
    !------------------------------------------------------------------

       IMPLICIT NONE
       CHARACTER(LEN=1), INTENT(IN) :: type
       REAL(KIND=8), INTENT(IN)     :: c
       INTEGER, INTENT(IN)          :: k
       INTEGER, INTENT(INOUT)       :: i1, i2, i3, i4
       CHARACTER(LEN=*), INTENT(IN) :: name

       INTEGER :: itmp

       if (i3 > i4) then
         itmp = i3; i3=i4; i4=itmp
         itmp = i1; i1=i2; i2=itmp
       end if

       if (i3 > n2 .or. (i3 <= n1 .and. i4 > n2))  return

       ! print '(A,1X,A,F12.9,5I3)','T, c,k,i1,i2,i3,i4', &
       !                              Type, c,k, i1,i2,i3,i4
       if (type == 'd') then
          call density(ns,ks,d,p(:,i1), p(:,i2),'s',ms)
          if( i3 <= n1 .and. i4 <= n2) then
             call add_rkm_d(k,c*d, hda(:,:,i3,i4))
          else if( i3 > n1 .and. i3 <= n2) then
             print *, 'Adding: d',c, i1,i2, 'to', i3,i4
             call add_rkm_d(k, c*d, hdb(:,:,i3,i4))
          else
            write(err,FMT= '(A,F12,8,4I3,A)')  &
                 'Error: not in hd range', c,k,i1,i2,i3,i4, name
            stop
          end if
        else if (type == 'x') then
          if ( (i3 > n1 .and. i3 <= n2) .and. (i4 >n1 .and. i4 <= n2)) then
             call density(ns,ns,dx, p(:,i1), p(:,i2),'x', ms)
             call add_rkm_x(k, c*dx, hx(:,:,i3,i4))
          end if
        else if (type == 's') then
          if ( (i3 > n1 .and. i3 <= n2) .and. (i4 >n1 .and. i4 <= n2)) then
             call density(ns,ns,dx, p(:,i1), p(:,i2),'x', ms)
             dx = dx + transpose(dx)
             call add_rkm_x(k, c*dx, hx(:,:,i3,i4))
          end if
       else
          print *, 'Error: unkown type', type, name
       end if

       END SUBROUTINE add_cont
     END SUBROUTINE some_nr



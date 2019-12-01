!=====================================================================
      subroutine solve_MCHF_equations
!=====================================================================
!
!  Solve the MCHF equations 
!------------------------------------------------------------------
      USE spline_param, ONLY: h,hmax,rmax,ns,ks,ml
      USE spline_galerkin, ONLY: sb,bb
      USE mchf_atomic_state
      USE mchf_inout
      USE orbitals
      USE mchf_param
      USE block_param
      USE CL_arguments

      IMPLICIT NONE
     Interface
       SUBROUTINE scf_nr(i,north, jorth, hfm, bb, hx, rhs, md, v, eii)
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
        REAL(KIND=8), DIMENSION(ns), INTENT(in)    ::  rhs
        REAL(KIND=8), DIMENSION(ns), INTENT(inout) :: v
        REAL(KIND=8), INTENT(out) :: eii
       END SUBROUTINE  scf_nr

       SUBROUTINE scf_matrix(i,hfm, hx, rhs) 
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
       END SUBROUTINE scf_matrix

      End Interface
      
      REAL(KIND=8), DIMENSION(ns,ns) :: hfm, hx
      REAL(KIND=8), DIMENSION(ns,ns) :: aa
      REAL(KIND=8), DIMENSION(ns) :: eigval
      REAL(KIND=8), DIMENSION(ns) ::  v, rhs

      REAL(KIND=8) :: sumi,c,eval, ans, end_diff, e_acc, ratio
      REAL(KIND=8), EXTERNAL :: a,b, bhl, bvmv, fk, gk, azl
      INTEGER :: ios, eof, i,j,m, mm, md, ll, k, it, nw, north, nd
      INTEGER, DIMENSION(nwf) :: jorth

      REAL(KIND=8)  :: et, dp, orb_diff, scf_diff, eii
      INTEGER       :: ipr, ic, lp, ip, max_it1, max_it2, nit
      LOGICAL, SAVE :: last, found, econv, converged, icycle, first=.true.

      Call Allocate_rka(0,kmax)

      ! getting started
      ipr = 0;
      last = .false.

      ! get the estimate of the initial energy
      et = etl_sum
      call update
      call diag(etotal)
      e_acc = abs(et/etotal - 1.d0)

      !  write header information
      write(scr,'(/7X,A,F18.10)') &
                    'Weighted total energy from initial estimates =', etotal
      write(err,'(/7X,A)') 'SCF Phase'
      print '(7X,A,1PD15.7)', 'Energy estimate accuracy:', e_acc

      ! determine the maximum number of iterations
      if ( .not. first) then
          max_it1 = 0
      else if (CL_max_it1 > 1) then
      !if (CL_max_it1 > 1) then
          max_it1 = CL_max_it1
      else if ( e_acc   < 1.d-05 .and. maxval(dpm) < 1.d-05  &
                 .and. acc == 1 ) then
         max_it1 = 2
         converged = .true.
      else if (e_acc < 1.d-3 .and. maxval(dpm) < 1.d-05) then
         max_it1 = 4
         dpm = 1.d-3
      else
         max_it1 = 4
         dpm = 1.d0
      end if
      call scf_phase
       
       nit = nit2
       ! update all orbitals simultaneously using a Newton Raphson method
      if ( nit == nwf) then
         !estimate the number or iterations for the ALL_NR_phase
        if (CL_max_it2 >0) then
           max_it2 = CL_max_it2
        else If (converged) then
           if ( nit < nwf)  then
               max_it2=0   !  this is the Level one
           else
               max_it2=20
           end if
           Print '(/10X,A,I2,A)', 'The first phase of Level',acc, ' has converged'
           Print '(10X,A,I2)', 'max_it2 for the second phase set to', max_it2
        else
           max_it2 = 12
           Print '(/10X,A,I2,A)', 'The first phase of Level',acc, ' has not converged'
        end if
        !max_it2 = 1
        call NR_phase
      end if

      ! finish the current level calculation
      if (acc == acc_max) then
         last = .true.
         call prepv(last)
         call write_bsw
         Call plot_bsw(plt)
      end if
      etl_sum = etotal
      call summry(orb_diff, scf_diff, end_diff)
      Call dealloc_integrals
      first = .false.

  CONTAINS

    !==================================================================
        SUBROUTINE scf_phase
    !==================================================================
    !   Update orbitals one at a time  using either an eigenvalue method
    !   or NR for a single orbital, including lagrange multipliers except
    !   for closed shells.
    !------------------------------------------------------------------

        IMPLICIT NONE
      
      converged = .false.
      nit = nit1
      vary = .true.
      if (nit < nwf)  vary(1:nwf-nit) = .false.
      
      et = etotal
      Do it = 1, max_it1
        write(6,'(/7X,A,I2.2/7X,A/)') 'Iteration 1.',it,'----------------'
        WRITE (6,'(11X,A,11X,A,15X,A,10X,A,8X,A)') &
            'nl', 'E(nl)','AZ', 'DPM', 'MAXR'
        call prepv(last)
        call lagrange
        Do ip = nwf-nit+1, nwf
           i = iord(ip)
           m = n(i)-l(i)
           north = 0
           !.. set orthogonality constraints by forming the list jorth
           do j = 1,nwf

             if (i == j .or. l(i) /= l(j)) cycle
             if (clsd(i) .and. clsd(j)) cycle  ! if both shells are closed

	     if (it == 1 .and. (in(i)==1 .or. in(j)==1) .and. j >=i &
                          .and. acc ==1) then
                cycle    ! orthogonal only to j < i
             else
                north = north+1
                jorth(north) = j
                if (j < i) m = m-1         ! index of eiv reduces only for j<i
             end if

           end do
           
           call scf_matrix(i,hfm, hx, rhs)
           !if (it ==1 .or. abs( qsum(i) -1.d0)  < 0.5 .or. dpm(i) > 0.5) then
           if (it ==1 .or. abs( qsum(i) -1.d0)  < 0.001 .or. dpm(i) > 0.5) then
	     call scf_eiv('V')
           else 
             md = ns+1+north
             v = p(:,i)                         ! needed for hf_nr
             call scf_nr(i,north,jorth,hfm,bb,hx,rhs,md,v,eii)
             e(i,i) = eii
           end if
           dpm(i) = maxval( abs(p(:,i)-v(:)))/maxval(abs(p(:,i)))
           do j = ns, ns/3,-1
              if (abs(v(j)) > end_tol ) exit
           end do
           maxr(i)=min(ns,j+1)
           v(j+1:ns) = 0.d0
           p(:,i) =v(:)             
           lp = l(i)+1
           az(i) = p(lp+1,i)*azl(z,h,ks,lp) 
           WRITE(scr,'(10X,A3,F20.12,F18.10,1P,D12.2,I8)') &
                 el(i), e(i,i), az(i), dpm(i), maxr(i)
         end do
         ! test convergence
         call update
         call diag(etotal)
         orb_diff =  maxval(abs(dpm(1:nwf)))
         scf_diff =  abs(et-etotal)/abs(etotal)
         do i= maxr(nwf),ns/2,-1
           if (abs(p(i,nwf))  > end_tol) exit
         end do
         print '(/(10X,A,T50,1P,D10.2,D12.2))', &
           'SCF convergence (diff vs tol)', scf_diff, scf_tol, &
           'Orbital convergence (diff vs tol)', orb_diff, orb_tol, &
           'Tail cut-off (orbital nwf)', p(i,nwf), end_tol
         write(scr,'(/7X,A,2F18.10)') 'Weighted total energy =', etotal
         write(err,'(7X,A,D12.5,3X,A,F18.10)') 'DeltaE =', etotal-et, &
                       'Weighted total energy = ' , etotal
         If ( orb_diff < orb_tol .or. scf_diff  < 1000*scf_tol) then
            converged = .true.
            exit
         else
             et = etotal
         End if
      End do
      END SUBROUTINE scf_phase

    !==================================================================
     SUBROUTINE NR_phase
    !==================================================================
    !   Update all orbitals simultaneousl using a NR method
    !------------------------------------------------------------------

      IMPLICIT NONE
     
      vary = .true.
      if (nit < nwf)  vary(1:nwf-nit) = .false.

      write(err,'(/7X,A)') 'ALL_NR  Phase'
      ! Compute total number of orthogonality conditions
      north = 0
      Do i = nwf-nit+1,nwf
        Do j = i+1,nwf
           if (clsd(i) .and. clsd(j)) cycle
           if ( l(i) == l(j) .and. abs(e(i,j)) >1.d-12) north = north +1
           !if ( l(i) == l(j) ) north = north +1
        End do
      End do

      Do it = 1,max_it2
         call prepv(last)
         write(6,'(/7X,A,I2.2/7X,A)') 'Iteration 2.',it,'----------------'
         Do ip = nwf-nit+1, nwf
           i = iord(ip)
           Do j = 1,nwf
              if ( j > nwf-nit .and.  i < j .and. abs(e(i,j)) >1.d-12) then
                !call rotate(i,j)
              end if
           End do
         End do

         call some_nr(0, nwf)

         call update
         call diag(etotal)

         ! test convergenc
         orb_diff =  maxval(abs(dpm(1:nwf)))
         scf_diff =  abs(et-etotal)/abs(etotal)
         do i= maxr(nwf),ns/2,-1
           if (abs(p(i,nwf))  > end_tol) exit
         end do
         end_diff = p(i,nwf)
         print '((10X,A,T50,1P,D10.2,D12.2))', &
           'SCF convergence (diff vs tol)', scf_diff, scf_tol, &
           'Orbital convergence (diff vs tol)', orb_diff, orb_tol, &
           'Tail cut-off (orbital nwf)', p(i,nwf), end_tol
         write(scr,'(/7X,A,F18.10)') 'Weighted total energy =', etotal
         write(err,'(7X,A,D12.4,3X,A,F18.10)')'DeltaE =', etotal-et, &
                      'Weighted total energy = ', etotal
         !stop
         If ( orb_diff < orb_tol .or. scf_diff  < scf_tol .or. &
                    etotal-et >1.d-2 ) then
             exit
         else
             et = etotal
         End if
      End do
     
      END SUBROUTINE NR_phase

    !==================================================================
        SUBROUTINE scf_eiv(type)
    !==================================================================
    !    Find the eigenvector of hfm for the m'th eigenvalue after 
    !    orthogonality has been applied. 
    !    For type = 'V', one eigenvector is selected
    !------------------------------------------------------------------

        IMPLICIT NONE
        CHARACTER(LEN=1), INTENT(IN) :: type

        !REAL(KIND=8), DIMENSION(ns,ns) :: aaa
        REAL (KIND=8), DIMENSION(ns,ns) :: ss
        REAL(KIND=8), DIMENSION(ns) ::  y, cv
        REAL(KIND=8), DIMENSION(3*ns) :: W
        REAL(KIND=8) :: c
        INTEGER :: j, jp, INFO
    
      Do jp = 1,north
         j = jorth(jp)
         call apply_orthogonality(hfm,rhs, p(:,j))
      End do

      !.. apply boundary conditions for matrices of size mm
      !   This way hfm will not be destroyed by the eigensolver
      aa = hfm(1:ns,1:ns)
      ss = bb
      CALL apply_bc (aa,'e')
      CALL apply_bc (ss,'n')

      eii  = dot_product(p(:,i), matmul(hfm,p(:,i)))
      print *, 'eii, rhs', eii, dot_product(p(:,i), rhs)

      eii = eii + dot_product(p(:,i), rhs)
      print *, 'eii+ rhs', eii 
        
      ! .. evaluates the eigenvalues and eigenvectors
      call dsygv(1,'V','L',ns,aa,ns,ss,ns,eigval,W,3*ns,INFO)
      ! .. the eigenvectors are now in aa
      if (INFO /= 0) then
        WRITE(UNIT=err,FMT='(A,I6)') 'Error in Eigenvalue routine', info
        STOP
      else
        WRITE(51,*) 'Eigenvalues or orbital matrix i',i
        WRITE(51,'(10F12.8)') eigval(2:ns-2)
      end if

      if (type /= 'V') return

      if (maxval(abs(rhs(2:ns-2))) < orb_tol ) then     
         v(1:ns) = aa(1:ns,1+m)
         if (v(ml) < 0.d0) v = -v
         e(i,i) = eigval(1+m)
      else
         
         v = 0.d0
         do jp = 2,ns-2
           cv(jp) = 0.d0
           ! because the orthogonality transformation moves the eigenvalue to 
           ! zero, it is necessary to skip near zero eigenvectors.
           if (abs(eigval(jp)) > scf_tol ) then
             cv(jp)= dot_product(aa(2:ns-2,jp), rhs(2:ns-2))/(eigval(jp)-eii)
             v = v - cv(jp)*aa(:, jp)
           endif
         end do

         v(1) = 0.d0; v(ns-1:ns) = 0.d0
         e(i,i) = eii
         call normalize(v)
      end if
         
      END SUBROUTINE  scf_eiv

    !==================================================================
      SUBROUTINE  write_bsw
    !==================================================================
    !    Ouput the orbitals
    !------------------------------------------------------------------
      USE spline_galerkin, ONLY: r1
      IMPLICIT NONE
      INTEGER          :: nel, jj, nodes
      REAL(KIND=8)     :: en
      CHARACTER(LEN=4) :: ell = '    '

      Do i = 1,nwf
        Write(ouf) el(i), Z, h, hmax, rmax,ks,ns,maxr(i), e(i,i), dpm(i)
        Write(ouf) p(:,i)
      End do

      ! this needs to be checked  to make sure another diagonalization
      ! does not change the energy!
      if (qsum(nwf) == 1) then
         call scf_matrix(nwf, hfm, hx, rhs)
         north =0
         m = n(nwf)-l(nwf)
         do j = 1,nwf-1
           if (l(nwf) /= l(j)) cycle
           if (clsd(nwf) .and. clsd(j)) cycle  ! if both shells are closed
           north = north+1
           jorth(north) = j
           m = m-1         ! index of eiv reduces only for j<i
         end do
         i = nwf
         call scf_eiv('N')
         nel = n(nwf)
         ell(2:4) = el(nwf)
         Do j = m+1, ns-3
            en = eigval(1+j)/qsum(nwf)
            v(1:ns) = aa(1:ns, 1+j)
            nodes=0
            do jj = 2,ns-3
               if (v(jj)*v(jj-1) < 0) nodes = nodes+1
            end do
            do jj = ns, ns/3,-1
              if (abs(v(jj)) > end_tol ) exit
            end do
            write( ell(1:3), FMT='(I3)') j+ n(nwf) -m
            Print '(7X,A,A5,F13.6,F12.6,I5)', ' Orbital: energy, <r>, nodes', &
                       ell, en, bvmv(ns,ks,r1,'s',v,v), nodes
            Write(ouf) ell, Z, h, hmax, rmax,ks,ns,jj, en
            Write(ouf) v
         end do
      end if

      End Subroutine write_bsw

   END subroutine solve_MCHF_equations

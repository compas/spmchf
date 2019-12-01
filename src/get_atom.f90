!======================================================================
    SUBROUTINE get_atom
!======================================================================
!
!    This routine gets information about the atomic state -- the
!    closed shells and the configuration from which it determines
!    the number of wave functions, and the electrons
!
!----------------------------------------------------------------------
       Use mchf_atomic_state
       Use orbitals
       USE mchf_inout
       USE mchf_param, ONLY: varied1, varied2, vstring
       USE block_param
       USE angular_data

       IMPLICIT NONE
       CHARACTER(LEN=60) :: atomic_state, form, str, who, string
       CHARACTER(LEN=20) :: tmp
       CHARACTER(LEN=3), DIMENSION(25) :: elc
       REAL(KIND=8), DIMENSION(25)    :: qsumc
       REAL(KIND=8)                   :: ss
       INTEGER :: nowf, nb, nch, ipos, len, ipos_B, ipos_E, keigv, &
                  i, j, itemp, ib, ie, ios, ip
       INTEGER, PARAMETER             :: Max_eig = 20
       LOGICAL                        :: connected
       !-----------------------------------------------
       !
       ! Get information about atom and Z
       write (scr, '(A)') ' ATOM Z (separated by blanks)'
       READ (isc, FMT=*) Atom, Z

       !  .. closed shells
       read (iuh, FMT=*) nclosd, str
       read (iuh, FMT='(25(1X,A3))') elc(1:nclosd)

       !  .. other orbitals
       read (iuh, FMT=*) nowf, str
       nwf = nclosd + nowf
       read (iuh, FMT='(18(1X,A3))') elc(nclosd + 1:nwf)

       !  block information
       read (iuh, FMT=*) nblock, idim, ncodim, who
       if (who /= 'snonh') then
          write (err, FMT=*) 'Data not suitable for simultaneous optimization'
          STOP
       end if

       call alloc_block
       !Print *, ' get_atom: Returned from alloc_block'
       ! .. get information for each block
       do nb = 1, nblock
          read (iuh, FMT=*) &
             term_bl(nb), ncfg_bl(nb), nze_bl(nb), cf_tot(nb), nze_max(nb)
          term_bl(nb) = adjustl(term_bl(nb))
       end do

       call get_eigenstates
       ! .. now we know nwf, allocate memory
       call allocate_atomic_state
       call alloc_wt

       el(1:nwf) = elc(1:nwf)

       ! .. store properites of the closed shells
       lmax = 0; SS = 0
       DO i = 1, nwf
          Call el_nl(el(i), n(i), l(i))
          lmax = max(l(i), lmax)
          if (i <= nclosd) then
             qsum(I) = 2*(2*L(I) + 1)
             clsd(i) = .true.
             s(i) = SS + qsum(i)/2
             SS = SS + qsum(i)
             if (ss > int(Z)) then
                write (err, FMT='(7X,A, F4.1,A,F4.1)') &
                   'Too many electrons in the core:   Z=', z, ',  NEL=', ss
                stop
             end if
          else
             clsd(i) = .false.
             s(i) = ss
             qsum(i) = 0.d0
          end if
       END DO
       kmax = 2*lmax

       ! Determine Lagrange multiplier structure
       ! and list of orbitals involved
       DO i = 1, nwf
          e(i, i) = 0.d0
          DO j = 1, i - 1
             e(i, j) = 0.d0
             IF (L(I) .EQ. L(J)) then
                IF (clsd(i) .AND. clsd(j)) then
                   E(I, J) = 1.D-12
                   E(J, I) = 1.d-12
                ElsE
                   E(I, J) = 1.D-5
                   E(J, I) = 1.D-5
                END IF
             END IF
          End do
       End DO

       form = "(17X,'Core =',5(1X,A3,'(',I4,')')/(23X,5(1X,A3,'(',I4,')')))"
       WRITE (UNIT=err, FMT="(//7X,'MCHF WAVE FUNCTIONS FOR  ',A6,' Z =',F5.1,/)") &
          ATOM, Z
       WRITE (UNIT=err, FMT=form) (el(i), int(qsum(i)), i=1, nclosd)
       WRITE (UNIT=err, FMT="(7X,A,10(1X,A3)/(10(1X,I3)/))") &
          'Other Orbitals =', (el(i), i=nclosd + 1, nwf)

       ! .. write information to summary file (unit=3)
       WRITE (UNIT=log, FMT="(//7X,'MCHF WAVE FUNCTIONS FOR  ',A6, &
              ' Z =',F5.1,/)") ATOM, Z
       WRITE (UNIT=log, FMT=form) (el(i), int(qsum(i)), i=1, nclosd)
       WRITE (UNIT=log, FMT="(7X,A,(10(1X,A3)))") 'Other Orbitals =', &
          (el(i), i=nclosd + 1, nwf)

       ! .. order the orbitals by n and l
       iord = 0
       DO i = 1, NWF
          IORD(i) = i
       END DO
       ! We only need shells of the same angular quantum number ordered
       DO i = 1, nwf - 1
          Do j = i + 1, nwf
             IF (N(iord(i)) > N(iord(j))) then
                itemp = iord(i)
                IORD(i) = IORD(j)
                IORD(j) = itemp
             END IF
          End DO
       END DO
       Do ip = 1, nwf
          i = iord(ip)
          print *, i, 'el(i)', n(i), l(i)
       end do
       ! .. determine orbitals to be varied, place at end
       call get_varied

    CONTAINS
!======================================================================
       SUBROUTINE get_eigenstates
!======================================================================
!
!    Determines the eigenvalues for each term (or block)
!
!----------------------------------------------------------------------
          IMPLICIT NONE
          REAL(KIND=8), DIMENSION(Max_eig, nblock) :: eigst_wt
          LOGICAL, DIMENSION(Max_eig, nblock) :: leig
          INTEGER :: niv

          ! .. determine eigenvalues for each block
          write (scr, FMT='(A,I5)') &
             ' The number of blocks (term and parity) is:', nblock
          write (scr, FMT='(A)') &
             ' Enter eigenstates and weights (in parentheses) '

          meig = 0; leig = .false.
          do nb = 1, nblock
             write (scr, '(1X,A)') term_bl(nb)
             read (isc, '(A)') str
             nch = 1
             len = len_trim(str)
             niv = 0
!>>>>
             do while (nch <= len)
                ipos = index(str(nch:len), ',')
                if (ipos .eq. 0) ipos = len + 2 - nch
                read (str(nch:nch + ipos - 2), *) tmp
                ipos_B = index(tmp, '(')

                if (ipos_B .eq. 0) then
                   read (str(nch:nch + ipos - 2), *) keigv
                   eigst_wt(keigv, nb) = 1.0
                else
                   ipos_E = index(tmp, ')')
                   read (str(nch:nch + ipos_B - 2), *) keigv
                   read (tmp(ipos_B + 1:ipos_E - 1), *) eigst_wt(keigv, nb)
                end if

                if (keigv .gt. Max_eig) then
                   write (0, *) 'Too high an eigenvalue requested:', &
                      'Maximum for current dimension is', max_eig
                   stop
                end if

                leig(keigv, nb) = .true.
                nch = nch + ipos
                niv = niv + 1
             end do
!>>>>
             niv_bl(nb) = niv
             nume(nb) = keigv
             meig = max(meig, keigv)
          end do

          call alloc_eigenstates
          ! Set weights, normalized
          eigst_weight = eigst_wt(1:meig, 1:nblock)/sum(eigst_wt(1:meig, 1:nblock))
          leigen = leig(1:meig, 1:nblock)

       END SUBROUTINE get_eigenstates

!======================================================================
       SUBROUTINE get_varied
!======================================================================
!
!    Determines the orbitals to be varied, possibly re-ordering
!
!----------------------------------------------------------------------
          USE cl_arguments
          IMPLICIT NONE
          LOGICAL :: done

          vary(1:nwf) = .false.
          INQUIRE (Unit=iup, OPENED=connected)
          If (.not. connected) then
             varied1 = 'all'; varied2 = 'all'
          END IF

          IF (varied1(1:3) == 'ALL' .OR. varied1(1:3) == 'all') THEN
             NIT1 = NWF
          ELSE IF (varied1(1:4) == 'NONE' .OR. varied1(1:4) == 'none') THEN
             NIT1 = 0
          ELSE IF (INDEX(varied1, '=') /= 0) THEN
             J = INDEX(varied1, '=')
             READ (varied1(J + 1:), *) NIT1
          END IF

          IF (varied2(1:3) == 'ALL' .OR. varied2(1:3) == 'all') THEN
             NIT2 = NWF
          ELSE IF (varied2(1:4) == 'NONE' .OR. varied2(1:4) == 'none') THEN
             NIT2 = 0
          ELSE IF (INDEX(varied2, '=') /= 0) THEN
             J = INDEX(varied2, '=')
             READ (varied2(J + 1:), *) NIT2
          END IF

       End subroutine get_varied

    END SUBROUTINE get_atom
!=======================================================================
    SUBROUTINE eptr(el, elsymb, iel, nwf)
!=======================================================================
!
!   Determines the position of the electron in the electron list
!   Zero if not found.
!
!----------------------------------------------------------------------
!
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: nwf
       CHARACTER(LEN=3), DIMENSION(nwf), INTENT(IN) :: el
       CHARACTER(LEN=3), INTENT(IN) :: elsymb
       INTEGER, INTENT(OUT) :: iel

       INTEGER :: n, i
       iel = 0
       do i = 1, nwf
          if (el(i) .eq. elsymb) then
             iel = i
             return
          endif
       end do
    end

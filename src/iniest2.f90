!***********************************************************************
!
      SUBROUTINE INIEST2(NMAX, NCF, NEIV) 
!
!  Serial version of iniestmpi.
!  Structure of the input sparse matrix hmx:
!    . It's a 1-d array
!    . Length: 1 to jcol(ncf)
!    . Number of non-zero elements for column j is:
!          jcol(j) - jcol(j-1) + 1
!    . Row index for element hmx(i) is irow(i)
!
!     NMAX -- largest matrix size for in memory
!     NCF  -- size of full CI
!     NEIV  -- number of selected eigenvalues
!  Xinghong He  98-10-28
!
!***********************************************************************
      USE dvdson_param
      USE angular_data, irow=> ih 
      USE block_param, jcol=> jptr
      IMPLICIT NONE
      INTEGER , INTENT(IN)       :: NMAX, NCF, NEIV

      INTEGER :: NS,  J, i, ii, NFOUND, INFO, eoff, lwork
      real(KIND=8), dimension(:,:), allocatable :: ap
      real(KIND=8), dimension(:), allocatable   :: eigval,work 

!-----------------------------------------------

      NS = MIN(NMAX,NCF) 
      allocate (ap(ns,ns))
      ap = 0.d0

!  Expand the sparse form to normal form for lower triangle matrix
!  MCHF format

      j = 1
      do ii = 1, jcol(ns)
       if (ii > jcol(j)) j = j+1
       ap(irow(ii),j) = hmx(ii)
      end do
    
      !print *, 'HMX:', hmx
      !print *, 'jcol', jcol
      !print *, 'AP:'
      !do i = 1, ns
      !print '((8F10.5))', ap(:, i)
      !end do
 
!  Merge ap from all nodes and then send to all nodes

      allocate(eigval(ns), work(3*ns))
      lwork = 3*ns

      !CALL DSPEVX ('Vectors also', 'In a range', 'Upper triangular', NS, AP, &
      !   -1., -1., 1, NEIV, 0.D0, NFOUND, EIGVAL, VEC, NS, WORK, IWORK, IFAIL, &
      !    INFO) 
      CALL DSYEV('V', 'L', ns, ap, ns, eigval, work, lwork, info)
     
      if (info /= 0) then
         print *, 'Error from call to DSPEVX in INIEST'
         print '(A/(20I4))', 'Information from DSYEV  INFO', info
         STOP
      end if
      !Print *, 'eigval:', eigval
      !Print *, 'eigvec'
      !do i = 1, ns
      !print '((8F10.5))', ap(:, i)
      !end do

 
!  Build the Basis.
 
      eigvec = 0
 
!  scatter the vectors of length NS to vectors of length NCF
!  Format of basis is that as a linear array as defined in DVDSON
 
      eoff = 0
      DO J = 1, NEIV 
         eigvec(eoff+1:eoff+ns) = ap(:,j)
         eoff = eoff +ncf
      END DO 
      
      eigvec(eoff+1:eoff+neiv) = eigval(1:neiv)
 
      deallocate(ap, eigval, work)
 
      RETURN  
      END SUBROUTINE INIEST2 

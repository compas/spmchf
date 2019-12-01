!-----------------------------------------------------------------------
      SUBROUTINE unpacki(int,i,kv,iel1,iel2,iel3,iel4)
!-----------------------------------------------------------------------
!
!     Given the type of integral, the index of the integral,
!     determine -- kv, iel1, iel2, iel3, iel4
!-----------------------------------------------------------------------
!  
        USE mchf_atomic_state, ONLY: nclosd
        USE angular_data,      ONLY: kval
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: int, i
        INTEGER, INTENT(OUT) :: kv, iel1, iel2, iel3, iel4
!
!      ... unpack the electron data
!      
          !print *, 'unpacki', int, i, kv
          kv = kval(i)
          if (int .le. 3 .or. int .eq. 6) then
              iel2 = mod(kv,64) + nclosd
              kv = kv/64
              iel1 = mod(kv,64) + nclosd
              kv = kv/64
          else if (int .eq. 4 .or. int .eq. 5) then
              iel4 = mod(kv,64) + nclosd
              kv = kv/64
              iel3 = mod(kv,64) + nclosd
              kv = kv/64
              iel2 = mod(kv,64) + nclosd
              kv = kv/64
              iel1 = mod(kv,64) + nclosd
              kv = kv/64
          end if
          end

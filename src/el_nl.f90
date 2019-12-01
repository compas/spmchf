!======================================================================
      SUBROUTINE el_nl(el,n,l) 
!======================================================================
!
!     Determine 'n' and 'l' from electron string
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(LEN=3) , INTENT(INOUT) :: el
      INTEGER, INTENT(OUT) :: n,l
      INTEGER, EXTERNAL :: lval
 
     el = adjustr(el)     ! right justify
     Read (el(1:2), *) n
     l = lval(el(3:3))

     END SUBROUTINE el_nl

!======================================================================
      INTEGER FUNCTION LVAL (SYMBOL) 
!======================================================================
!
!    Look up the l-values associated with the symbol
!----------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER , INTENT(IN) :: SYMBOL 
      INTEGER :: LOCATE 
      CHARACTER(LEN=40) :: & 
                SET='spdfghiklmnopqrstuvwxSPDFGHIKLMNOPQRSTUVWX' 
 
      LOCATE = INDEX(SET,SYMBOL) 
      IF (LOCATE <= 20) THEN 
         LVAL = LOCATE - 1 
      ELSE 
         LVAL = LOCATE - 20 
      ENDIF 
      RETURN  
      END FUNCTION LVAL 


      

!====================================================================
    MODULE mchf_inout
!====================================================================
!
!   contains input/output parameters
!
!--------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INTEGER(4) :: err = 0 !   error (system default for screen)
       INTEGER(4) :: log = 3 !   mchf.log file
       INTEGER(4) :: isc = 5 !   standard input (screen)
       INTEGER(4) :: scr = 6 !   screen output
       INTEGER(4) :: iih = 11 !   ih.nb.lst file
       INTEGER(4) :: iic = 12 !   ico.lst file
       INTEGER(4) :: icl = 13 !   c.lst file
       INTEGER(4) :: iuc = 20 !   input unit for cfg.inp
       INTEGER(4) :: iuf = 21 !   input unit for bsw.inp
       INTEGER(4) :: iup = 22 !   input unit for mchf_param
       INTEGER(4) :: iuy = 23 !   yint.lst file
       INTEGER(4) :: iuh = 24 !   header file cfg.h
       INTEGER(4) :: ouf = 31 !   output unit for bsw.out
       INTEGER(4) :: oul = 32 !   output unit for bsw.l
       INTEGER(4) :: plt = 33 !   plot.dat file (for plotting)

    CONTAINS
!====================================================================
       SUBROUTINE open_files
!====================================================================
!   Open files to be used
!--------------------------------------------------------------------
          IMPLICIT NONE
          Logical :: old

          OPEN (UNIT=log, FILE='spmchf.log', STATUS='UNKNOWN', POSITION='asis')

          INQUIRE (FILE='cfg.h', EXIST=OLD)
          IF (OLD) then
             OPEN (UNIT=IUH, FILE='cfg.h', STATUS='OLD', FORM='FORMATTED')
          Else
             STOP 'OPEN_FILES: cfg.h file not found'
          End if

          INQUIRE (FILE='cfg.inp', EXIST=OLD)
          IF (OLD) then
             OPEN (UNIT=IUC, FILE='cfg.inp', STATUS='OLD', FORM='FORMATTED')
          Else
             STOP 'OPEN_FILES: cfg.inp file not found'
          End if

          INQUIRE (FILE='bsw.inp', EXIST=OLD)
          IF (OLD) &
             OPEN (UNIT=IUF, FILE='bsw.inp', STATUS='OLD', FORM='UNFORMATTED')

          INQUIRE (FILE='mchf_param', EXIST=OLD)
          IF (OLD) OPEN (UNIT=iup, FILE='mchf_param', STATUS='OLD', &
                         FORM='FORMATTED')

          INQUIRE (FILE='yint.lst', EXIST=OLD)
          IF (OLD) then
             OPEN (UNIT=iuy, FILE='yint.lst', STATUS='OLD', FORM='UNFORMATTED')
          Else
             STOP 'OPEN_FILES: yint.lst file not found'
          End if

          INQUIRE (FILE='c.lst', EXIST=OLD)
          IF (OLD) then
             OPEN (UNIT=icl, FILE='c.lst', STATUS='OLD', FORM='UNFORMATTED')
          Else
             STOP 'OPEN_FILES: c.lst file not found'
          End if

          INQUIRE (FILE='ico.lst', EXIST=OLD)
          IF (OLD) then
             OPEN (UNIT=iic, FILE='ico.lst', STATUS='OLD', FORM='UNFORMATTED')
          Else
             STOP 'OPEN_FILES: ico.lst file not found'
          End if

          OPEN (UNIT=OUF, FILE='bsw.out', STATUS='UNKNOWN', FORM='UNFORMATTED')
          OPEN (UNIT=plt, FILE='plot.dat', STATUS='UNKNOWN')
       END subroutine open_files

!====================================================================
    END MODULE mchf_inout

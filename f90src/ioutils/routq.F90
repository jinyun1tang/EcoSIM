SUBROUTINE routq(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE OPENS CHECKPOINT FILES AND READS
!     FILE NAMES FOR PLANT SPECIES AND MANAGEMENT
!
      use data_kind_mod, only : r8 => SHR_KIND_R8
  use fileUtil, only : open_safe
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use PlantMngmtDataType
  use EcoSIMHistMod
  use GridDataType
  use EcoSIMConfig
  implicit none
  integer, intent(in) :: NT,NE,NTX,NEX,NHW,NHE,NVN,NVS


  character(len=*), parameter :: mod_filename = __FILE__
  integer :: NPP(JY,JX)
  CHARACTER(len=16) :: DATAA(JP,JY,JX),DATAB(JP,JY,JX)
  CHARACTER(len=16) :: OUTX,OUTC,OUTM,OUTR,OUTQ
  CHARACTER(len=4) :: CHARY
  CHARACTER(len=2) :: CLIMATE
  integer :: IDATE,IYR,NX,NY,NZ,NN,NH1,NH2,NV1,NV2,NS


! begin_execution

!
! OPEN CHECKPOINT FILES FOR PLANT VARIABLES
!
  IF(is_first_year)THEN
    D9999: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        DO NZ=1,5
          IDAY0(NZ,NY,NX)=-1E+06
          IYR0(NZ,NY,NX)=-1E+06
          IDAYH(NZ,NY,NX)=1E+06
          IYRH(NZ,NY,NX)=1E+06
        enddo
      enddo
    ENDDO D9999
    IF(is_restart_run)THEN
      IDATE=IDATA(9)
    ELSE
      IDATE=IDATA(3)
    ENDIF
!   open checkpoint files for i/o
    WRITE(CHARY,'(I4)')IDATE
    OUTX='P'//DATA1(1)(1:2)//CHARY(1:4)
    OUTC='C'//DATA1(1)(1:2)//CHARY(1:4)
    OUTM='M'//DATA1(1)(1:2)//CHARY(1:4)
    OUTR='R'//DATA1(1)(1:2)//CHARY(1:4)
    OUTQ='Q'//DATA1(1)(1:2)//CHARY(1:4)
    call OPEN_safe(26,outdir,outx,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(27,outdir,outc,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(28,outdir,outm,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(29,outdir,outr,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(30,outdir,outq,'UNKNOWN',mod_filename,__LINE__)
  ENDIF
!
! READ PLANT MANAGEMENT FILE NAMES FOR EACH GRID CELL
!
  if(lverb)write(*,*)'plant management file: ',DATAC(10,NE,NEX)

  IF(DATAC(10,NE,NEX).NE.'NO')THEN
    NN=0
    call OPEN_safe(14,PREFIX,DATAC(10,NE,NEX),'OLD',mod_filename,__LINE__)
50  READ(14,*,END=1000)NH1,NV1,NH2,NV2,NS
!   NN=NN+1
!   maximum number of NS is JP, which currently is set to 5.
    NN=1
    D4995: DO NX=NH1,NH2
      D4990: DO NY=NV1,NV2
        NP0(NY,NX)=NS
        D4985: DO NZ=1,NS
          LSG(NZ,NY,NX)=NN
        ENDDO D4985
      ENDDO D4990
    ENDDO D4995

    IF(NS.GT.0)THEN
!DATAX(NZ): pft file defines pft properties
!DATAY(NZ): management file, seeding information
!IETYP=Koppen climate zone:-2=phytotron
      READ(14,*)(DATAX(NZ),DATAY(NZ),NZ=1,NS)
      D4975: DO NX=NH1,NH2
        D4970: DO NY=NV1,NV2
          D4965: DO NZ=1,NS
            IF(IETYP(NY,NX).GT.0)THEN
              WRITE(CLIMATE,'(I2)')IETYP(NY,NX)
!the type of pft is specified by genra+Koppen climate zone
              DATAX(NZ)=DATAX(NZ)(1:4)//CLIMATE
            ENDIF
          ENDDO D4965
        ENDDO D4970
      ENDDO D4975
    ENDIF

    IF(.not. is_restart_run)THEN
      D8995: DO NX=NH1,NH2
        D8990: DO NY=NV1,NV2
          NP(NY,NX)=NS
          NPP(NY,NX)=0
!DATAP(NZ,NY,NX) and DATAM(NZ,NY,NX) are to be read in readqmod.F90
          D100: DO NZ=1,NP(NY,NX)
            DATAP(NZ,NY,NX)=DATAX(NZ)
            DATAM(NZ,NY,NX)=DATAY(NZ)
          ENDDO D100
          D101: DO NZ=NP(NY,NX)+1,5
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          ENDDO D101
        ENDDO D8990
      ENDDO D8995
    ELSE
!
!     READ PLANT SPECIES NAMES FROM EARLIER RUN IF NEEDED
!
      REWIND(30)
8000  CONTINUE

      D9995: DO NX=NHW,NHE
        D9990: DO NY=NVN,NVS
          READ(30,90,END=1001)IDATE,IYR,NPP(NY,NX) &
            ,(DATAZ(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NPP(NY,NX))
90        FORMAT(2I4,1I3,5(A16,I4))

        ENDDO D9990
      ENDDO D9995


      IF(IDATE.LT.IDAYR.OR.IYR.LT.IYRR)THEN
        GO TO 8000
      ELSEIF(IDATE.GE.IDAYR.AND.IYR.EQ.IYRR)THEN
!
!       MATCH PREVIOUS AND CURRENT PLANT SPECIES
!
        D7975: DO NX=NHW,NHE
          D7970: DO NY=NVN,NVS
            D7965: DO NZ=1,NS
              DATAA(NZ,NY,NX)=DATAX(NZ)
              DATAB(NZ,NY,NX)=DATAY(NZ)
            ENDDO D7965
          ENDDO D7970
        ENDDO D7975

        D7995: DO NX=NH1,NH2
          D7990: DO NY=NV1,NV2
            NP(NY,NX)=MAX(NS,NPP(NY,NX))
            D195: DO NN=1,NP(NY,NX)
              DATAP(NN,NY,NX)='NO'
              DATAM(NN,NY,NX)='NO'
            ENDDO D195
            IF(NPP(NY,NX).GT.0)THEN
              D200: DO NN=1,NPP(NY,NX)
                D205: DO NZ=1,NS
                  IF(DATAZ(NN,NY,NX).EQ.DATAX(NZ).AND.IFLGC(NN,NY,NX).EQ.1)THEN
                    DATAP(NN,NY,NX)=DATAX(NZ)
                    DATAM(NN,NY,NX)=DATAY(NZ)
                    DATAA(NZ,NY,NX)='NO'
                    DATAB(NZ,NY,NX)='NO'
                  ENDIF
                ENDDO D205
              ENDDO D200
!
!             ADD NEW PLANT SPECIES
!
              D250: DO NN=1,NP(NY,NX)
!               WRITE(*,2223)'250',NX,NY,NZ,NN,NP(NY,NX),NS
!                2,(DATAB(NZ,NY,NX),NZ=1,NS),DATAP(NN,NY,NX),DATAM(NN,NY,NX)
                IF(DATAP(NN,NY,NX).EQ.'NO')THEN
                  D255: DO NZ=1,NS
!                   WRITE(*,2223)'255',NX,NY,NZ,NN,NP(NY,NX),NS
!                    2,DATAM(NN,NY,NX),DATAB(NZ,NY,NX)
!2223                FORMAT(A8,6I4,10A16)
                    IF(DATAA(NZ,NY,NX).NE.'NO')THEN
                      DATAP(NN,NY,NX)=DATAA(NZ,NY,NX)
                      DATAM(NN,NY,NX)=DATAB(NZ,NY,NX)
                      DATAA(NZ,NY,NX)='NO'
                      DATAB(NZ,NY,NX)='NO'
                      EXIT
                    ENDIF
                  ENDDO D255
                ENDIF
              ENDDO D250
              D201: DO NZ=NP(NY,NX)+1,5
                DATAP(NZ,NY,NX)='NO'
                DATAM(NZ,NY,NX)='NO'
              ENDDO D201
            ELSE
              D265: DO NZ=1,NS
                DATAP(NZ,NY,NX)=DATAX(NZ)
                DATAM(NZ,NY,NX)=DATAY(NZ)
              ENDDO D265
              D270: DO NZ=NS+1,5
                DATAP(NZ,NY,NX)='NO'
                DATAM(NZ,NY,NX)='NO'
              ENDDO D270
            ENDIF
!
!           SET NUMBER OF PLANT SPECIES
!
            NN=5
            D202: DO NZ=JP,1,-1
              IF(DATAP(NZ,NY,NX).EQ.'NO')THEN
                NN=NN-1
              ELSE
                EXIT
              ENDIF
            ENDDO D202
            NP(NY,NX)=NN
            NP0(NY,NX)=NN
          ENDDO D7990
        ENDDO D7995
      ENDIF
    ENDIF
    GO TO 50
1001 CONTINUE
    D6995: DO NX=NHW,NHE
      D6990: DO NY=NVN,NVS
        NP(NY,NX)=NS
        D300: DO NZ=1,NP(NY,NX)
          DATAP(NZ,NY,NX)=DATAX(NZ)
          DATAM(NZ,NY,NX)=DATAY(NZ)
        ENDDO D300
        D301: DO NZ=NP(NY,NX)+1,5
          DATAP(NZ,NY,NX)='NO'
          DATAM(NZ,NY,NX)='NO'
        ENDDO D301
      ENDDO D6990
    ENDDO D6995
    GO TO 50
1000  CLOSE(14)
  ELSE
    D5995: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NP(NY,NX)=0
        DO NZ=1,NP0(NY,NX)
          DATAP(NZ,NY,NX)='NO'
          DATAM(NZ,NY,NX)='NO'
        enddo
      enddo
    ENDDO D5995
  ENDIF
  RETURN
END subroutine routq

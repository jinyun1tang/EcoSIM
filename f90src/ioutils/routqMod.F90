module routqMod
!
! DESCRIPTION
! code to read plant management info.
  use GridConsts
  use fileUtil, only : open_safe
  use FlagDataType
  use EcoSIMCtrlDataType
  use PlantMngmtDataType
  use EcoSIMHistMod
  use GridDataType
  use EcoSIMConfig

  implicit none

  save
  character(len=*),private, parameter :: mod_filename = __FILE__

  public :: routq
  contains

  SUBROUTINE routq(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE OPENS CHECKPOINT FILES AND READS
!     FILE NAMES FOR PLANT SPECIES AND MANAGEMENT
!
      use data_kind_mod, only : r8 => SHR_KIND_R8
  implicit none
  integer, intent(in) :: NT,NE,NTX,NEX,NHW,NHE,NVN,NVS

  CHARACTER(len=16) :: OUTX,OUTC,OUTM,OUTR,OUTQ
  CHARACTER(len=4) :: CHARY
  integer :: IDATE,NY,NX,NZ


! begin_execution

!
! OPEN CHECKPOINT FILES FOR PLANT VARIABLES
!
  IF(is_first_year)THEN
    D9999: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        DO NZ=1,JP
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
    print*,'file 30: ',outq
    call OPEN_safe(26,outdir,outx,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(27,outdir,outc,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(28,outdir,outm,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(29,outdir,outr,'UNKNOWN',mod_filename,__LINE__)
    call OPEN_safe(30,outdir,outq,'UNKNOWN',mod_filename,__LINE__)
  ENDIF

  call ReadPlantInfo_ascii(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)

  END subroutine routq

!------------------------------------------------------------------------------------------

  subroutine ReadPlantInfo_ascii(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)

  USE EcoSIMHistMod
  implicit none
  integer, intent(in) :: NT,NE,NTX,NEX,NHW,NHE,NVN,NVS
  integer :: IDATE
  integer :: NPP(JY,JX)
  integer :: IYR,NX,NY,NZ,NN,NH1,NH2,NV1,NV2,NS
  CHARACTER(len=2) :: CLIMATE
  CHARACTER(len=16) :: DATAA(JP,JY,JX),DATAB(JP,JY,JX)
!
! READ PLANT MANAGEMENT FILE NAMES FOR EACH GRID CELL
!
!  if(lverb)
  write(*,*)'plant management file: ',DATAC(10,NE,NEX)

  NPP=0
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
!the following code assumes there is only one grid specified by
!the coordinates NH1,NV1,NH2,NV2

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
      print*,'redxx:  ', (DATAX(NZ),DATAY(NZ),NZ=1,NS)
    ENDIF

    print*,'is_restart_run',is_restart_run
    IF(.not. is_restart_run)THEN
! there was no chechk point file read in, so update pft info
! from input file
      D8995: DO NX=NH1,NH2
        D8990: DO NY=NV1,NV2
          NP(NY,NX)=NS
          NPP(NY,NX)=0
!DATAP(NZ,NY,NX) and DATAM(NZ,NY,NX) are to be read in readqmod.F90
          D100: DO NZ=1,NP(NY,NX)
            DATAP(NZ,NY,NX)=DATAX(NZ)
            DATAM(NZ,NY,NX)=DATAY(NZ)
          ENDDO D100

          D101: DO NZ=NP(NY,NX)+1,JP
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          ENDDO D101
        ENDDO D8990
      ENDDO D8995
    ELSE
! is a continued run from checkpoint file
!
!     READ PLANT SPECIES NAMES FROM EARLIER RUN IF NEEDED
!
      print*,'rewind'
      REWIND(30)
8000  CONTINUE

      print*,'here 8000 cont'
!  recover DATAZ from checkpoint file
      D9995: DO NX=NHW,NHE
        D9990: DO NY=NVN,NVS
          READ(30,90,END=1001)IDATE,IYR,NPP(NY,NX) &
            ,(DATAZ(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NPP(NY,NX))
90        FORMAT(2I4,1I3,5(A16,I4))
          print*,IDATE,IYR,NPP(NY,NX),(DATAZ(NZ,NY,NX),IFLGC(NZ,NY,NX),NZ=1,NPP(NY,NX))
        ENDDO D9990
      ENDDO D9995

      print*,IDAYR,IYRR
      IF(IDATE.LT.IDAYR.OR.IYR.LT.IYRR)THEN
        GO TO 8000
      ELSEIF(IDATE.GE.IDAYR.AND.IYR.EQ.IYRR)THEN
!
!       MATCH PREVIOUS AND CURRENT PLANT SPECIES
!
        print*,'match pft',NS,NHW,NHE,NVN,NVS
        D7975: DO NX=NHW,NHE
          D7970: DO NY=NVN,NVS
            D7965: DO NZ=1,NS
              print*,'mat ',DATAX(NZ),DATAY(NZ)
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

              D201: DO NZ=NP(NY,NX)+1,JP
                DATAP(NZ,NY,NX)='NO'
                DATAM(NZ,NY,NX)='NO'
              ENDDO D201
            ELSE

              D265: DO NZ=1,NS
                DATAP(NZ,NY,NX)=DATAX(NZ)
                DATAM(NZ,NY,NX)=DATAY(NZ)
              ENDDO D265
              D270: DO NZ=NS+1,JP
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
    print*,'goto 50:1'
    print*,''
    GO TO 50
1001 CONTINUE

    print*,'NxxxS',NS
    D6995: DO NX=NHW,NHE
      D6990: DO NY=NVN,NVS
        NP(NY,NX)=NS
        D300: DO NZ=1,NP(NY,NX)
          DATAP(NZ,NY,NX)=DATAX(NZ)
          DATAM(NZ,NY,NX)=DATAY(NZ)
        ENDDO D300
        D301: DO NZ=NP(NY,NX)+1,JP
          DATAP(NZ,NY,NX)='NO'
          DATAM(NZ,NY,NX)='NO'
        ENDDO D301
      ENDDO D6990
    ENDDO D6995
    print*,'goto 50:2'
    GO TO 50
1000  CLOSE(14)
    print*,'close14 ',IGO
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
  end subroutine ReadPlantInfo_ascii
!------------------------------------------------------------------------------------------


  subroutine ReadPlantInfoNC(NT,NE,NTX,NEX,NHW,NHE,NVN,NVS)

  USE EcoSIMCtrlMod, only : pft_mgmt_in
  use netcdf
  use ncdio_pio
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: NT,NE,NTX,NEX,NHW,NHE,NVN,NVS
  integer :: IDATE
  integer :: NPP(JY,JX)
  integer :: IYR,NX,NY,NZ,NN,NH1,NH2,NV1,NV2,NS
  CHARACTER(len=2) :: CLIMATE
  CHARACTER(len=16) :: DATAA(JP,JY,JX),DATAB(JP,JY,JX)
  integer :: ntopou,NTOPO
  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear
  integer :: ret

  if (len_trim(pft_mgmt_in)==0)then
    D5995: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NP(NY,NX)=0
        DO NZ=1,NP0(NY,NX)
          DATAP(NZ,NY,NX)='NO'
          DATAM(NZ,NY,NX)='NO'
        enddo
      enddo
    ENDDO D5995
  else
    call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

    call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)
    if(readvar)then
      call endrun('fail to find pft_dflag in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
      'in '//trim(mod_filename))
!    call ncd_getvint(pftinfo_nfid,'pft_dflag',pft_dflag)

    if (pft_dflag==-1)then
!not pft data
      DO NX=NHW,NHE
        DO  NY=NVN,NVS
          NP(NY,NX)=0
          DO NZ=1,NP0(NY,NX)
            DATAP(NZ,NY,NX)='NO'
            DATAM(NZ,NY,NX)='NO'
          enddo
        enddo
      ENDDO
    else
      if(pft_dflag==0)then
        iyear=1
      else
        iyear=1+IGO
      endif
      !obtain the number of topographic units
      ntopou=get_dim_len(pftinfo_nfid, 'ntopou')

      DO NTOPO=1,ntopou
        call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)

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
  ! a vegetated topo unit
  !  input data is of shape 'pft_type(year, ntopou, maxpfts, nchar1)'

        ENDIF

      ENDDO
    ENDIF
    call ncd_pio_closefile(pftinfo_nfid)
  endif
  end subroutine ReadPlantInfoNC

end module routqMod

module ReadManagementMod

  use data_kind_mod, only : r8 => SHR_KIND_R8
  use abortutils , only : endrun
  use fileUtil   , only : open_safe
  use minimathmod, only : isLeap
  use GridConsts
  use FlagDataType
  use FertilizerDataType
  use ClimForcDataType
  use SoilWaterDataType
  use LandSurfDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
implicit none
  private
  character(len=*), parameter :: mod_filename = __FILE__

  public :: ReadManagementFiles
  contains

!------------------------------------------------------------------------------------------
  subroutine ReadTillageFile(FileTillage,NH1,NH2,NV1,NV2)
  implicit none
  character(len=*), intent(in) :: FileTillage
  integer, intent(in) :: NH1,NH2,NV1,NV2

  integer :: NY,NX
  integer  :: IPLOW,IDY
  real(r8) :: DPLOW,DY
  integer  :: LPY,IDY1,IDY2,IDY3
  call OPEN_safe(10,PREFIX,DATA1(8),'OLD',mod_filename,__LINE__)

!
!     DY=date DDMMYYYY
!     IPLOW,DPLOW=intensity,depth of disturbance
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     DCORP=intensity (fire) or depth (tillage,drainage) of disturbance
!
  do while(.TRUE.)
    READ(10,*,END=305)DY,IPLOW,DPLOW
    LPY=0
    IDY1=INT(DY/1.0E+06)
    IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
    IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    DO 8995 NX=NH1,NH2
      DO 8990 NY=NV1,NV2
        ITILL(IDY,NY,NX)=IPLOW
        DCORP(IDY,NY,NX)=DPLOW
8990  CONTINUE
8995  CONTINUE
  enddo
305   CONTINUE
  CLOSE(10)
  end subroutine ReadTillageFile
!------------------------------------------------------------------------------------------

  subroutine ReadIrrigationFile(FileIrrig,NH1,NH2,NV1,NV2)

  implicit none
  character(len=*), intent(in) :: FileIrrig
  integer, intent(in) :: NH1,NH2,NV1,NV2
  integer :: NY,NX,LPY,J,JEN
  integer :: IDY1,IDY2,IDY3,I,IDY
  integer :: IDYS,IHRS,IFLGVX,JST,IDYE,IHRE
  real(r8):: DY
  real(r8) :: DST,DEN,CIRRX,RR,FIRRX
  real(r8) :: DIRRX,PHQX,CKAQX,RRH,WDPTHI,CCLQX
  real(r8) :: CMGQX,CNAQX,CSOQX
  real(r8) :: CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX
  call OPEN_safe(2,PREFIX,FileIrrig,'OLD',mod_filename,__LINE__)

  IF(FileIrrig(1:4).EQ.'auto')THEN
!
!       AUTOMATED IRRIGATION
!
!       DST,DEN=start,end dates,hours DDMMHHHH
!       IFLGVX=flag for irrigation criterion,0=SWC,1=canopy water potl
!       FIRRX=depletion of SWC from CIRRX to WP(IFLGV=0),or minimum canopy
!       water potential(IFLGV=1), to trigger irrigation
!       CIRRX= fraction of FC to which irrigation will raise SWC
!       DIRRX= depth to which water depletion and rewatering is calculated
!       WDPTHI=depth at which irrigation is applied
!       PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,
!       CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
!       concentration in irrigation water
!
    READ(2,*,END=105)DST,DEN,IFLGVX,FIRRX,CIRRX,DIRRX,WDPTHI &
      ,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX &
      ,CSOQX,CCLQX
    LPY=0
!   idy1: day, idy2:mon, idy3:hour
    IDY1=INT(DST/1.0E+06)
    IDY2=INT(DST/1.0E+04-IDY1*1.0E+02)
    IDY3=INT(DST-(IDY1*1.0E+06+IDY2*1.0E+04))
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDYS=IDY1
    else
      IDYS=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    IHRS=IDY3
    LPY=0
    IDY1=INT(DEN/1.0E+06)
    IDY2=INT(DEN/1.0E+04-IDY1*1.0E+02)
    IDY3=INT(DEN-(IDY1*1.0E+06+IDY2*1.0E+04))
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDYE=IDY1
    else
      IDYE=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    IHRE=IDY3
!
!   TRANSFER INPUTS TO MODEL ARRAYS
!
    DO 7965 NX=NH1,NH2
      DO 7960 NY=NV1,NV2
        IFLGV(NY,NX)=IFLGVX
        IIRRA(1,NY,NX)=IDYS
        IIRRA(2,NY,NX)=IDYE
        IIRRA(3,NY,NX)=INT(IHRS/100)
        IIRRA(4,NY,NX)=INT(IHRE/100)
        FIRRA(NY,NX)=FIRRX
        CIRRA(NY,NX)=CIRRX
        DIRRA(1,NY,NX)=DIRRX
        DIRRA(2,NY,NX)=WDPTHI
        DO 220 I=1,366
          PHQ(IDY,NY,NX)=PHQX
          CN4Q(IDY,NY,NX)=CN4QX/14.0
          CNOQ(IDY,NY,NX)=CNOQX/14.0
          CPOQ(IDY,NY,NX)=CPOQX/31.0
          CALQ(IDY,NY,NX)=CALQX/27.0
          CFEQ(IDY,NY,NX)=CFEQX/55.8
          CCAQ(IDY,NY,NX)=CCAQX/40.0
          CMGQ(IDY,NY,NX)=CMGQX/24.3
          CNAQ(IDY,NY,NX)=CNAQX/23.0
          CKAQ(IDY,NY,NX)=CKAQX/39.1
          CSOQ(IDY,NY,NX)=CSOQX/32.0
          CCLQ(IDY,NY,NX)=CCLQX/35.5
220     CONTINUE
7960  CONTINUE
7965  CONTINUE
  ELSE
!
!       SCHEDULED IRRIGATION
!
    do while(.TRUE.)
!
!     DY,RR,JST,JEN=date DDMMYYYY,amount (mm),start and end hours
!     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,
!     CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
!     concentration in irrigation water
!
      READ(2,*,END=105)DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX &
          ,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      LPY=0
      IDY1=INT(DY/1.0E+06)
      IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
      IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
      IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
      IF(IDY2.EQ.1)then
        IDY=IDY1
      else
        IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      endif
      RRH=RR/(JEN-(JST-1))
      DO 8965 NX=NH1,NH2
        DO 8960 NY=NV1,NV2
          DO 2535 J=1,24
            IF(J.GE.JST.AND.J.LE.JEN)RRIG(J,IDY,NY,NX)=RRH/1000.0
2535      CONTINUE
!
!           TRANSFER INPUTS TO MODEL ARRAYS
!
          PHQ(IDY,NY,NX)=PHQX
          CN4Q(IDY,NY,NX)=CN4QX/14.0
          CNOQ(IDY,NY,NX)=CNOQX/14.0
          CPOQ(IDY,NY,NX)=CPOQX/31.0
          CALQ(IDY,NY,NX)=CALQX/27.0
          CFEQ(IDY,NY,NX)=CFEQX/55.8
          CCAQ(IDY,NY,NX)=CCAQX/40.0
          CMGQ(IDY,NY,NX)=CMGQX/24.3
          CNAQ(IDY,NY,NX)=CNAQX/23.0
          CKAQ(IDY,NY,NX)=CKAQX/39.1
          CSOQ(IDY,NY,NX)=CSOQX/32.0
          CCLQ(IDY,NY,NX)=CCLQX/35.5
          WDPTH(IDY,NY,NX)=WDPTHI
8960    CONTINUE
8965  CONTINUE
  enddo
  ENDIF
105 CONTINUE
  CLOSE(2)
  end subroutine ReadIrrigationFile

!------------------------------------------------------------------------------------------

  subroutine ReadFertlizerFile(FertFile,NH1,NH2,NV1,NV2)

  implicit none
  character(len=*), intent(in) :: FertFile
  integer,intent(in) :: NH1,NH2,NV1,NV2
  integer :: NY,NX,LPY
  integer :: IDY1,IDY2,IDY3,IDY
  real(r8) :: FDPTHI,DY
  real(r8) :: Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB,Z4A
  real(r8) :: PMA,PMB,PHA,CAS,CAC
  real(r8) :: RSC1,RSN1,RSP1,RSC2,RSN2,RSP2
  real(r8) :: ROWX
  integer  :: IR0,IR1,IR2
!
!     DY=date DDMMYYYY
!     *A,*B=broadcast,banded fertilizer application
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!     PM*,PH*=Ca(H2PO4)2,apatite
!     CAC,CAS=CaCO3,CaSO4
!     *1,*2=litter,manure amendments
!     RSC,RSN,RSC=amendment C,N,P content
!     FDPTHI=application depth
!     ROWX=band row width
!     IRO,IR1,IR2=fertilizer,litter,manure type
!
  call OPEN_safe(8,PREFIX,FertFile,'OLD',mod_filename,__LINE__)
  do while(.TRUE.)
    READ(8,*,END=85)DY,Z4A,Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB &
      ,PMA,PMB,PHA,CAC,CAS,RSC1,RSN1,RSP1,RSC2,RSN2,RSP2,FDPTHI &
      ,ROWX,IR0,IR1,IR2
    LPY=0
    IDY1=INT(DY/1.0E+06)
    IDY2=INT(DY/1.0E+04-IDY1*1.0E+02)
    IDY3=INT(DY-(IDY1*1.0E+06+IDY2*1.0E+04))
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif

    DO 8985 NX=NH1,NH2
      DO 8980 NY=NV1,NV2
!
!         ENTER AMENDMENTS INTO MODEL ARRAYS
!
!         NH4,NH3,UREA,NO3 BROADCAST (A) AND BANDED (B)
!
        FERT(1,IDY,NY,NX)=Z4A
        FERT(2,IDY,NY,NX)=Z3A
        FERT(3,IDY,NY,NX)=ZUA
        FERT(4,IDY,NY,NX)=ZOA
        FERT(5,IDY,NY,NX)=Z4B
        FERT(6,IDY,NY,NX)=Z3B
        FERT(7,IDY,NY,NX)=ZUB
        FERT(8,IDY,NY,NX)=ZOB
!
!         MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE BROADCAST (A)
!         AND BANDED (B)
!
        FERT(9,IDY,NY,NX)=PMA
        FERT(10,IDY,NY,NX)=PMB
        FERT(11,IDY,NY,NX)=PHA
!
!         LIME AND GYPSUM
!
        FERT(12,IDY,NY,NX)=CAC
        FERT(13,IDY,NY,NX)=CAS
!
!         PLANT AND ANIMAL RESIDUE C, N AND P
!
        FERT(14,IDY,NY,NX)=RSC1
        FERT(15,IDY,NY,NX)=RSN1
        FERT(16,IDY,NY,NX)=RSP1
        FERT(17,IDY,NY,NX)=RSC2
        FERT(18,IDY,NY,NX)=RSN2
        FERT(19,IDY,NY,NX)=RSP2
!
!         DEPTH AND WIDTH OF APPLICATION
!
        FDPTH(IDY,NY,NX)=FDPTHI
        ROWI(IDY,NY,NX)=ROWX
!
!         TYPE OF FERTILIZER,PLANT OR ANIMAL RESIDUE
!
        IYTYP(0,IDY,NY,NX)=IR0
        IYTYP(1,IDY,NY,NX)=IR1
        IYTYP(2,IDY,NY,NX)=IR2
8980    CONTINUE
8985  CONTINUE
  enddo
85 CONTINUE
  CLOSE(8)
  end subroutine ReadFertlizerFile

!------------------------------------------------------------------------------------------

  subroutine ReadManagementFiles(FileManage)

  implicit none
  character(len=*), intent(in) :: FileManage

  integer :: NH1,NV1,NH2,NV2
!
!   NH1,NV1,NH2,NV2=N,W and S,E corners of landscape unit
!   DATA1(8),DATA1(5),DATA1(6)=disturbance,fertilizer,irrigation files
!   PREFIX=path for files in current or higher level directory

  call OPEN_safe(13,PREFIX,FileManage,'OLD',mod_filename,__LINE__)
  do while(.TRUE.)
    READ(13,*,END=200)NH1,NV1,NH2,NV2
    READ(13,*)DATA1(8),DATA1(5),DATA1(6)
!
!   READ TILLAGE INPUT FILE
!
    IF(DATA1(8).NE.'NO')THEN
      call ReadTillageFile(DATA1(8),NH1,NH2,NV1,NV2)
    ENDIF
!
!   READ FERTLIZER INPUT FILE
!
    IF(DATA1(5).NE.'NO')THEN
      call ReadFertlizerFile(DATA1(5),NH1,NH2,NV1,NV2)
    ENDIF
!
!   READ IRRIGATION INPUT FILE
!
    IF(DATA1(6).NE.'NO')THEN
      call ReadIrrigationFile(DATA1(6),NH1,NH2,NV1,NV2)
    ENDIF
  enddo
200 CONTINUE
  CLOSE(13)
  end subroutine ReadManagementFiles
end module ReadManagementMod

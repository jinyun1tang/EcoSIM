module ReadManagementMod
  !!
  ! Description
  ! code to read Soil management info
  !
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils , only : endrun
  use fileUtil   , only : open_safe,int2str
  use minimathmod, only : isLeap
  use GridConsts
  use FlagDataType
  use FertilizerDataType
  use ClimForcDataType
  use EcoSIMCtrlMod, only : lverb, first_topou
  use SoilWaterDataType
  use LandSurfDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
  USE ncdio_pio
  use netcdf
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ReadManagementFiles
  contains

!------------------------------------------------------------------------------------------
  subroutine ReadTillageFile(soilmgmt_nfid,FileTillage,NH1,NH2,NV1,NV2)
  implicit none
  character(len=*), intent(in) :: FileTillage
  integer, intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid

  integer :: NY,NX
  integer  :: IPLOW,IDY
  real(r8) :: DPLOW,DY
  integer  :: LPY,IDY1,IDY2,IDY3
  integer :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=24) :: tillf(12)
  integer :: ntill

  ntill=get_dim_len(soilmgmt_nfid,'ntill')
  if(ntill>12)then
    call endrun('Not enough memory size for array tillf in '//trim(mod_filename),__LINE__)
  endif
  call check_var(soilmgmt_nfid, FileTillage, vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find '//trim(FileTillage)//' in '//trim(mod_filename), __LINE__)
  endif

  call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, tillf),&
      trim(mod_filename))

!
!     DY=date DDMMYYYY
!     IPLOW,DPLOW=intensity,depth of disturbance
!     ITILL=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
!     DCORP=intensity (fire) or depth (tillage,drainage) of disturbance
!
  kk=1
  do while(len_trim(tillf(kk))>0)

    read(tillf(kk),'(I2,I2,I4)')IDY1,IDY2,IDY3
    LPY=0
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    read(tillf(kk),*)DY,IPLOW,DPLOW
    if(lverb)then
      print*,tillf(kk)
      print*,idy1,idy2,idy3,IPLOW,DPLOW
    endif
    D8995: DO NX=NH1,NH2
      D8990: DO NY=NV1,NV2
        iSoilDisturbType_col(IDY,NY,NX) = IPLOW
        DepzCorp_col(IDY,NY,NX)                = DPLOW
      ENDDO D8990
    ENDDO D8995
    kk=kk+1
  enddo
  end subroutine ReadTillageFile
!------------------------------------------------------------------------------------------

  subroutine ReadIrrigationFile(soilmgmt_nfid,FileIrrig,NH1,NH2,NV1,NV2)

  implicit none
  character(len=*), intent(in) :: FileIrrig
  integer, intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid

  integer :: NY,NX,LPY,J,JEN
  integer :: IDY1,IDY2,IDY3,I,IDY
  integer :: IDYS,IHRS,IFLGVX,JST,IDYE,IHRE
  real(r8):: DY
  real(r8) :: DST,DEN,CIRRX,RR,FIRRX
  real(r8) :: DIRRX,PHQX,CKAQX,RRH,WDPTHI,CCLQX
  real(r8) :: CMGQX,CNAQX,CSOQX
  real(r8) :: CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX
  integer :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=128) :: irrigf(24)
  integer :: nirrig

  nirrig=get_dim_len(soilmgmt_nfid,'nirrig')
  if(nirrig>24)then
    call endrun("not enough size for array irrigf in "//trim(mod_filename), __LINE__)
  endif

  call check_var(soilmgmt_nfid, FileIrrig, vardesc, readvar)

  if(.not. readvar)then
    call endrun('fail to find '//trim(FileIrrig)//' in '//trim(mod_filename), __LINE__)
  endif

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

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf(1)),&
      trim(mod_filename))

    READ(irrigf(1),*)DST,DEN,IFLGVX,FIRRX,CIRRX,DIRRX,WDPTHI &
      ,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX &
      ,CSOQX,CCLQX
    READ(irrigf(1),'(I2,I2,I4)')IDY1,IDY2,IDY3

    IF(lverb)then
      print*,irrigf(1)
      print*,IDY1,IDY2,IDY3,DEN,IFLGVX,FIRRX,CIRRX,DIRRX,WDPTHI &
        ,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX &
        ,CSOQX,CCLQX
    endif

    LPY=0
!   idy1: day, idy2:mon, idy3:hour
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
    D7965: DO NX=NH1,NH2
      D7960: DO NY=NV1,NV2
        IFLGV(NY,NX)   = IFLGVX
        IIRRA(1,NY,NX) = IDYS
        IIRRA(2,NY,NX) = IDYE
        IIRRA(3,NY,NX) = INT(IHRS/100)
        IIRRA(4,NY,NX) = INT(IHRE/100)
        FIRRA(NY,NX)   = FIRRX
        CIRRA(NY,NX)   = CIRRX
        DIRRA(1,NY,NX) = DIRRX
        DIRRA(2,NY,NX) = WDPTHI
        D220: DO I     = 1, 366
          PHQ(IDY,NY,NX)                                    = PHQX
          NH4_irrig_mole_conc(IDY,NY,NX)                    = CN4QX/14.0_r8
          NO3_irrig_mole_conc(IDY,NY,NX)                    = CNOQX/14.0_r8
          H2PO4_irrig_mole_conc(IDY,NY,NX)                  = CPOQX/31.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Al,IDY,NY,NX)  = CALQX/27.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Fe,IDY,NY,NX)  = CFEQX/55.8_r8
          trcsalt_irrig_mole_conc_col(idsalt_Ca,IDY,NY,NX)  = CCAQX/40.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Mg,IDY,NY,NX)  = CMGQX/24.3_r8
          trcsalt_irrig_mole_conc_col(idsalt_Na,IDY,NY,NX)  = CNAQX/23.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_K,IDY,NY,NX)   = CKAQX/39.1_r8
          trcsalt_irrig_mole_conc_col(idsalt_SO4,IDY,NY,NX) = CSOQX/32.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Cl,IDY,NY,NX)  = CCLQX/35.5_r8
        ENDDO D220
      ENDDO D7960
    ENDDO D7965
  ELSE
!
!       SCHEDULED IRRIGATION
!
    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf),&
      trim(mod_filename))

    do kk=1,24
      if(len_trim(irrigf(kk))==0)exit
!
!     DY,RR,JST,JEN=date DDMMYYYY,amount (mm),start and end hours
!     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,
!     CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
!     concentration in irrigation water
!
      READ(irrigf(kk),*)DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX &
          ,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      READ(irrigf(kk),'(I2,I2,I4)')IDY1,IDY2,IDY3

      IF(lverb)then
        write(*,*)irrigf(kk)
        write(*,*)IDY1,IDY2,IDY3,DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX &
            ,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      endif
      LPY=0
      IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
      IF(IDY2.EQ.1)then
        IDY=IDY1
      else
        IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      endif
      RRH=RR/(JEN-(JST-1))
      D8965: DO NX=NH1,NH2
        D8960: DO NY=NV1,NV2
          D2535: DO J=1,24
            IF(J.GE.JST.AND.J.LE.JEN)RRIG(J,IDY,NY,NX)=RRH/1000.0_r8
          ENDDO D2535
!
!           TRANSFER INPUTS TO MODEL ARRAYS
!
          PHQ(IDY,NY,NX)                                    = PHQX
          NH4_irrig_mole_conc(IDY,NY,NX)                    = CN4QX/14.0_r8
          NO3_irrig_mole_conc(IDY,NY,NX)                    = CNOQX/14.0_r8
          H2PO4_irrig_mole_conc(IDY,NY,NX)                  = CPOQX/31.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Al,IDY,NY,NX)  = CALQX/27.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Fe,IDY,NY,NX)  = CFEQX/55.8_r8
          trcsalt_irrig_mole_conc_col(idsalt_Ca,IDY,NY,NX)  = CCAQX/40.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Mg,IDY,NY,NX)  = CMGQX/24.3_r8
          trcsalt_irrig_mole_conc_col(idsalt_Na,IDY,NY,NX)  = CNAQX/23.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_K,IDY,NY,NX)   = CKAQX/39.1_r8
          trcsalt_irrig_mole_conc_col(idsalt_SO4,IDY,NY,NX) = CSOQX/32.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Cl,IDY,NY,NX)  = CCLQX/35.5_r8
          WDPTH(IDY,NY,NX)=WDPTHI
        ENDDO D8960
      ENDDO D8965
    enddo
  ENDIF
  end subroutine ReadIrrigationFile

!------------------------------------------------------------------------------------------

  subroutine ReadFertlizerFile(soilmgmt_nfid,FertFile,NH1,NH2,NV1,NV2)

  implicit none
  character(len=*), intent(in) :: FertFile
  integer,intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid

  integer :: NY,NX,LPY
  integer :: IDY1,IDY2,IDY3,IDY
  real(r8) :: FDPTHI,DY
  real(r8) :: Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB,Z4A
  real(r8) :: PMA,PMB,PHA,CAS,CAC
  real(r8) :: RSC1,RSN1,RSP1,RSC2,RSN2,RSP2
  real(r8) :: ROWX
  integer  :: IR0,IR1,IR2
  integer  :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=128) :: fertf(12)
  integer :: nfert

  nfert=get_dim_len(soilmgmt_nfid,'nfert')
  if(nfert>12)then
    call endrun('Not enough size for array fertf in '//trim(mod_filename), __LINE__)
  endif
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

  call check_var(soilmgmt_nfid, FertFile, vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find '//trim(FertFile)//' in '//trim(mod_filename), __LINE__)
  endif

  call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, fertf),&
      trim(mod_filename))

  do kk=1,12
    if(len_trim(fertf(kk))==0)exit

    READ(fertf(kk),*)DY,Z4A,Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB &
      ,PMA,PMB,PHA,CAC,CAS,RSC1,RSN1,RSP1,RSC2,RSN2,RSP2,FDPTHI &
      ,ROWX,IR0,IR1,IR2

    LPY=0
    IDY1=INT(DY/1.0E+06_r8)
    !return for bad fertilization data
    if(IDY1==0)return
    IDY2=INT(DY/1.0E+04_r8-IDY1*1.0E+02_r8)
    IDY3=INT(DY-(IDY1*1.0E+06_r8+IDY2*1.0E+04_r8))       
     
    IF(LVERB)then
      print*,fertf(kk)
      PRINT*,IDY1,IDY2,IDY3,Z4A,Z3A,ZUA,ZOA,Z4B,Z3B,ZUB,ZOB &
        ,PMA,PMB,PHA,CAC,CAS,RSC1,RSN1,RSP1,RSC2,RSN2,RSP2,FDPTHI &
        ,ROWX,IR0,IR1,IR2
    endif
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif

    D8985: DO NX=NH1,NH2
      D8980: DO NY=NV1,NV2
!
!         ENTER AMENDMENTS INTO MODEL ARRAYS
!
!         NH4,NH3,UREA,NO3 BROADCAST (A) AND BANDED (B)
!
        FERT(ifert_nh4,IDY,NY,NX)       = Z4A        !NH4 broadcast
        FERT(ifert_nh3,IDY,NY,NX)       = Z3A        !NH3 broadcast
        FERT(ifert_urea,IDY,NY,NX)      = ZUA        !Urea broadcast
        FERT(ifert_no3,IDY,NY,NX)       = ZOA        !NO3 broadcast
        FERT(ifert_nh4_band,IDY,NY,NX)  = Z4B        !NH4 band
        FERT(ifert_nh3_band,IDY,NY,NX)  = Z3B        !NH3 band
        FERT(ifert_urea_band,IDY,NY,NX) = ZUB        !Urea band
        FERT(ifert_no3_band,IDY,NY,NX)  = ZOB        !NO3 band
!
!         MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE BROADCAST (A)
!         AND BANDED (B)
!
        FERT(ifert_P_Ca_H2PO4_2,IDY,NY,NX)      = PMA        !PM broadcast
        FERT(ifert_P_Ca_H2PO4_2_band,IDY,NY,NX) = PMB        !PM band
        FERT(ifert_P_apatite,IDY,NY,NX)         = PHA        !
!
!         LIME AND GYPSUM
!
        FERT(ifert_Ca_lime,IDY,NY,NX)   = CAC
        FERT(ifert_Ca_gypsum,IDY,NY,NX) = CAS
!
!         PLANT AND ANIMAL RESIDUE C, N AND P
!
        FERT(ifert_plant_resC,IDY,NY,NX)  = RSC1
        FERT(ifert_plant_resN,IDY,NY,NX)  = RSN1
        FERT(ifert_plant_resP,IDY,NY,NX)  = RSP1
        FERT(ifert_plant_manuC,IDY,NY,NX) = RSC2
        FERT(ifert_plant_manuN,IDY,NY,NX) = RSN2
        FERT(ifert_plant_manuP,IDY,NY,NX) = RSP2
!
!         DEPTH AND WIDTH OF APPLICATION
!
        FDPTH(IDY,NY,NX) = FDPTHI    !depth
        ROWI(IDY,NY,NX)  = ROWX       !width
!
!         TYPE OF FERTILIZER,PLANT OR ANIMAL RESIDUE
!
        IYTYP(0,IDY,NY,NX)=IR0
        IYTYP(1,IDY,NY,NX)=IR1
        IYTYP(2,IDY,NY,NX)=IR2
      ENDDO D8980
    ENDDO D8985
  enddo
  end subroutine ReadFertlizerFile

!------------------------------------------------------------------------------------------

  subroutine ReadManagementFiles(yeari)
  !
  !DESCRIPTION
  !read in management information for soil: tillage
  !fertilizer and irrigation
  use EcoSIMCtrlMod, only : soil_mgmt_in,Lirri_auto
  implicit none
  integer, intent(in) :: yeari

  integer :: NH1,NV1,NH2,NV2
  type(file_desc_t) :: soilmgmt_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: ntopou
  integer :: NTOPO
  character(len=10) :: fertf
  character(len=10) :: tillf
  character(len=10) :: irrigf
  integer :: iyear,year
!
!   NH1,NV1,NH2,NV2=N,W and S,E corners of landscape unit
!   DATA1(8),DATA1(5),DATA1(6)=disturbance,fertilizer,irrigation files
!   PREFIX=path for files in current or higher level directory

  call ncd_pio_openfile(soilmgmt_nfid, soil_mgmt_in, ncd_nowrite)

  iyear=1
  DO while(.true.)
    call ncd_getvar(soilmgmt_nfid,'year',iyear,year)
    if(year==yeari)exit
    iyear=iyear+1
  ENDDO

  ntopou=get_dim_len(soilmgmt_nfid, 'ntopou')
  if(ntopou==0)return
  if(first_topou)ntopou=1  
  DO NTOPO=1,ntopou
    call ncd_getvar(soilmgmt_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(soilmgmt_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(soilmgmt_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(soilmgmt_nfid,'NV2',ntopo,NV2)

    if(any((/NH1,NH2,NV1,NV2/)<0))THEN
      call endrun('something wrong in NHX or NHX indices on file '//trim(soil_mgmt_in)&
        //' in '//trim(mod_filename), __LINE__)
    endif

    call check_var(soilmgmt_nfid, 'fertf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find fertf in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, fertf, &
      start = (/1,ntopou,iyear/),count = (/len(fertf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
!!!!
    call check_var(soilmgmt_nfid, 'tillf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find tillf in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, tillf, &
      start = (/1,ntopou,iyear/),count = (/len(tillf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
!!!!
    call check_var(soilmgmt_nfid, 'irrigf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find irrigf in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf, &
      start = (/1,ntopou,iyear/),count = (/len(irrigf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
!
!   READ TILLAGE INPUT FILE
!
    IF(trim(tillf).NE.'NO')THEN
      if(lverb)print*,'ReadTillageFile'
      call ReadTillageFile(soilmgmt_nfid,tillf,NH1,NH2,NV1,NV2)
    ENDIF
!
!   READ FERTLIZER INPUT FILE
!
    IF(trim(fertf).NE.'NO')THEN
      if(lverb)print*,'ReadFertlizerFile'
      call ReadFertlizerFile(soilmgmt_nfid,fertf,NH1,NH2,NV1,NV2)
    ENDIF
!
!   READ IRRIGATION INPUT FILE
!
    IF(trim(irrigf).NE.'NO')THEN
      Lirri_auto=irrigf(1:4)=='auto'
      if(lverb)print*,'ReadIrrigationFile'
      call ReadIrrigationFile(soilmgmt_nfid,irrigf,NH1,NH2,NV1,NV2)
    ENDIF
  ENDDO
  call ncd_pio_closefile(soilmgmt_nfid)

  end subroutine ReadManagementFiles
end module ReadManagementMod

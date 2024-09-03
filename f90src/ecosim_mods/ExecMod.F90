module ExecMod
!!
! Description:
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : endrun, padr, print_info
  use GridConsts
  use FlagDataType
  use EcoSimConst
  use EcoSIMCtrlMod     , only : etimer
  use EcoSIMCtrlDataType
  use EcoSimSumDataType
  use GridDataType
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), SAVE :: TLW,TLH,TLO,TLC,TLN,TLP,TLI
  real(r8) :: DIFFQ,DIFFH,DIFFO,DIFFC,DIFFN,DIFFP,DIFFI

  public :: exec
  contains

  SUBROUTINE exec(I)
!!
! Description:
! THIS SUBROUTINE TAKES MASS BALANCE VARIABLES CALCULATED
! IN 'REDIST' AND PERFORMS MASS BALANCE CHECKS AT THE END
! OF EACH DAY OF THE MODEL RUN, AND ERRORS OF > 1UG ARE FLAGGED.
!
  implicit none

  integer, intent(in) :: I

! begin_execution

! CALCULATE MASS BALANCES FOR WATER, HEAT, O2, C, N, P AND SOLUTES
!
  IF(I.EQ.IBEGIN.OR.I.EQ.ISTART.OR.I.EQ.ILAST+1)THEN
    TLW=WatMassStore_lnd-CRAIN+CRUN+CEVAP+QH2OLoss_lnds
!    if(tlw/=tlw)then
!      call print_info('tlw/=tlw',(/padr('WatMassStore_lnd',10),padr('CRAIN',10), &
!      padr('CRUN',10),padr('CEVAP',10),padr('QH2OLoss_lnds',10)/), &
!      (/WatMassStore_lnd,CRAIN,CRUN,CEVAP,QH2OLoss_lnds/))
!    endif
    TLH=HeatStore_lnd-HEATIN_lnd+HeatOut_lnds
    TLO=TSoilO2G_lnd-SurfGas_O2_lnd+OXYGOU
    TLC=LitRMStoreLndscap(ielmc)+POMHumStoreLndscap(ielmc)+TGasC_lnd-SurfGas_CO2_lnd+TOMOU_lnds(ielmc)-TORGF-Litrfall_lnds(ielmc)
    TLN=LitRMStoreLndscap(ielmn)+POMHumStoreLndscap(ielmn)+TGasN_lnd+TDisolNH4_lnd+tNO3_lnd-SurfGas_N2_lnd-TZIN+TOMOU_lnds(ielmn)-TORGN-Litrfall_lnds(ielmn)
    TLP=LitRMStoreLndscap(ielmp)+POMHumStoreLndscap(ielmp)+TDisolPi_lnd-TPIN+TOMOU_lnds(ielmp)-TORGP-Litrfall_lnds(ielmp)
    TLI=TION-TIONIN+TIONOU
  ENDIF
!
! CALCULATE DEVIATION SINCE MASS BALANCE WAS LAST RESET
!
  IF(etimer%its_time_to_diag())THEN
    DIFFQ=(WatMassStore_lnd-CRAIN+CRUN+CEVAP+QH2OLoss_lnds-TLW)/TAREA
    DIFFH=(HeatStore_lnd-HEATIN_lnd+HeatOut_lnds-TLH)/TAREA
    DIFFO=(TSoilO2G_lnd-SurfGas_O2_lnd+OXYGOU-TLO)/TAREA
    DIFFC=(LitRMStoreLndscap(ielmc)+POMHumStoreLndscap(ielmc)+TGasC_lnd-SurfGas_CO2_lnd+TOMOU_lnds(ielmc)-TORGF-Litrfall_lnds(ielmc)-TLC)/TAREA
    DIFFN=(LitRMStoreLndscap(ielmn)+POMHumStoreLndscap(ielmn)+TGasN_lnd+TDisolNH4_lnd+tNO3_lnd-SurfGas_N2_lnd-TZIN+TOMOU_lnds(ielmn) &
      -TORGN-Litrfall_lnds(ielmn)-TLN)/TAREA
    DIFFP=(LitRMStoreLndscap(ielmp)+POMHumStoreLndscap(ielmp)+TDisolPi_lnd-TPIN+TOMOU_lnds(ielmp)-TORGP-Litrfall_lnds(ielmp)-TLP)/TAREA
    DIFFI=(TION-TIONIN+TIONOU-TLI)/TAREA
    WRITE(*,212)I,iYearCurrent
    WRITE(18,213)I,iYearCurrent,DIFFQ,DIFFH,DIFFO,DIFFC,DIFFN &
      ,DIFFP,DIFFI
    if(diffq/=diffq)then
      write(*,*)'DIFFQ=',DIFFQ
      write(*,*)'WatMassStore_lnd=',WatMassStore_lnd
      write(*,*)'CRAIN=',CRAIN
      write(*,*)'CRUN=',CRUN
      write(*,*)'CEVAP=',CEVAP
      write(*,*)'QH2OLoss_lnds=',QH2OLoss_lnds
      write(*,*)'TLW=',TLW
      write(*,*)'TAREA=',TAREA
      call endrun(msg='NaN encounterd in '//mod_filename,line=__LINE__)
    endif
212   FORMAT('NOW EXECUTING DAY',I6,'   OF YEAR',I6)
213   FORMAT(2I6,10F16.6)
!
!     FLAG DEVIATIONS > 1UG OR 1 J IN ENTIRE MODEL LANDSCAPE,
!     RESET MASS BALANCE
!
    IF(ABS(DIFFC).GT.ppmc)THEN
      WRITE(18,191)I,iYearCurrent
191   FORMAT('CARBON BALANCE LOST ON DAY, YEAR',2I4)
      TLC=LitRMStoreLndscap(ielmc)+POMHumStoreLndscap(ielmc)+TGasC_lnd-SurfGas_CO2_lnd+TOMOU_lnds(ielmc)-TORGF-Litrfall_lnds(ielmc)
    ENDIF
    IF(ABS(DIFFN).GT.ppmc)THEN
      WRITE(18,192)I,iYearCurrent
192   FORMAT('NITROGEN BALANCE LOST ON DAY, YEAR',2I4)
      TLN=LitRMStoreLndscap(ielmn)+POMHumStoreLndscap(ielmn)+TGasN_lnd+TDisolNH4_lnd+tNO3_lnd-SurfGas_N2_lnd-TZIN+TOMOU_lnds(ielmn)-TORGN-Litrfall_lnds(ielmn)
    ENDIF
    IF(ABS(DIFFP).GT.ppmc)THEN
      WRITE(18,193)I,iYearCurrent
193   FORMAT('PHOSPHORUS BALANCE LOST ON DAY, YEAR',2I4)
      TLP=LitRMStoreLndscap(ielmp)+POMHumStoreLndscap(ielmp)+TDisolPi_lnd-TPIN+TOMOU_lnds(ielmp)-TORGP-Litrfall_lnds(ielmp)
    ENDIF
    IF(ABS(DIFFQ).GT.ppmc)THEN
      WRITE(18,194)I,iYearCurrent
194   FORMAT('WATER BALANCE LOST ON DAY, YEAR',2I4)
      TLW=WatMassStore_lnd-CRAIN+CRUN+CEVAP+QH2OLoss_lnds
    ENDIF
    IF(ABS(DIFFH).GT.ppmc)THEN
      WRITE(18,195)I,iYearCurrent
195   FORMAT('THERMAL BALANCE LOST ON DAY, YEAR',2I4)
      TLH=HeatStore_lnd-HEATIN_lnd+HeatOut_lnds
    ENDIF
    IF(ABS(DIFFO).GT.ppmc)THEN
      WRITE(18,196)I,iYearCurrent
196   FORMAT('OXYGEN BALANCE LOST ON DAY, YEAR',2I4)
      TLO=TSoilO2G_lnd-SurfGas_O2_lnd+OXYGOU
    ENDIF
    IF(ABS(DIFFI).GT.ppmc)THEN
      WRITE(18,197)I,iYearCurrent
197   FORMAT('ION BALANCE LOST ON DAY, YEAR',2I4)
      TLI=TION-TIONIN+TIONOU
    ENDIF
  ENDIF
  IF(IDAYR.LT.0)THEN
    IDAYR=LYRX+IDAYR
  ELSE
    IDAYR=I
  ENDIF
  IOLD=I
  IMNG=0
  NYR=0
  RETURN

  END subroutine exec
end module ExecMod

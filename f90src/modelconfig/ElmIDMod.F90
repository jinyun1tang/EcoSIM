module ElmIDMod
!
! DESCRIPTION:
!  Chemical element ids
implicit none
  integer, parameter :: ielmc=1    !carbon element
  integer, parameter :: ielmn=2    !nitrogen element
  integer, parameter :: ielmp=3    !phosphorus element
  integer, parameter :: npelms=3   !totally three elements
! photosynthesis
  integer, parameter :: ic4_photo=4
  integer, parameter :: ic3_photo=3
! fertilizer
  integer, parameter :: ifert_nh4=1
  integer, parameter :: ifert_nh3=2
  integer, parameter :: ifert_urea=3
  integer, parameter :: ifert_no3=4
  integer, parameter :: ifert_nh4_band=5
  integer, parameter :: ifert_nh3_band=6
  integer, parameter :: ifert_urea_band=7
  integer, parameter :: ifert_no3_band=8
  integer, parameter :: ifert_P_Ca_H2PO4_2=9
  integer, parameter :: ifert_P_Ca_H2PO4_2_band=10
  integer, parameter :: ifert_P_apatite=11
  integer, parameter :: ifert_Ca_lime=12
  integer, parameter :: ifert_Ca_gypsum=13
  integer, parameter :: ifert_plant_resC=14
  integer, parameter :: ifert_plant_resN=15
  integer, parameter :: ifert_plant_resP=16
  integer, parameter :: ifert_plant_manuC=17
  integer, parameter :: ifert_plant_manuN=18
  integer, parameter :: ifert_plant_manuP=19
! root order
  integer, parameter :: fstroot =1
  integer, parameter :: sndroot =2
! root profile
  integer, parameter :: irtshallow=0
! mycorrhizae
  integer, parameter :: mycornone=1
  integer, parameter :: mycorarbu=2   !arbuscular

  integer, parameter :: ifineroot=1
  integer, parameter :: imycorrhz=2

  integer, parameter :: iliving_branch=0

  integer, parameter :: itrue=1
  integer, parameter :: ifalse=0

  integer, parameter :: ialive=0
  integer, parameter :: idead =1

  integer, parameter :: ihv_noaction=0
  integer, parameter :: ihv_terminate=1
  integer, parameter :: ihv_tmareseed=2

  integer, parameter :: iplt_annual=0
  integer, parameter :: iplt_preanu=1

  integer, parameter :: iplt_bryophyte=0
  integer, parameter :: iplt_grasslike=1
  integer, parameter :: iplt_treelike=2


end module ElmIDMod

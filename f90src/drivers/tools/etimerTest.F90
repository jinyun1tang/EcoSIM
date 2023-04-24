PROGRAM main

  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod     , only : etimer
  use ecosim_time_mod   , only : getdow,ecosim_time_dat_type,get_steps_from_ymdhs
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__

  type(ecosim_time_dat_type)  :: etime_dat

  real(r8) :: a(3)

  integer :: nyr, J,mon
  character(len=16) :: ymdhs

  print*,'1900,14 ',getdow(1900,14)
  print*,'1999,14 ',getdow(1999,14)
  print*,'1699,14 ',getdow(1699,14)
  print*,'2099,14 ',getdow(2099,14)

  call etimer%Init(year0=1804,nyears=12)

  call etimer%setClock(dtime=3600._r8,nelapstep=0)

  DO while(.true.)
    do while(.true.)
      call etimer%get_ymdhs(ymdhs)
      
      DO J=1,24
        call etimer%update_time_stamp()
      ENDDO

      if(etimer%its_a_new_month())print*,ymdhs,etimer%get_curr_day()
      if(etimer%its_a_new_year())exit
    enddo
    if(etimer%its_time_to_exit())exit
  enddo

  call get_steps_from_ymdhs('18810315060000',3600,etime_dat,1880) 

  a=(/1.,2.,3./)
  print*,'18810315060000 ',sum(a(1:0))
  print*,etime_dat
end program main
PROGRAM main

  use data_kind_mod     , only : r8 => SHR_KIND_R8
  use EcoSIMCtrlMod     , only : etimer
  use ecosim_time_mod   , only : getdow
  implicit none

  character(len=*), parameter :: mod_filename = __FILE__

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

      if(etimer%its_a_new_month())print*,ymdhs,etimer%get_cur_day()
      if(etimer%its_a_new_year())exit
    enddo
    if(etimer%its_time_to_exit())exit
  enddo
end program main
program main
  
  use run_angio_m
  use global_m

  implicit none

  call get_command_argument(1,sim_id)
  call run_angio(sim_id)
 
  

end program main



program main
  
  use run_angio_m
  use mpi

  implicit none

  integer :: rank, size, ierr
  character(5) :: input, dir
  character(1) :: n2s1
  character(2) :: n2s2


  call MPI_INIT(ierr)

  call MPI_COMM_RANK(mpi_comm_world, rank, ierr)
  call MPI_COMM_SIZE(mpi_comm_world, size, ierr)



  if(rank<10) then
     write(n2s1,'(I1)')rank
     input = 'inp0'//n2s1
     dir = 'dir0'//n2s1     
  else
     write(n2s2,'(I2)') rank
     input = 'inp'//n2s2
     dir = 'dir'//n2s2
  end if

  call system('mkdir '//trim(dir))

  call MPI_BARRIER(mpi_comm_world, ierr)


  call run_angio(input)
 
  

end program main



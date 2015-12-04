program main

  use get_data_m

  implicit none
  integer :: Lsize(3), np, boundary_points, ndim, np_part
  integer :: lines_phi, lines_phis, initial_file, last_file, delta, ifile, ip
  integer :: i,j,k, output
  real, allocatable ::  phi(:), phis(:)
  integer, allocatable :: lxyz(:,:), lxyz_inv(:,:,:)
  real :: f, volume, length, diameter, M_PI
  logical :: periodic 

  character (8) :: temp
  character (3) :: dir_name
  character (5) :: file_id5
  character (6) :: file_id6
  character (100) :: looong_string
  
  write(*,*) "initial file, last_file, delta, dir_name:"
  read(*,*) initial_file, last_file, delta, dir_name
  
  M_PI = 3.14159265359
  Lsize(1) = 50
  Lsize(2) = 50
  Lsize(3) = 30
  
  periodic = .true.
  boundary_points = 2
  
  np = 8 * Lsize(1) * Lsize(2) * Lsize(3)
  np_part = 8*(Lsize(1)+2*boundary_points)*(Lsize(2)+2*boundary_points)*(Lsize(3)+2*boundary_points)
  
  ALLOCATE(lxyz(np_part,1:3))
  ALLOCATE(lxyz_inv(-Lsize(1)-boundary_points:Lsize(1)+boundary_points, &
         -Lsize(2)-boundary_points:Lsize(2)+boundary_points, &
         -Lsize(3)-boundary_points:Lsize(3)+boundary_points))  
  
  ALLOCATE(phi(np))
  ALLOCATE(phis(np))
  
  write(*,*) "---------------------------------------------------------"
  write(*,*) "                      get_data_m module                 "
  write(*,*) "                                                        " 
  write(*,*) " This program calculate the mean diameter and the       "  
  write(*,*) " the number of branches for each time step of a         "
  write(*,*) " simulation and plot two graphs in function of the time."
  write(*,*) "---------------------------------------------------------"

  call space_init(Lsize, lxyz, lxyz_inv, boundary_points, np, ndim, periodic)


  write(*,'(A)') "Calculating . . ."
  open(UNIT=15,FILE='data'//dir_name//'.dat')
  do ifile=initial_file, last_file, delta
     
     phis(:) = 0.0
     
     if(ifile>=100000) then
        write(file_id6,'(I6)') ifile
        call system('wc -l '//dir_name//'/phi'//file_id6//'.xyz'//' > lines_phi.aux')
        call system('wc -l '//dir_name//'/phis'//file_id6//'.xyz'//' > lines_phis.aux')
     else
        write(file_id5,'(I5)') ifile
        call system('wc -l '//dir_name//'/phi'//file_id5//'.xyz'//' > lines_phi.aux')
        call system('wc -l '//dir_name//'/phis'//file_id5//'.xyz'//' > lines_phis.aux')
     end if
     
     open(UNIT=100, FILE='lines_phi.aux')
     read(100,*) lines_phi, temp
     close(100)

     open(UNIT=200, FILE='lines_phis.aux')
     read(200,*) lines_phis, temp
     close(200)

     if(ifile>=100000) then 
        open(UNIT=ifile, FILE=dir_name//'/phis'//file_id6//'.xyz')
     else
        open(UNIT=ifile, FILE=dir_name//'/phis'//file_id5//'.xyz')
     end if
     
     do ip=1, lines_phis
        read(ifile,*) i, j, k, f
        phis(lxyz_inv(i,j,k)) = 1.0
     end do
     
     close(ifile)
     
     call run_get_data(phis,Lsize,lxyz,lxyz_inv,np,periodic,output)
     
     diameter = 2.00*sqrt(lines_phi/(M_PI*lines_phis))
     
     write(15,*)  ifile, output, diameter
     
  end do
  close(15)
  
  
  open(UNIT=16,FILE='plot_file.aux')
  write(16,*) "set xlabel 'time step'"
  write(16,*) "set ylabel 'branches'"
  write(16,*) "set key bottom"
  write(16,*) "set grid"
  write(16,*) "set term pdf"
  write(16,*) "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1"
  write(16,*) "set output 'branches"//dir_name//".pdf'"
  write(16,*) "plot 'data"//dir_name//".dat' using 1:2 title 'branches' w linespoint ls 7"
  write(16,*) "set ylabel 'diameter'"
  write(16,*) "set output 'diameter"//dir_name//".pdf'"  
  write(16,*) "plot 'data"//dir_name//".dat' using 1:3 title 'diameter' w linespoint ls 2"
  write(16,*) "set term wxt"
  close(16)
  
  call system('gnuplot < plot_file.aux')
  
  write(*,*) "-------------------------------------------------------"
  write(*,*) " "
  write(*,*) "data file: data"//dir_name//".dat" 
  write(*,*) "plotting..."
  write(*,*) "           branches"//dir_name//".pdf... done!"
  write(*,*) "           diameter"//dir_name//".pdf... done!"  
  write(*,*) " "
  write(*,*) "------------------------------------------------------"
  write(*,*) "End program. "
  
  DEALLOCATE(phi)
  DEALLOCATE(phis)
  DEALLOCATE(lxyz)
  DEALLOCATE(lxyz_inv)


end program main



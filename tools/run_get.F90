program main

  use get_data_m

  implicit none
  integer :: Lsize(3), np, boundary_points, ndim, np_part
  integer :: lines_phi, lines_phis, initial_file, last_file, delta, ifile, ip
  integer :: i,j,k, output, genus
  real, allocatable ::  phi(:), phis(:)
  integer, allocatable :: lxyz(:,:), lxyz_inv(:,:,:)
  real :: f, volume, length, diameter, M_PI
  logical :: periodic

  character (8) :: temp
  character (3) :: dir_name
  character (4) :: file_id4
  character (5) :: file_id5
  character (6) :: file_id6
  character (100) :: looong_string

  write(*,*) "initial file, last_file, delta, dir_name:"
  read(*,*) initial_file, last_file, delta, dir_name

  M_PI = 3.14159265359
  Lsize(1) = 100
  Lsize(2) = 100
  Lsize(3) = 50

  periodic = .false.
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
  write(*,*) " simulation and plot two graphs as a function of time.  "
  write(*,*) "---------------------------------------------------------"

  call space_init(Lsize, lxyz, lxyz_inv, boundary_points, np, ndim, periodic)



  write(*,'(A)') "Calculating . . ."
  open(UNIT=15,FILE='data'//dir_name//'.dat')
  do ifile=initial_file, last_file, delta

     phis(:) = 0.0
     write(*,*) ifile
     if(ifile>=100000) then
        write(file_id6,'(I6)') ifile
        call system('wc -l '//dir_name//'/phi'//file_id6//'.xyz'//' > lines_phi.aux')
        call system('wc -l '//dir_name//'/phis'//file_id6//'.xyz'//' > lines_phis.aux')
     elseif(ifile>=10000.and.ifile<100000) then
        write(file_id5,'(I5)') ifile
        call system('wc -l '//dir_name//'/phi'//file_id5//'.xyz'//' > lines_phi.aux')
        call system('wc -l '//dir_name//'/phis'//file_id5//'.xyz'//' > lines_phis.aux')
     else
       write(file_id4,'(I4)') ifile
       call system('wc -l '//dir_name//'/phi'//file_id4//'.xyz'//' > lines_phi.aux')
       call system('wc -l '//dir_name//'/phis'//file_id4//'.xyz'//' > lines_phis.aux')
     end if

     open(UNIT=100, FILE='lines_phi.aux')
     read(100,*) lines_phi, temp
     close(100)

     open(UNIT=200, FILE='lines_phis.aux')
     read(200,*) lines_phis, temp
     close(200)

     if(ifile>=100000) then
        open(UNIT=ifile, FILE=dir_name//'/phis'//file_id6//'.xyz')
     elseif(ifile>=10000.and.ifile<100000) then
        open(UNIT=ifile, FILE=dir_name//'/phis'//file_id5//'.xyz')
     else
       open(UNIT=ifile, FILE=dir_name//'/phis'//file_id4//'.xyz')
     end if

     phis(:) = 0.d0
      ! do k = -50,49
      !   phis(lxyz_inv(0,0,k)) = 1.0
      ! end do
      !
      ! do i = 1, 10
      !   phis(lxyz_inv(i,0,10)) = 1.0
      !   phis(lxyz_inv(i,0,-10)) = 1.0
      !   phis(lxyz_inv(-i,0,-20)) = 1.0
      !   phis(lxyz_inv(i,i,30)) = 1.0
      !   phis(lxyz_inv(i,i,-20)) = 1.0
      !   phis(lxyz_inv(i,i,i)) = 1.0
      !   do j=5,10
      !     phis(lxyz_inv(i,j,30)) = 1.0
      !   end do
      ! end do
      !
      ! do i=-10,10
      !   phis(lxyz_inv(10,0,i)) = 1.0
      ! end do
      ! phis(lxyz_inv(10,0,10)) = 0.0
      ! phis(lxyz_inv(10,0,-10)) = 0.0
      ! phis(lxyz_inv(0,0,10)) = 0.0
      ! phis(lxyz_inv(0,0,-10)) = 0.0
      ! phis(lxyz_inv(0,0,-20)) = 0.0
      ! phis(lxyz_inv(0,0,1)) = 0.0
       open(UNIT=600, file='phis150000.xyz')
     do ip=1, lines_phis

        read(ifile,*) i, f
        phis(i) = 1.0
      !if(phis(ip)>0)   write(600,*) lxyz(ip,1:3), phis(ip)
      write(600,*) lxyz(i,1:3),f
        !read(ifile,*) i, j, k, f
        !phis(lxyz_inv(i,j,k)) = 1.0
     end do
     close(600)
     close(ifile)


     call run_get_data(phis,Lsize,lxyz,lxyz_inv,np,periodic,output,genus)

     diameter = 2.00*sqrt(lines_phi/(M_PI*lines_phis))

     write(15,*)  ifile, output, diameter, genus

  end do
  close(15)

!  go to 1
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
  write(16,*) "set ylabel 'Genus'"
  write(16,*) "set key bottom"
  write(16,*) "set output 'genus"//dir_name//".pdf'"
  write(16,*) "plot 'data"//dir_name//".dat' using 1:4 title 'Genus' w linespoint ls 2"
  write(16,*) "set term wxt"
  close(16)

  call system('gnuplot < plot_file.aux')

  write(*,*) "-------------------------------------------------------"
  write(*,*) " "
  write(*,*) "data file: data"//dir_name//".dat"
  write(*,*) "plotting..."
  write(*,*) "           branches"//dir_name//".pdf... done!"
  write(*,*) "           diameter"//dir_name//".pdf... done!"
  write(*,*) "           genus"//dir_name//".pdf...    done!"
  write(*,*) " "
  write(*,*) "------------------------------------------------------"
  write(*,*) "End program. "
! 1 write(*,*) "branchs:", output
  DEALLOCATE(phi)
  DEALLOCATE(phis)
  DEALLOCATE(lxyz)
  DEALLOCATE(lxyz_inv)


end program main

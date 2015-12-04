!! Copyright (C) 2015 M. Moreira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! 

module run_angio_m  

  use thinning_m
 
  
  implicit none
  
  private
  
  public :: run_angio
  
  type tip_cell_t
     integer :: ip    ! ip global
     real :: x, y, z ! global mesh position 
     real :: phi     ! phi_c inside tipcell
  end type tip_cell_t
  
  type mesh_t
     real :: phi, T, mu, lapl_phi, lapl_T, lapl_mu
     real :: alpha_p, consumption
     integer :: source
  end type mesh_t
  
  type sec_mesh_t
     real :: phi
  end type sec_mesh_t
  
  contains
  
    subroutine run_angio()

      implicit none
      
      
      ! begin parameters
      real :: cell_radius, diffusion_const, interface_width, vegf_p, vegf_c, diff_oxy_length, vegf_rate, &
           vegf_source_conc, prolif_rate, vessel_radius, dt, chi, vegf_grad_min, vegf_grad_max
      integer :: tstep, iseed, dr(3)
      character(len=3) :: dir_name 
      character(len=10) :: file_name, file_id
      real :: notch_distance, depletion_weight
      ! mesh variables
      integer, allocatable :: lxyz(:,:), lxyz_inv(:,:,:),  grid_cell_domain(:)
      integer :: Lsize(1:3), boundary_points
      integer :: np, np_part, ip, n_max_tipc
      real, allocatable :: gg(:,:), lapl(:), f(:)
      type(mesh_t), allocatable :: cell(:)
      type(sec_mesh_t), allocatable :: necrotic_tissue(:)
      type(tip_cell_t), allocatable :: tipc(:)
      ! vegf related variables
      integer :: n_source = 0, source_max, n_source_initial
      integer, allocatable :: vegf_xyz(:,:)
      ! misc
      integer :: nstep, counter = 0, ndim, output_period, extra_steps
      real :: hs, time_init, time_end, ctime, phi_max
      logical :: periodic
      ! thinning
      logical :: thinning
      real, allocatable :: phis(:)
      ! np_vegf_s  - number of points on the spherical surface representing the limit signaling of vegf source
      ! np_tip_s   - number of points on the spherical surface representing tip cell limits
      ! np_tip_all - all points inside tip cell
      ! n_tipcell - number of activated tip cells
      ! n_source  - number of vegf sources
      integer :: np_vegf_s, np_tip_s, np_tip_all, n_tipcell = 0
      integer, allocatable :: vegf_s(:,:), tip_s(:,:), tip_all(:,:)

      !temporary
      real :: phi_med, t_med

      ! initializing parameters
      call  parameters_init(cell_radius, diffusion_const, interface_width, vegf_p, vegf_c, diff_oxy_length,&
      vegf_rate, vegf_source_conc, prolif_rate, vessel_radius, tstep, dt, chi, Lsize, dr, dir_name, iseed,&
      boundary_points, source_max, vegf_grad_min, vegf_grad_max, depletion_weight, output_period, extra_steps, &
      n_max_tipc, thinning, periodic)

	
      ! notch distance
      notch_distance = (4.d0*cell_radius)**2
      ! minimum vegf gradient for branch ** 2
      vegf_grad_min = vegf_grad_min**2
      ! number of points in the mesh
      np = 8*Lsize(1)*Lsize(2)*Lsize(3) ! number of points
      np_part = 8*(Lsize(1)+2*boundary_points)*(Lsize(2)+2*boundary_points)*(Lsize(3)+2*boundary_points) ! number of points plus boundary points
      ! phi max inside etc
      phi_max = (prolif_rate*vegf_p*cell_radius)/(2.d0*chi*vegf_grad_max)
      
	
      ! allocating matrices and vectors
      ALLOCATE(lxyz(np_part,1:3))
      ALLOCATE(lxyz_inv(-Lsize(1)-boundary_points:Lsize(1)+boundary_points, &
                        -Lsize(2)-boundary_points:Lsize(2)+boundary_points, &
                        -Lsize(3)-boundary_points:Lsize(3)+boundary_points))
      ALLOCATE(cell(np))
      ALLOCATE(tipc(n_max_tipc))
      ALLOCATE(necrotic_tissue(np))
      ALLOCATE(vegf_xyz(1:source_max, 1:3))
      ALLOCATE(gg(1:np,1:3))
      ALLOCATE(lapl(1:np))
      ALLOCATE(f(1:np))
      ALLOCATE(phis(1:np))
      ! surfaces and points from spheres
      ALLOCATE(vegf_s(1:10000,1:3))
      ALLOCATE(tip_s(1:300,1:3))
      ALLOCATE(tip_all(1:500,1:3))
      


      ! initializing tip cell points     
      call spherical_surface(cell_radius, tip_s, np_tip_s, 0.2)  ! spherical surface increments of tip cell      
      call gen_cell_points(cell_radius, tip_all, np_tip_all)     ! all points of tip cell

      ! initializing points for the vegf sources 
!      call gen_cell_points(diff_oxy_length, vegf_s, np_vegf_s)
     call spherical_surface(diff_oxy_length, vegf_s, np_vegf_s, 0.2)  

      ! initializing space matrices
      call space_init(Lsize, lxyz, lxyz_inv, boundary_points, np, ndim, periodic)

      ! initializing simulation box
      call simul_box_init(cell, necrotic_tissue, tip_all, cell_radius, lxyz, lxyz_inv, Lsize, vegf_xyz, &
      vegf_source_conc, source_max, n_source, np, diff_oxy_length, vessel_radius, np_tip_all, ndim, iseed)

      ! saving the number of initial activated vegf sources
      n_source_initial = n_source

      ! Calculating distances between all points of the mesh and vegf sources
      ALLOCATE(grid_cell_domain(1:np))
      call gen_grid_cell_domain(grid_cell_domain, vegf_xyz, Lsize, lxyz, np, n_source, cell_radius, ndim)


      ! printing header 
      call print_header(Lsize, n_source, dir_name, periodic) 
  
      ! Cahn-Hilliard steps only for the necrotic tissue 
      call cahn_hilliard(necrotic_tissue, 100, dt, dr, np, interface_width, lxyz, lxyz_inv, dir_name)
  
      write(*,'(A)') "Running Cahn-Hilliard steps for the hypoxic tissue... done!" 

      ! init from file
      ! caution: this feature can't load the tip cells
      ! just the fields Phi and T
      ! use on your own risk!
      ! write(*,'(A)') "Loading previous system..."
      ! file_id = ' 73000'
      !                     field, ...,    ..., ... , directory, ...
      ! call init_from_file(cell, lxyz_inv, np, np_phi, 'db6' ,  file_id)


      write(*,'(A)') "Initiating the core program..."      

      ! saving the initial condition
      OPEN (UNIT=12,FILE=trim(dir_name//'/phii.xyz'))
      OPEN (UNIT=22,FILE=trim(dir_name//'/ti.xyz'))
      do ip=1, np
         if(cell(ip)%phi>0.d0) then
            write(12,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%phi
         end if
         write(22,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%T
      end do
      close(12)
      close(22)



      nstep = 0
      do while(nstep<=tstep)
         nstep = nstep + 1

         ! the program will run extra_steps steps after the last source deactivation
         if(n_source.eq.0) then
            tstep = nstep + extra_steps
            n_source = -1
         end if
         
       
         call CPU_TIME(time_init)
 
         ! calculating gradient of vegf
         call dderivatives_grad(cell, gg, np, lxyz, lxyz_inv, dr)
 
         ! ETC search
         call etc_search(cell, tipc, gg, Lsize, lxyz, lxyz_inv, notch_distance, n_tipcell, tip_s, &
                         vegf_c, nstep, vegf_grad_min, np, cell_radius, n_max_tipc, np_tip_s, periodic)

         ! calculating laplacian of phi
         f(1:np) = cell(1:np)%phi
         call dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)
         cell(1:np)%lapl_phi = lapl(1:np)


         ! Chemical potential:  phi**3 - phi - epsilon*laplacian(phi ) 
         do ip = 1, np
            cell(ip)%mu = cell(ip)%phi * cell(ip)%phi *cell(ip)%phi &
                 - cell(ip)%phi - interface_width*cell(ip)%lapl_phi&
                 + depletion_weight*(( cell(ip)%phi**2 - 1)*( ( 1.d0 + necrotic_tissue(ip)%phi )**2) &
                 *(necrotic_tissue(ip)%phi-2 )*(0.1875))

            ! Travasso 2D
            !+ depletion_weight*(cell(ip)%phi**2 - 1)*(( 1.d0 + necrotic_tissue(ip)%phi)**2) &
            !*(necrotic_tissue(ip)%phi-2 )*(0.1875) ) 
            ! Nonomura 3D
            !+ depletion_weight*(cell(ip)%phi-5.0)*(cell(ip)%phi + 1.0)*(( 1.d0 + necrotic_tissue(ip)%phi)**2) &
            !*(necrotic_tissue(ip)%phi-2.0 )*(0.0625)
         end do
       
         ! Calculating laplacian of mu 
         f(1:np) = cell(1:np)%mu
         call dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)
         cell(1:np)%lapl_mu = lapl(1:np)
    
         ! Calculating laplacian of T
         f(1:np) = cell(1:np)%T
         call dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)
         cell(1:np)%lapl_T = lapl(1:np)

         ! Temporal evolution phi(t+dt) and T(t+dt)
         
         do ip = 1, np

            ! calculating alpha_p
            if(cell(ip)%T>=vegf_p) then
               cell(ip)%alpha_p =  prolif_rate*vegf_p
            else
               cell(ip)%alpha_p = prolif_rate*cell(ip)%T
            end if
            
            hs = heaviside(cell(ip)%phi)
            
            if(cell(ip)%source<1) then
               cell(ip)%consumption = vegf_rate*cell(ip)%T*cell(ip)%phi*hs
               cell(ip)%T = cell(ip)%T + dt*(diffusion_const*cell(ip)%lapl_T - cell(ip)%consumption)
            end if
            
            cell(ip)%phi = cell(ip)%phi + dt*(cell(ip)%lapl_mu + cell(ip)%alpha_p*cell(ip)%phi*hs) 

         end do
        
   
         if (n_tipcell>0) then
            
            call etc_move(cell, tipc, gg, lxyz, lxyz_inv, Lsize, vegf_c, vegf_grad_min, vegf_grad_max, tip_s, &
                  tip_all, n_tipcell, phi_max, nstep, dt, chi, cell_radius, np_tip_all, grid_cell_domain, n_source_initial,&
                  vegf_xyz, ndim, np_tip_s, np, periodic)
           
         end if ! if n_tipcell > 0	 

         call source_deactivate(cell, vegf_xyz, n_source, vegf_s, lxyz, lxyz_inv, np_vegf_s, Lsize, periodic)

         call CPU_TIME(time_end)

         ctime = ctime + (time_end - time_init)
         
         if(nstep.eq.100) then
            ctime = (ctime*(tstep-nstep) )/6000.d0
            
            if( ctime>60.d0) then
               write(*,'(A,F10.2)') "Estimated time (hour): ",ctime/60.d0
            else
               write(*,'(A,F10.2)') "Estimated time (min): ",ctime
            end if
         end if
         
     

         ! output
         
         counter = counter + 1 

         if(n_tipcell.eq.n_max_tipc) then
            counter = output_period
         end if

         if(counter.eq.output_period) then
            write(*,*) nstep
            counter = 0

            if(thinning) then
               do ip=1, np
                  if(cell(ip)%phi>=0) then
                     phis(ip) = 1.d0
                  else
                     phis(ip) = -1.d0
                  end if
               end do
               
               call thinning_run(phis, lxyz, lxyz_inv, lsize, np)
            end if
            
            write(file_name,'(I6)') nstep
            OPEN (UNIT=nstep,FILE=dir_name//'/phi'//trim(file_name)//'.xyz')
            OPEN (UNIT=nstep+1,FILE=dir_name//'/t'//trim(file_name)//'.xyz')
            if(thinning) OPEN (UNIT=nstep+2,FILE=dir_name//'/phis'//trim(file_name)//'.xyz')
            do ip=1, np
               
               if(phis(ip)>0.and.thinning) then
                  write(nstep+2,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), phis(ip)
               end if
               
               if(cell(ip)%phi>0) then
                  write(nstep,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%phi
               end if
               write(nstep+1,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%T
            end do
            if(thinning) close(nstep+2)
            close(nstep+1)
            close(nstep)
         end if
         ! end of the output

         if(n_tipcell.eq.n_max_tipc) EXIT
         
     
      end do
      

      OPEN (UNIT=333,FILE=trim(dir_name//'/phi.xyz'))
      OPEN (UNIT=222,FILE=trim(dir_name//'/t.xyz'))
      do ip=1, np
         if(cell(ip)%phi>0.d0) then
            write(333,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%phi
         end if
         write(222,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), cell(ip)%T
      end do
      close(333)
      close(222)
      
      
      DEALLOCATE(lxyz)
      DEALLOCATE(lxyz_inv)
      DEALLOCATE(cell)
      DEALLOCATE(tipc)
      DEALLOCATE(necrotic_tissue)
      DEALLOCATE(vegf_xyz)
      DEALLOCATE(gg)
      DEALLOCATE(lapl)
      DEALLOCATE(f)
      DEALLOCATE(phis)
      DEALLOCATE(vegf_s)
      DEALLOCATE(tip_s)
      DEALLOCATE(tip_all)
      DEALLOCATE(grid_cell_domain)
      
     
      
    end subroutine run_angio


    subroutine init_from_file(cell, lxyz_inv, np, np_phi, dir_name, file_id)

      implicit none
      ! external
      type(mesh_t), allocatable, intent(inout) :: cell(:)
      integer, intent(in) :: np, np_phi
      integer, allocatable, intent(in) :: lxyz_inv(:,:,:)
      character(3), intent(in) :: dir_name
      character(10), intent(in) :: file_id
      ! private
      integer :: r(3), ip
      real :: phi_temp
      
      open(UNIT=100, FILE=dir_name//'/phi'//trim(file_id)//'.xyz')
      open(UNIT=200, FILE=dir_name//'/t'//trim(file_id)//'.xyz')
      cell(:)%phi = -1.d0
      do ip=1, np
         read(200,*) r(1:3), cell(ip)%T
      end do

      do ip=1, np_phi
         read(100,*) r(1:3), phi_temp

         cell(lxyz_inv(r(1),r(2),r(3) ) )%phi = phi_temp
      end do

      close(100)
      close(200)
      
    end subroutine init_from_file

    subroutine gen_grid_cell_domain(grid_cell_domain, vegf_xyz, Lsize, lxyz, np, n_source, cell_radius, ndim)
   
      implicit none

      integer, allocatable, intent(inout) :: grid_cell_domain(:)
      real, intent(in) :: cell_radius
      integer, allocatable, intent(in) :: lxyz(:,:), vegf_xyz(:,:)
      integer, intent(in) :: np, n_source, Lsize(3), ndim 
      integer :: i, j, unitary(3)
      real :: distance(3), distance_abs
      
      unitary(1:3) = 1
      grid_cell_domain(1:np) = -1

      do i=1, n_source
         do j=1, np
            distance(1:3) = abs(vegf_xyz(i,1:3) - lxyz(j,1:3))
            distance(1:ndim) = min(distance(1:ndim), 2*Lsize(1:ndim)-unitary(1:ndim)-abs(distance(1:ndim)))
            distance_abs = sqrt(distance(1)**2 + distance(2)**2 + distance(3)**2 )
            if(distance_abs<=2.d0*cell_radius) then
               grid_cell_domain(j) = i
            end if
         end do
      end do
      
    end subroutine gen_grid_cell_domain
     
    
    subroutine cahn_hilliard(field, nstep, dt, dr, np, interface_width, lxyz, lxyz_inv, dir_name)
      
      implicit none
      
      type(sec_mesh_t), allocatable, intent(inout) :: field(:)
      integer,allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) :: nstep, np, dr(3)
      real, intent(in) :: interface_width, dt
      character(len=3), intent(in) :: dir_name
      ! internal variables
      real, allocatable :: f(:), lapl(:)
      real, allocatable :: mu(:)
      integer :: tstep, ip
      
      ALLOCATE(f(1:np))
      ALLOCATE(mu(1:np))
      ALLOCATE(lapl(1:np))
      
      do tstep=1, nstep
         
         ! calculating laplacian of phi
         f(1:np) = field(1:np)%phi
         call dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)
         
         ! Chemical potential:  phi**3 - phi - epsilon*laplacian(phi )
         do ip = 1, np
            mu(ip) = field(ip)%phi * field(ip)%phi *field(ip)%phi &
                 - field(ip)%phi - interface_width*lapl(ip)			 
         end do
         
         ! Calculating laplacian of mu 
         
         call dderivatives_lapl(mu, lapl, np, dr, lxyz, lxyz_inv)

         do ip=1, np
            field(ip)%phi = field(ip)%phi + dt*lapl(ip)
         end do
         
      end do
      
      OPEN (UNIT=50,FILE=dir_name//'/nec.out')
      do ip=1, np
         if(field(ip)%phi>=0.d0) then 	
            write(50,'(I10,I10,I10,F10.2)') lxyz(ip,1:3), field(ip)%phi 	
         end if
      end do
      close(50)	 
      
      DEALLOCATE(f)
      DEALLOCATE(mu)
      DEALLOCATE(lapl)	
      
    end subroutine cahn_hilliard
    
    subroutine etc_search(cell, tipc, gg, Lsize, lxyz, lxyz_inv, notch_distance, n_tipcell, tip_s, &
         vegf_c, nstep, vegf_grad_min, np, cell_radius, n_max_tipc, np_tip_s, periodic)
      
      implicit none
      ! external variables
      type(mesh_t),allocatable,  intent(in) :: cell(:)
      type(tip_cell_t),allocatable,  intent(inout) :: tipc(:)
      real,allocatable, intent(in) :: gg(:,:)
      integer, intent(in) :: np, n_max_tipc
      integer,allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      real, intent(in) :: vegf_c, notch_distance, vegf_grad_min, cell_radius
      integer, intent(in) :: nstep, Lsize(1:3), np_tip_s
      integer,allocatable, intent(in) ::  tip_s(:,:)
      integer, intent(inout) ::  n_tipcell
      logical, intent(in) :: periodic
      ! internal
      integer :: i, j, k, l, dx, dy, dz, dxyz, ip, ip2
      real :: grad_T, hs, temp, temp_phi,  M_Pi
      logical :: signal, activated

      !np_tip_s = 25
      M_Pi = 3.14159265359

      activated = .false.

	  
      do ip=1, np
         
		 
         if(cell(ip)%phi>0.9) then ! phi=1
            
            if(cell(ip)%T>vegf_c) then ! T>Tc
               
               grad_T = gg(ip,1)*gg(ip,1) + gg(ip,2)*gg(ip,2) + gg(ip,3)*gg(ip,3)
               
               if(grad_T>vegf_grad_min) then ! grad(T) > Gm
         
                  signal = .true.
                  if(n_tipcell>0) then
                     do l = 1, n_tipcell
                        
                        dx = lxyz(ip,1) - tipc(l)%x
                        dy = lxyz(ip,2) - tipc(l)%y
                        dz = lxyz(ip,3) - tipc(l)%z
                        
                        dxyz = dx*dx + dy*dy + dz*dz
                        
                        if (dxyz<notch_distance) then
                           signal = .false.
                           EXIT
                        end if
                     end do
                  end if
                  
                  if(signal) then
                     
                     !  all cell borders have phi>0 ?
                     temp = 0.d0
                     
                     do l = 1, np_tip_s ! number of points given by spherical surface
                        
                        i = lxyz(ip,1)  + tip_s(l,1)
                        j = lxyz(ip,2)  + tip_s(l,2)
                        k = lxyz(ip,3)  + tip_s(l,3)
                                                

                        ! boundary condition in z-axis for non periodic system
                        temp = k
                        hs = heaviside(temp)

                        if(.not.periodic) then
                           if(abs(k)>Lsize(3)-hs -int(1.5*cell_radius)) then
                              signal = .false.
                              EXIT
                           end if
                        end if
                        
                        temp = 0.d0 
                        ip2 = lxyz_inv(i,j,k) ! new ip global point
                        
                        
                        if(cell(ip2)%phi<0.9) then
                           signal = .false.
                           EXIT
                        end if ! spherical verification phi=1
                     end do
                     
                     if (signal) then
                        
                        ! comparing phi between actual candidate and
                        ! the last selected
                        
                        if ( cell(ip)%phi>temp_phi ) then
                           
                           tipc(n_tipcell+1)%x = lxyz(ip,1)
                           tipc(n_tipcell+1)%y = lxyz(ip,2)
                           tipc(n_tipcell+1)%z = lxyz(ip,3)
                           
                           tipc(n_tipcell+1)%phi = cell(ip)%phi
                              
                           tipc(n_tipcell+1)%ip = ip
                           
                           temp_phi = cell(ip)%phi
                           
                           activated = .true.
                           
                        end if ! tip_phi(np_local) > tip_phi(np_local -1)
                        
                     end if  ! all cell radius with phi = 1
                  end if ! signal
               end if  ! grad(T) > Gm
            end if  ! T>Tc
         end if  ! phi = 1.0
      end do  ! activating tip cell
      
      if(activated) then
         n_tipcell = n_tipcell + 1
         write(*,'(A,I10,I10)') "activated - n_tip, nstep:", n_tipcell, nstep
         write(*,'(I10,I10,I10)') lxyz(tipc(n_tipcell)%ip,1:3)
         if(n_tipcell.eq.n_max_tipc) then
            write(*,'(A)') "WARNING: Reached the maximum number of ETCs. The program will be terminated to avoid SEGFAULT."
            write(*,'(A)') "The last state of the system will be saved."
         end if
      end if


    end subroutine etc_search


    subroutine etc_move(cell, tipc, gg, lxyz, lxyz_inv, Lsize, vegf_c, vegf_grad_min, vegf_grad_max, tip_s, &
               tip_all, n_tipcell, phi_max, nstep, dt, chi, cell_radius, np_tip_all, grid_cell_domain, n_source, &
			   vegf_xyz, ndim, np_tip_s, np, periodic)

      implicit none

      ! external variables
      type(mesh_t), allocatable, intent(inout) :: cell(:)
      type(tip_cell_t), allocatable, intent(inout) :: tipc(:)
      real, allocatable, intent(in) :: gg(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), vegf_xyz(:,:),  grid_cell_domain(:)
      real, intent(in) :: vegf_c, vegf_grad_min, phi_max, vegf_grad_max, dt, chi, cell_radius
      integer, intent(in) :: nstep, Lsize(1:3), np_tip_all, n_source, ndim, np_tip_s, np
      integer, allocatable, intent(in) ::  tip_s(:,:), tip_all(:,:)
      integer, intent(inout) ::  n_tipcell
      logical, intent(in) :: periodic
      ! internal
      integer :: jdim
      real, allocatable :: phi_boolean(:), tip_points(:)
      real :: grad_T, hs, temp_phi, temp, dx_vec(3), relative_pos(1:3), GG_projection, GG_rejection(1:3), grad(1:3),&
           relative_image(3), unitary(3), abs_distance
      integer :: i, j, l, k, m, ip, ip2, sinal 

      unitary(:)=1
      ALLOCATE(phi_boolean(1:np))
      ALLOCATE(tip_points(1:np))
      ! Deactivating tip cells
      phi_boolean(:) = cell(:)%phi
      tip_points(:) = -1.d0

      ! deactivating tip cells
      
      do i=1, n_tipcell
         ip = tipc(i)%ip
         grad_T = gg(ip,1)*gg(ip,1) + gg(ip,2)*gg(ip,2) + gg(ip,3)*gg(ip,3)
         
         if(cell(ip)%T<vegf_c.or.grad_T<vegf_grad_min) then
            write(*,'(A,I10,I10)') "deactivated - n_tip, nstep:", n_tipcell-1, nstep
            write(*,'(I10,I10,I10)') lxyz(ip,1:3) 
            tipc(i)%ip = tipc(n_tipcell)%ip
            tipc(i)%x = tipc(n_tipcell)%x
            tipc(i)%y = tipc(n_tipcell)%y
            tipc(i)%z = tipc(n_tipcell)%z
            tipc(i)%phi = tipc(n_tipcell)%phi
            n_tipcell = n_tipcell - 1
         end if
      end do

      do ip2 = 1, n_tipcell
         
         
         ! Integrating tipcells velocity
         ! Euler Method
         
         ip = tipc(ip2)%ip
         
         ! avoiding ETCs and Hypoxic Cells superposition
         grad(1:3) = gg(ip,1:3)
      
         if(grid_cell_domain(ip)> 0) then
            
            relative_pos(1:3) =  vegf_xyz(grid_cell_domain(ip),1:3)  - lxyz(ip,1:3)
            relative_image(1:ndim) = 2.d0*Lsize(1:ndim) -unitary(1:ndim)- abs(relative_pos(1:ndim))        

            do jdim=1, ndim
               if (abs(relative_pos(jdim)).gt.abs(relative_image(jdim))) then
                  relative_pos(jdim) = relative_image(jdim)
               end if
            end do

            abs_distance = sqrt(relative_pos(1)**2+ relative_pos(2)**2  + relative_pos(3)**2)

            relative_pos(1:3) = relative_pos(1:3)/abs_distance

            GG_projection =  grad(1)*relative_pos(1) + grad(2)*relative_pos(2) + grad(3)*relative_pos(3)     
     
            relative_pos(1:3) = relative_pos(1:3)*GG_projection

            GG_rejection(1:3) = grad(1:3) - relative_pos(1:3)

            grad(1:3) = GG_rejection(1:3)

         end if
        
    
         grad_T = grad(1)*grad(1) + grad(2)*grad(2) + grad(3)*grad(3)
         grad_T = grad_T**(0.5)


         if(grad_T<vegf_grad_min) grad_T = vegf_grad_min

         hs = heaviside(grad_T - vegf_grad_max)
         temp = 1.0 + (vegf_grad_max/grad_T - 1.0)*hs

         dx_vec(1) = chi * grad(1)*temp*dt ! dx_i (i-component of increment position)
         dx_vec(2) = chi * grad(2)*temp*dt
         dx_vec(3) = chi * grad(3)*temp*dt
        
         
         tipc(ip2)%x = tipc(ip2)%x + dx_vec(1)
         tipc(ip2)%y = tipc(ip2)%y + dx_vec(2)
         tipc(ip2)%z = tipc(ip2)%z + dx_vec(3)
        
         ! boundary condiditions
         
         if(tipc(ip2)%x>Lsize(1) - 1 ) then
            tipc(ip2)%x = tipc(ip2)%x - 2.d0*Lsize(1) + 1.d0
         else if(tipc(ip2)%x< -Lsize(1)) then
            tipc(ip2)%x = tipc(ip2)%x + 2.d0*Lsize(1) - 1.d0
         end if
         
         if(tipc(ip2)%y>Lsize(2) - 1 ) then
            tipc(ip2)%y = tipc(ip2)%y - 2.d0*Lsize(2) + 1.d0
         else if(tipc(ip2)%y< -Lsize(2)) then
            tipc(ip2)%y = tipc(ip2)%y + 2.d0*Lsize(2) - 1.d0
         end if
         
         if(periodic) then
            if(tipc(ip2)%z>Lsize(3) - 1 ) then
               tipc(ip2)%z = tipc(ip2)%z - 2.d0*Lsize(3) + 1.d0
            else if(tipc(ip2)%z< -Lsize(3)) then
               tipc(ip2)%z = tipc(ip2)%z + 2.d0*Lsize(3) - 1.d0
            end if
         else
            if(tipc(ip2)%z>=Lsize(3)-1-1.5*cell_radius .or. tipc(ip2)%z<=-Lsize(3)+1.5*cell_radius) then
               tipc(ip2)%z = tipc(ip2)%z - dx_vec(3)
             end if
         end if


         tipc(ip2)%ip = lxyz_inv(int(anint(tipc(ip2)%x)),&
              int( anint(tipc(ip2)%y)),&
              int( anint(tipc(ip2)%z))) ! new ip global

         ! Calculating Phi_c inside tipcell
           
         tipc(ip2)%phi =  (cell(int(tipc(ip2)%ip))%alpha_p*cell_radius)/&
              (2.d8*chi*grad_T*abs(temp)) ! phi_c


         if(phi_boolean(int(tipc(ip2)%ip))<0) then
            cell(int( tipc(ip2)%ip ))%phi = tipc(ip2)%phi
         end if
           
         do l=1, np_tip_all
            
            i = int(anint(tipc(ip2)%x  + tip_all(l,1)))
            j = int(anint(tipc(ip2)%y  + tip_all(l,2)))
            k = int(anint(tipc(ip2)%z  + tip_all(l,3)))
            
            m = lxyz_inv(i,j,k) ! new ip global point
           
            if(phi_boolean(m)<0) then 
      
               if(cell(m)%phi<phi_max) then  
                  cell(m)%phi = tipc(ip2)%phi
               else 
                  cell(m)%phi = phi_max
               end if
      
            end if
            
            tip_points(m) = 1.0
            
         end do

        
      end do
      
      do ip=1, np
         if(tip_points(ip)>0) then
            cell(ip)%T = cell(ip)%T +  dt*cell(ip)%consumption
         end if
      end do

      DEALLOCATE(phi_boolean)    
      DEALLOCATE(tip_points)

    end subroutine etc_move
    
    
    
    subroutine source_deactivate(cell, vegf_xyz, n_source, vegf_s, lxyz, lxyz_inv, np_vegf_s, Lsize, periodic)
      
      implicit none
      
      type(mesh_t), allocatable, intent(inout) :: cell(:)
      integer, allocatable, intent(inout) :: vegf_xyz(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), vegf_s(:,:)
      integer, intent(in) :: Lsize(3), np_vegf_s
      integer, intent(inout) :: n_source
      logical, intent(in) :: periodic
      ! internal
      integer :: n_source_o, sinal, r(3), i, j, ip, ip_source, temp(3)
      logical :: sair

  !    integer :: remove_list
      ! deactivating vegf sources

      vegf_copy(:,:) = vegf_xyz(:,:)

      n_source_o = n_source
      
      sair = .false.
      do i=1, n_source_o
         
         ip_source = lxyz_inv(vegf_xyz(i,1),vegf_xyz(i,2),vegf_xyz(i,3))
         
         do j=1, np_vegf_s
            
            r(1) = vegf_xyz(i,1) + vegf_s(j,1)
            r(2) = vegf_xyz(i,2) + vegf_s(j,2)
            r(3) = vegf_xyz(i,3) + vegf_s(j,3)
            
            ip = lxyz_inv(r(1),r(2),r(3))
            
            if( cell(ip)%phi>0.0) then
               
               cell(ip_source)%source = -1

               temp(1:3) = vegf_xyz(i,1:3)
               vegf_xyz(i,1:3) = vegf_xyz(n_source,1:3)
               vegf_xyz(n_source,1:3) = temp(1:3)

               n_source = n_source - 1
               sair = .TRUE.
               EXIT
            end if

         end do

         if(sair) EXIT

      end do

      
    end subroutine source_deactivate
    
    subroutine dderivatives_grad(f, gradient, np, lxyz, lxyz_inv, dr)
      
      implicit none
      
      type(mesh_t), allocatable, intent(in) :: f(:)
      real, allocatable, intent(inout) :: gradient(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) :: np,  dr(3)
      integer :: ip, i, j, k
  
      do ip=1, np

         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)
         
         
         gradient(ip,1) = (f( lxyz_inv(i+dr(1),j,k) )%T -&
              f( lxyz_inv(i-dr(1),j,k) )%T )/(2.d0*dr(1))
         
         gradient(ip,2) = (f( lxyz_inv(i,j+dr(2),k) )%T -&
              f( lxyz_inv(i,j-dr(2),k) )%T )/(2.d0*dr(2))
         
         gradient(ip,3) = (f( lxyz_inv(i,j,k+dr(3)) )%T  -&
              f( lxyz_inv(i,j,k-dr(3)) )%T )/(2.d0*dr(3))
         
      end do
      
    end subroutine dderivatives_grad
    
    
    subroutine dderivatives_lapl(f, lapl, np, dr, lxyz, lxyz_inv)
      
      implicit none

      real, allocatable, intent(inout) :: f(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      real, allocatable, intent(inout) :: lapl(:)
      integer, intent(in) :: np, dr(3)
      integer :: i, j, k, ip
      
      do ip=1, np
         
         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)
         
         lapl(ip) = (f(lxyz_inv(i+dr(1),j,k)) + f( lxyz_inv(i-dr(1),j,k))-2.d0*f(ip))/(real(dr(1)**2)) + &
              (f(lxyz_inv(i,j+dr(2),k)) + f( lxyz_inv(i,j-dr(2),k))-2.d0*f(ip))/(real(dr(2)**2)) + &
              (f(lxyz_inv(i,j,k+dr(3))) + f( lxyz_inv(i,j,k-dr(3)))-2.d0*f(ip))/(real(dr(3)**2))
         
      end do
      
      
    end subroutine dderivatives_lapl
    
    
    subroutine simul_box_init(cell, necrotic_tissue, necrotic_all, cell_radius, lxyz, lxyz_inv, Lsize,&
         vegf_xyz, vegf_source_conc, source_max, n_source, np, diff_oxy_length, vessel_radius, np_necrotic_all, ndim, iseed)
      
      implicit none
      
      ! input variables
      type(mesh_t), allocatable, intent(inout) :: cell(:)
      type(sec_mesh_t), allocatable, intent(inout) :: necrotic_tissue(:)
      integer, allocatable, intent(inout) :: vegf_xyz(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(in) :: source_max, np, iseed, np_necrotic_all, necrotic_all(:,:), Lsize(3), ndim
      real, intent(in) :: vegf_source_conc, diff_oxy_length, vessel_radius, cell_radius
      ! output
      integer, intent(inout) :: n_source
      ! internal only
      integer :: ip_source, ip, i, j, k, l, m, sinal, i_nec, j_nec, k_nec, attempt
      real :: r_n, dr, ij_radius, frontier_radius, distance(1:3), distance_abs
      logical :: superposition
      
      frontier_radius = real(min(Lsize(1),Lsize(2) )) - 2.d0*cell_radius

      ! initializing vessel and ECM
      do ip = 1, np
         
         i = lxyz(ip,1)
         j = lxyz(ip,2)
         k = lxyz(ip,3)
         
         necrotic_tissue(ip)%phi = -1.d0

         ij_radius = sqrt( real((i)**2 + (j)**2))

         if ( ij_radius <= vessel_radius ) then
            cell(ip)%phi = +1.d0
            cell(ip)%T = 0.d0
            cell(ip)%source = 0
         else
            cell(ip)%phi = -1.d0   
            cell(ip)%source = 0
            cell(ip)%T = 0.d0
         end if

         ! initiating the VEGF distribution on the ECM outside the sources
         ! T(R) = Ts*exp(Ri(R-d-Rc)/10.0)

         if( ij_radius>vessel_radius .and.  ij_radius <= vessel_radius+diff_oxy_length) then
            cell(ip)% T = vegf_source_conc*exp(0.1 *  (real(ij_radius)-25.0)  )
         end if

         if(ij_radius> vessel_radius+diff_oxy_length) then
            cell(ip)% T = vegf_source_conc
         end if
         
      end do



      n_source = 0
      attempt = 0
      ! random distribution
      do while (attempt<1000) 

         attempt = attempt + 1
         
         do ip_source=1, source_max
            
            r_n = ran2(iseed)
            ip = int(anint(r_n*np))
            
            i = lxyz(ip,1)
            j = lxyz(ip,2)
            k = lxyz(ip,3)  

            ij_radius = sqrt( real((i)**2 + (j)**2))

            if ( ij_radius > vessel_radius+diff_oxy_length &
                 .and. abs(k)<Lsize(3)-cell_radius .and. &
                 ij_radius < frontier_radius ) then
               
               ! verifying superposition 
               superposition = .false.
               
               do l=1, n_source

                  ! minimum image method

                  distance(1:3) = abs(vegf_xyz(l,1:3) - lxyz(ip,1:3))
                  distance(1:ndim) = min(distance(1:ndim), 2*Lsize(1:ndim)-1-abs(distance(1:ndim)))
                  distance = distance*distance
                  distance_abs = sqrt(sum(distance))

                  !dr = sqrt(real((i-vegf_xyz(l,1))**2 + (j-vegf_xyz(l,2))**2 + (k-vegf_xyz(l,3))**2))
                  if(distance_abs<4.d0*cell_radius+1.d0) then
                     superposition = .true.
                     EXIT
                  end if
               end do
               
               if(.not.superposition) then 
                  
                  ! creating necrotic tissue in the secondary mesh
                  
                  
                  do l=1, np_necrotic_all
                     
                     i_nec = i  + necrotic_all(l,1)
                     j_nec = j  + necrotic_all(l,2)
                     k_nec = k  + necrotic_all(l,3)
                     
                     m = lxyz_inv(i_nec,j_nec,k_nec) ! new ip global point
                     
                     necrotic_tissue(m)%phi = 1.d0
                     cell(m)%T = vegf_source_conc
                     
                  end do
                  
                  cell(ip)%T = vegf_source_conc
                  
                  cell(ip)%source = 1
                  
                  n_source = n_source + 1	
                  
                  vegf_xyz(n_source,1) = i
                  vegf_xyz(n_source,2) = j
                  vegf_xyz(n_source,3) = k
                  
                  
               end if
               
            end if
            
            
         end do
      end do
      

    end subroutine simul_box_init
  
    
    subroutine space_init(Lsize, lxyz, lxyz_inv, boundary_points, np, ndim, periodic)
      
      implicit none
      
      ! input/output variables
      integer, intent(in) ::  Lsize(1:3), boundary_points
      integer, allocatable, intent(inout) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(out) :: ndim
      ! internal variables
      integer :: i, j, k, l, m, n, ip, ip_part, np
      real :: hs(1:3)
      logical :: boundary, periodic

      ip = 0
      ip_part = np

      ! bulk points 
      ! allocated from 1 to np

      do i=-Lsize(1), Lsize(1)-1
         do j=-Lsize(2), Lsize(2)-1
            do k=-Lsize(3), Lsize(3)-1
               
               ip = ip + 1
               
               lxyz(ip,1) = i
               lxyz(ip,2) = j
               lxyz(ip,3) = k
               
               lxyz_inv(i,j,k) = ip                 
               
            end do
         end do
      end do
     
      ! boundary points
      ! allocated from np to np_part

      do i=-Lsize(1)-boundary_points, Lsize(1)-1+boundary_points
         do j=-Lsize(2)-boundary_points, Lsize(2)-1+boundary_points
            do k=-Lsize(3)-boundary_points, Lsize(3)-1+boundary_points
               
               l = i
               m = j
               n = k

               boundary = .false.

               hs(1) = heaviside(real(i))
               hs(2) = heaviside(real(j))
               hs(3) = heaviside(real(k))

               if( abs(i)>Lsize(1)-hs(1)) then
                  boundary = .true.               
                  l = i - SIGN(1,i)*(2*Lsize(1)-hs(1))
               end if

               if( abs(j)>Lsize(2)-hs(2)) then
                  boundary = .true.               
                  m = j  - SIGN(1,j)*(2*Lsize(2)-hs(2))
               end if

               if( abs(k)>Lsize(3)-hs(3)) then
                  boundary = .true.
                  if(periodic) then
                     n = k - SIGN(1,k)*(2*Lsize(3)-hs(3))
                  else
                     n = k - SIGN(1,k)
                  end if
               end if
               
               if(boundary) then
                  ip_part = ip_part + 1
                  
                  lxyz(ip_part,1) = l
                  lxyz(ip_part,2) = m
                  lxyz(ip_part,3) = n
                  
                  lxyz_inv(i,j,k) = lxyz_inv(l, m, n)
               end if
               
               if(periodic) then
                  ndim = 3
               else
                  ndim = 2
               end if
               
               ! the updating of boundaries should be
               ! cell(ip_part) = cell( lxyz_inv( lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ) )
               ! because: 
               ! lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ->(l, m, n) in the bulk
               ! then lxyz_inv(m,n,l) will give the updated value in the bulk
               ! for the respective boundary position cell(ip_part)
               
            end do
         end do
      end do
      
    end subroutine space_init
    
    
    subroutine parameters_init(cell_radius, diffusion_const, interface_width, vegf_p, vegf_c, diff_oxy_length,&
         vegf_rate, vegf_source_conc, prolif_rate, vessel_radius, tstep, dt, chi, Lsize, dr, dir_name, iseed,&
         boundary_points, source_max, vegf_grad_min, vegf_grad_max, depletion_weight, output_period, extra_steps, &
         n_max_tipc, thinning,periodic)
      
      implicit none 
      
      real, intent(inout) :: cell_radius, diffusion_const, interface_width, vegf_p, vegf_c, diff_oxy_length, vegf_rate, &
           vegf_source_conc, prolif_rate, vessel_radius, dt, chi, vegf_grad_min, vegf_grad_max, depletion_weight
      integer, intent(inout) :: tstep, Lsize(3), iseed, boundary_points, source_max, dr(3), output_period, n_max_tipc, &
           extra_steps
      character(len=3), intent(inout) :: dir_name
      character(len=255) :: temp
      logical :: periodic, thinning
      
      OPEN (UNIT=1,FILE='input_file')
      read(1,*) cell_radius, temp ! R_c - Cell Radius
      read(1,*) diffusion_const, temp ! D - Ang. Fac. Diffusion Constant
      read(1,*) interface_width, temp ! Eps - Interface Witdh
      read(1,*) vegf_p, temp ! T_p - Ang. fac. conc. for highest proliferation
      read(1,*) vegf_c, temp ! T_c - Critical Ang. Fac. for branching
      read(1,*) diff_oxy_length, temp ! diff_oxy - Min. Oxygen Diffusion Radius
      read(1,*) vegf_rate, temp ! alpha_T - Ang. Fac. Consumption Rate
      read(1,*) vegf_source_conc, temp ! Ts - Ang. Fac. on the source
      read(1,*) prolif_rate, temp ! Alphap - Proliferation rate
      read(1,*) vessel_radius, temp ! radius - Initial Vessel Radius
      read(1,*) tstep, temp ! tstep - Total time step
      read(1,*) dt, temp ! dt - Time increment
      read(1,*) chi, temp ! chi - Chemotaxis
      read(1,*) Lsize(1:3), dr(1:3), temp !Box Length - x,y,z, dr
      read(1,*) dir_name, temp ! Simulation name 
      read(1,*) iseed, temp ! Initial Seed for RAN2
      read(1,*) boundary_points, temp ! boundary points
      read(1,*) source_max, temp ! max number of vegf sources
      read(1,*) vegf_grad_min, temp ! minimum vegf gadient for branching
      read(1,*) vegf_grad_max, temp ! vegf for max velocity
      read(1,*) depletion_weight, temp ! repulsive force coefficient
      read(1,*) output_period, temp ! the period which the data will be written to the output
      read(1,*) extra_steps, temp ! number of extra steps after the desactivation of all hypoxic cells
      read(1,*) n_max_tipc, temp ! max number of tip cells
      read(1,*) thinning, temp ! thinning on/off
      read(1,*) periodic, temp ! boundary conditions
      CLOSE(1)
      call system('mkdir '//dir_name)

      OPEN (UNIT=2,FILE=dir_name//'/parameters'//dir_name//'.init')
      write(2,'(F10.2,A)') cell_radius, " cell_radius" ! R_c - Cell Radius
      write(2,'(F10.2,A)') diffusion_const, " diffusion_const" ! D - Ang. Fac. Diffusion Constant
      write(2,'(F10.2,A)') interface_width, " interface_width" ! Eps - Interface Witdh
      write(2,'(F10.2,A)') vegf_p, " vegf_p" ! T_p - Ang. fac. conc. for highest proliferation
      write(2,'(F10.3,A)') vegf_c,  " vegf_c" ! T_c - Critical Ang. Fac. for branching
      write(2,'(F10.2,A)') diff_oxy_length, " diff_oxy_length" ! diff_oxy - Min. Oxygen Diffusion Radius
      write(2,'(F10.2,A)') vegf_rate, " alpha_T" ! alpha_T - Ang. Fac. Consumption Rate
      write(2,'(F10.2,A)') vegf_source_conc, " vegf_source_conc" ! Ts - Ang. Fac. on the source
      write(2,'(F10.2,A)') prolif_rate, " prolif_rate" ! Alphap - Proliferation rate
      write(2,'(F10.2,A)') vessel_radius, " vessel_radius" ! radius - Initial Vessel Radius
      write(2,'(I10,A)')   tstep, " tstep" ! tstep - Total time step
      write(2,'(F10.5,A)') dt, " dt" ! dt - Time increment
      write(2,'(F10.2,A)') chi, " chi" ! chi - Chemotaxis
      write(2,'(I3, I3, I3, I3, I3, I3,A)') Lsize(1:3), dr(1:3), " box_length_xyz_dr" !Box Length - x,y,z, dr
      write(2,'(A,A)') dir_name, " dir_name" ! Simulation name 
      write(2,'(I10,A)') iseed, " iseed" ! Initial Seed for RAN2
      write(2,'(I10,A)') boundary_points, " bounary_points"! boundary points
      write(2,'(I10,A)') source_max, " source_max" ! max number of vegf sources
      write(2,'(F10.2,A)') vegf_grad_min, " vegf_grad_min" ! minimum vegf gadient for branching
      write(2,'(F10.2,A)') vegf_grad_max, " vegf_grad_max" ! vegf for max velocity
      write(2,'(F10.2,A)') depletion_weight, " depletion_weight" ! repulsive force coefficient
      write(2,'(I10,A)') output_period, " output_period" ! the period which the data will be written to the output
      write(2,'(I10,A)') extra_steps, " extra_steps"
      write(2,'(I10,A)') n_max_tipc, " n_max_tipc" ! max number of tip cells
      write(2,'(L1,A)') thinning, " thinning" ! thinning on/off
      write(2,'(L1,A) ') periodic, " periodic" ! boundary conditions
      CLOSE(2)      
    end subroutine parameters_init
    

    
    subroutine print_header(Lsize, n_source, dir_name, periodic)

      implicit none
      integer, intent(in) :: Lsize(1:3), n_source
      character(len=255) :: cwd, hostname
      character(len=32) :: username
      character(8)  :: date
      character(10) :: time
      character(5)  :: zone
      character(3)  :: dir_name
      integer, dimension(8) :: values
      logical :: periodic

      call date_and_time(date,time,zone,values)
      call date_and_time(DATE=date,ZONE=zone)
      call date_and_time(TIME=time)
 
      call hostnm(hostname) 
      call getcwd(cwd)
      call getlog(username) 
      write(*,'(A)') "                                Running Angio"
      write(*,'(A)') "       "
      write(*,'(A)') "Version        :       4.1.s"
      write(*,'(A,A)') "Locate         :       ", trim(cwd)
      write(*,'(A,A)') "User           :       ", trim(username)
      write(*,'(A)') "Developer      :       Moreira, M."
      write(*,'(A)') "       "
      write(*,'(A,A)') "                      The code is running in ", trim(hostname)
      write(*,'(A)') "       "
      write(*,'(A)') "       "
      write(*,'(A,2X,A,2X,A)') "             Calculation started on", date(7:8)//"/"//date(5:6)//"/"//date(1:4),&
           time(1:2)//":"//time(3:4)//":"//time(5:6)
      write(*,'(A)') "       "
      write(*,'(A)') "************************************ Grid *************************************"
      write(*,'(A)') "Simulation Box:"
      write(*,'(A,I3,A,I3,A,I3,A)') "  Lengths = (",Lsize(1),",", Lsize(2),",", Lsize(3), ")"
      write(*,'(A)') "  the code will run in 3 dimension(s)."
      if(periodic) then
         write(*,'(A)') "  the code will treat the system as periodic in 3 dimension(s)."
      else
         write(*,'(A)') "  the code will treat the system as periodic in 2 dimension(s)."
      end if
      write(*,'(A)') "*******************************************************************************"
      write(*,'(A,A)') "Simulation ID: ", dir_name
      write(*,'(A,I10)') "VEGF Sources:", n_source
    end subroutine print_header
    

    function heaviside(x)
      
      implicit none
      
      real :: x, heaviside
      
      if ( x<0.d0) then
         heaviside = 0.d0
      else
         heaviside = 1.d0
      end if
    end function heaviside


    subroutine spherical_surface(Rc, R, ndim, delta)
      
      implicit none
      
      real, intent(in) :: Rc, delta
      integer, intent(inout) :: ndim
      integer, allocatable, intent(inout) :: R(:,:)
      real :: theta, phi, M_Pi
      character(len=10) :: dr, ds, sdx, sdy, sdz
      integer :: i, dx, dy, dz
      M_Pi = 3.14159265359
      ! Note:
      ! No intrinsic exists to convert between a numeric value and a formatted character 
      ! string representation – for instance, given the CHARACTER value '154',
      ! obtaining an INTEGER or REAL value with the value 154, or vice versa.
      ! Instead, this functionality is provided by internal-file I/O, 
      ! as in the following example: 
    
      ! Convert a string to a numeric value
      ! read (string,'(I10)') value
      ! print *, value
          
      ! Convert a value to a formatted string
      ! write (string2,'(I10)') value
      ! print *, string2


      i = 0
      theta = 0.0
      phi = 0.0      
      
      do while (theta<=M_Pi) 
         phi = 0.0
         do while (phi<=2.d0*M_Pi)
            dx =  int(anint(Rc * sin(theta) * cos(phi)) )
            dy =  int(anint(Rc * sin(theta) * sin(phi) ))
            dz =  int(anint(Rc * cos(theta)) )
            write(sdx,'(I2)') abs(dx)
            write(sdy,'(I2)') abs(dy)
            write(sdz,'(I2)') abs(dz)

            dr = trim(sdx)//trim(sdy)//trim(sdz)

            if (dr .ne. ds) then
               i = i + 1
              
               R(i,1) = dx 
               R(i,2) = dy  
               R(i,3) = dz  

               ds = dr
            end if
            phi = phi + delta
         end do
         theta = theta + delta
      end do
 
      i = i +1
	  R(i,1) = 0
	  R(i,2) = 0
	  R(i,3) = -Rc
      ndim = i

    end subroutine spherical_surface
    


    subroutine gen_cell_points(Rc,R,ndim)
    
      implicit none
      
      real, intent(in) :: Rc
      integer, intent(inout) :: ndim
      integer, allocatable, intent(inout) :: R(:,:)
      real :: dx, i, j, k,  M_Pi
      integer :: counter
      M_Pi = 3.14159265359
      counter = 0
      i = -Rc
      do while(i<=Rc)
         j =-Rc
         do while (j<=Rc)
            k =-Rc
            do while(k<=Rc)
               dx = sqrt(i**2 + j**2+ k**2 )
               if(dx <= Rc) then
                  counter = counter + 1
                  R(counter,1) = int(anint(i))
                  R(counter,2) = int(anint(j))
                  R(counter,3) = int(anint(k))        
                  
               end if
               k = k + 1.d0
            end do
            j = j +1.d0
         end do
         i = i + 1.d0
      end do

      ndim = counter

    end subroutine gen_cell_points


    function ran2(idum) 
      integer :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV 
      real :: ran2,AM,EPS,RNMX 
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, & 
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, & 
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS) 
      integer :: idum2,j,k,iv(NTAB),iy 
      save iv,iy,idum2 
      data idum2/123456789/, iv/NTAB*0/, iy/0/ 
      if (idum.le.0) then 
         idum=max(-idum,1) 
         idum2=idum 
         do 11 j=NTAB+8,1,-1 
            
            k=idum/IQ1 
            idum=IA1*(idum-k*IQ1)-k*IR1 
            if (idum.lt.0) idum=idum+IM1 
            if (j.le.NTAB) iv(j)=idum 
11          continue 
            iy=iv(1) 
         endif 
         k=idum/IQ1 
         idum=IA1*(idum-k*IQ1)-k*IR1 
         if (idum.lt.0) idum=idum+IM1 
         k=idum2/IQ2 
         idum2=IA2*(idum2-k*IQ2)-k*IR2 
         if (idum2.lt.0) idum2=idum2+IM2 
         j=1+iy/NDIV 
         iy=iv(j)-idum2 
         iv(j)=idum 
         if(iy.lt.1)iy=iy+IMM1 
         ran2=min(AM*iy,RNMX) 
         return 
       end function ran2 
       

  end module run_angio_m
  

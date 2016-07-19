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

  use global_m
  use thinning_m
  use init_m
  use dderivatives_m
  use etc_dynamics_m
  use source_m
  use blood_flow_m
  use misc_m
  use congraph_m

  implicit none

  private

  public :: run_angio


  contains

    subroutine run_angio()

      implicit none


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
      ! flow
      calculate_flow = .true.
      calc_flow_period = 50
      flow_count = 0

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
      ALLOCATE(graph(1:np))
      ! surfaces and points from spheres
      ALLOCATE(vegf_s(1:10000,1:3))
      ALLOCATE(tip_s(1:300,1:3))
      ALLOCATE(tip_all(1:500,1:3))
      ALLOCATE(flow(1:np))
      ALLOCATE(flow_full(1:np))

      ! initializing tip cell points
      call spherical_surface(cell_radius, tip_s, np_tip_s, 0.2)  ! spherical surface increments of tip cell
      call gen_cell_points(cell_radius, tip_all, np_tip_all)     ! all points of tip cell

      ! initializing points for the vegf sources
      ! initializing points to fill with blood flow
      ALLOCATE(d2sphere(1:3000))
      ALLOCATE(sphere(1:3000,1:3))
      call gen_cell_points(8.0, sphere, np_sphere, d2sphere)
      ! initializing points for the vegf sources
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
      ! file_id = '10000'
      !                     field, ...,    ..., ... , directory, ...
      ! call init_from_file(cell, lxyz_inv, np, 702,'tst' ,  file_id)

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
      calc_flow_period = output_period

      do while(nstep<=tstep)
         nstep = nstep + 1

         ! the program will run extra_steps steps after the last source deactivation
         if(n_source.eq.0) then
            tstep = nstep + extra_steps
            n_source = -1
         end if

         if(n_tipcell>0) then
            calc_flow_period = 50
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

         flow_count = flow_count + 1
 	 if(calculate_flow) then

	   if(flow_count.eq.calc_flow_period .or. counter+1.eq.output_period) then
                flow(:) = -1.d0
                flow_count = 0
                do ip=1, np
                   if(cell(ip)%phi>=0) then
                      phis(ip) = 1.d0
                   else
                      phis(ip) = -1.d0
                   end if
                end do

                call thinning_run(phis, lxyz, lxyz_inv, lsize, np)
                do i=-5,5
                  do j=1-5,5
                    if(phis(lxyz_inv(i,j,-Lsize(3))).gt.0) initial_node = lxyz_inv(i,j,-Lsize(3))
                  end do
                end do
                call verify_graph_connection(phis, graph, initial_node, lxyz, lxyz_inv)
                call flow_calc(graph, flow, lxyz, lxyz_inv, Lsize, np, nstep)
                call fill_vessels(flow_full, cell%phi, Lsize, lxyz, lxyz_inv, flow, d2sphere, sphere, np, np_sphere)
                call source_deactivate_flow(cell, vegf_xyz, n_source, vegf_s, lxyz, lxyz_inv,&
                     np_vegf_s, Lsize, periodic, flow_full, vegf_all, np_vegf_all)


	   end if
 	 else
               call source_deactivate(cell, vegf_xyz, n_source, vegf_s, lxyz, lxyz_inv, np_vegf_s, Lsize, periodic)
 	 endif

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

            write(file_name,'(I6)') nstep
            OPEN (UNIT=nstep,FILE=dir_name//'/phi'//trim(file_name)//'.xyz')
            OPEN (UNIT=nstep+1,FILE=dir_name//'/t'//trim(file_name)//'.xyz')
            if(thinning) OPEN (UNIT=nstep+2,FILE=dir_name//'/phis'//trim(file_name)//'.xyz')
            do ip=1, np

               if(phis(ip)>0.and.thinning) then
                  write(nstep+2,'(I10,I10,I10,F10.2,F10.4)') lxyz(ip,1:3), phis(ip), flow(ip)
               end if

               if(cell(ip)%phi>0) then

                 write(nstep,'(I10,I10,I10,F10.2,F10.4)') lxyz(ip,1:3), cell(ip)%phi, flow_full(ip)
               end if
               if(lxyz(ip,3).eq.0) then
                  write(nstep+1,'(I10,I10,F10.3)') lxyz(ip,1:2), cell(ip)%T
               end if
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
      DEALLOCATE(d2sphere)
      DEALLOCATE(sphere)
      DEALLOCATE(vegf_all)
      DEALLOCATE(flow)
      DEALLOCATE(flow_full)
      DEALLOCATE(graph)

    end subroutine run_angio

  end module run_angio_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

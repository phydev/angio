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

module init_m

  use global_m
  use misc_m


  implicit none

  private

  public :: simul_box_init, space_init, parameters_init, print_header, init_from_file

  contains


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
                  l = i - SIGN(1,i)*(2*Lsize(1))!- SIGN(1,i)*heaviside(-real(i))
               end if

               if( abs(j)>Lsize(2)-hs(2)) then
                  boundary = .true.
                  m = j  - SIGN(1,j)*(2*Lsize(2))!- SIGN(1,j)*heaviside(-real(j))
               end if

               if( abs(k)>Lsize(3)-hs(3)) then
                  boundary = .true.
                  if(periodic) then
                     n = k - SIGN(1,k)*(2*Lsize(3))!- SIGN(1,k)*heaviside(-real(k))
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
      write(2,'(F10.3,A)') vegf_grad_min, " vegf_grad_min" ! minimum vegf gadient for branching
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
      write(*,'(A)') "                                Running Angio with Blood Flow"
      write(*,'(A)') "       "
      write(*,'(A)') "Version        :       5.6.s (August 1, 2016)"
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
      cell(:)%phi = -1.d0
      do ip=1, np_phi
         read(100,*) r(1:3), phi_temp
         cell(lxyz_inv(r(1),r(2),r(3) ) )%phi = phi_temp
      end do
      close(100)
      !open(UNIT=200, FILE=dir_name//'/t'//trim(file_id)//'.xyz')
      !do ip=1, np
      !   read(200,*) r(1:3), cell(ip)%T
      !end do
      !close(200)

    end subroutine init_from_file

  end module init_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

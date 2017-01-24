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

module etc_dynamics_m

  use global_m
  use misc_m

  implicit none

  private

  public :: etc_search, etc_move

  contains


     subroutine etc_search(cell, tipc, gg, Lsize, lxyz, lxyz_inv, notch_distance, n_tipcell, tip_s, &
         vegf_c, nstep, vegf_grad_min, np, cell_radius, n_max_tipc, np_tip_s, grid_cell_domain, vegf_xyz, &
         ndim, periodic)

      implicit none
      ! external variables
      type(mesh_t), allocatable, intent(in) :: cell(:)
      type(tip_cell_t), allocatable, intent(inout) :: tipc(:)
      real, allocatable, intent(in) :: gg(:,:)
      integer, intent(in) :: np, n_max_tipc
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), grid_cell_domain(:), vegf_xyz(:,:)
      real, intent(in) :: vegf_c, notch_distance, vegf_grad_min, cell_radius
      integer, intent(in) :: nstep, Lsize(1:3), np_tip_s, ndim
      integer, allocatable, intent(in) ::  tip_s(:,:)
      integer, intent(inout) ::  n_tipcell
      logical, intent(in) :: periodic
      ! internal
      real :: relative_pos(1:3), unitary(1:3), relative_image(1:3), abs_distance
      integer :: i, j, k, l, dx, dy, dz, dxyz, ip, ip2, jdim
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
                        ! minimum image method
                        dx = min(abs(dx), 2*Lsize(1)-1-abs(dx))
                        dy = min(abs(dy), 2*Lsize(2)-1-abs(dy))
                        if(ndim.eq.3) dz = min(abs(dz), 2*Lsize(3)-1-abs(dz))

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
                           if(abs(k)>Lsize(3)-hs-int(2.0*cell_radius)) then
                              signal = .false.
                              EXIT
                           end if
                        end if

                        temp = 0.d0
                        ip2 = lxyz_inv(i,j,k) ! new ip global point


                        if(cell(ip2)%phi<0.0) then
                           signal = .false.
                           EXIT
                        end if ! spherical verification phi=1
                     end do

                     if (signal) then
                       ! avoiding ETCs and Hypoxic Cells superposition
                       abs_distance = 1000
                       if(grid_cell_domain(ip)> 0) then

                          relative_pos(1:3) =  vegf_xyz(grid_cell_domain(ip),1:3)  - lxyz(ip,1:3)
                          relative_image(1:ndim) = 2.d0*Lsize(1:ndim) -unitary(1:ndim)- abs(relative_pos(1:ndim))

                          do jdim=1, ndim
                             if (abs(relative_pos(jdim)).gt.abs(relative_image(jdim))) then
                                relative_pos(jdim) = relative_image(jdim)
                             end if
                          end do

                          abs_distance = sqrt(relative_pos(1)**2+ relative_pos(2)**2  + relative_pos(3)**2)
                        end if



                        ! comparing phi between actual candidate and
                        ! the last selected

                        if ( cell(ip)%phi>temp_phi .and. abs_distance.gt.2*cell_radius) then

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


         if(grad_T<sqrt(vegf_grad_min)) grad_T = sqrt(vegf_grad_min)

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
            if(tipc(ip2)%z>=Lsize(3)-1-2.0*cell_radius .or. tipc(ip2)%z<=-Lsize(3)+2.0*cell_radius) then
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


  end module etc_dynamics_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

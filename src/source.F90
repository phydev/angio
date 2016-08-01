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

module source_m

  use global_m
  use misc_m

  implicit none

  private

  public :: source_deactivate, fill_vessels, source_deactivate_flow

  contains


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
      logical :: deactivated !sair
      real :: cutoff_check

      ! deactivating vegf sources

      n_source_o = n_source

      !sair = .false.

      deactivated = .false.

      do i=1, n_source_o

         ip_source = lxyz_inv(vegf_xyz(i,1),vegf_xyz(i,2),vegf_xyz(i,3))

         do j=1, np_vegf_s

            r(1) = vegf_xyz(i,1) + vegf_s(j,1)
            r(2) = vegf_xyz(i,2) + vegf_s(j,2)
            r(3) = vegf_xyz(i,3) + vegf_s(j,3)

            ip = lxyz_inv(r(1),r(2),r(3))


            if( cell(ip)%phi>0.0) then

               deactivated = .true.
               cell(ip_source)%source = -1

               !temp(1:3) = vegf_xyz(i,1:3)
               !vegf_xyz(i,1:3) = vegf_xyz(n_source,1:3)
               !vegf_xyz(n_source,1:3) = temp(1:3)

               !n_source = n_source - 1
               !sair = .TRUE.
               !EXIT
            end if

         end do

         if(.not.deactivated) cell(ip_source)%source = 1

         !if(sair) EXIT

      end do


    end subroutine source_deactivate


    subroutine fill_vessels(flow_full, phi, Lsize, lxyz, lxyz_inv, flow, d2sphere, sphere, np, nps)

      implicit none

      integer, intent(in) :: np, nps, Lsize(3)
      real, intent(in) :: phi(1:np)
      real, allocatable, intent(inout) :: flow_full(:)
      real, allocatable, intent(in) :: flow(:)
      integer, allocatable, intent(in) :: sphere(:,:), d2sphere(:), lxyz(:,:), lxyz_inv(:,:,:)

      ! intern variables
	  real :: hs(1:3)
      integer :: ip, ips, d2temp,ip2, r(1:3)

      flow_full(:) = flow(:)
      do ip=1, np

         if(phi(ip) .ge. 0.d0 .and. flow(ip) .lt. 0.d0) then

            d2temp = 100000

            do ips = 1, nps

               r(1:3) = lxyz(ip,1:3) + sphere(ips,1:3)

               ip2 = lxyz_inv(r(1),r(2),r(3))

               if(flow(ip2) .gt. 0.d0 .and. d2sphere(ips) .lt. d2temp) then
                  d2temp = d2sphere(ips)
                  flow_full(ip) =  flow(ip2)
               end if

            end do
         end if
      end do

    end subroutine fill_vessels

    subroutine source_deactivate_flow(cell, vegf_xyz, n_source, vegf_s, lxyz, lxyz_inv,&
         np_vegf_s, Lsize, periodic, flow)!, vegf_all, np_vegf_all)

      implicit none

      type(mesh_t), allocatable, intent(inout) :: cell(:)
      real, allocatable,  intent(in) :: flow(:)
      integer, allocatable, intent(in) :: vegf_xyz(:,:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:), vegf_s(:,:)!, vegf_all(:,:)
      integer, intent(in) :: Lsize(3), np_vegf_s!, np_vegf_all
      integer, intent(inout) :: n_source
      logical, intent(in) :: periodic
      ! internal
      integer :: n_source_o, sinal, r(3), i, j, ip, ip_source, temp(3)
      logical :: deactivated !sair
      real :: cutoff_check

      ! deactivating vegf sources

      n_source_o = n_source

      !sair = .false.


      do i=1, n_source_o
         deactivated = .false.

         ip_source = lxyz_inv(vegf_xyz(i,1),vegf_xyz(i,2),vegf_xyz(i,3))

         do j=1, np_vegf_s

            r(1) = vegf_xyz(i,1) + vegf_s(j,1)
            r(2) = vegf_xyz(i,2) + vegf_s(j,2)
            r(3) = vegf_xyz(i,3) + vegf_s(j,3)

            ip = lxyz_inv(r(1),r(2),r(3))


            if(10.0*flow(ip)>0.d8) then

               deactivated = .true.
               cell(ip_source)%source = -1

               !write(*,*) "hypoxic cell deactivated (x,y,z), n_source:", lxyz(ip_source,1:3), n_source -1
               !write(*,*) "xyz,flow:", lxyz(ip,1:3), flow(ip)
               !write(*,*) sqrt(real( (lxyz(ip,1)-lxyz(ip_source,1))**2+ (lxyz(ip,2)-lxyz(ip_source,2))**2 + (lxyz(ip,3)-lxyz(ip_source,3))**2))
               !temp(1:3) = vegf_xyz(i,1:3)
               !vegf_xyz(i,1:3) = vegf_xyz(n_source,1:3)
               !vegf_xyz(n_source,1:3) = temp(1:3)

               !n_source = n_source - 1
               !sair = .TRUE.
               !EXIT
            end if

         end do

         if(.not.deactivated) then

            cell(ip_source)%T = 1.0
            cell(ip_source)%source = 1

            !do j=1, np_vegf_all

            !   r(1) = vegf_xyz(i,1) + vegf_all(j,1)
            !   r(2) = vegf_xyz(i,2) + vegf_all(j,2)
            !   r(3) = vegf_xyz(i,3) + vegf_all(j,3)

            !   ip = lxyz_inv(r(1),r(2),r(3))

            !   cell(ip)%T = 1.d0
            !   cell(ip)%source = 1.d0

            !end do
         else

            !do j=1, np_vegf_all

            !   r(1) = vegf_xyz(i,1) + vegf_all(j,1)
            !   r(2) = vegf_xyz(i,2) + vegf_all(j,2)
            !   r(3) = vegf_xyz(i,3) + vegf_all(j,3)

            !   ip = lxyz_inv(r(1),r(2),r(3))

            cell(ip)%source = -1.d0
            !end do

         end if

         !if(sair) EXIT

      end do


    end subroutine source_deactivate_flow


  end module source_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

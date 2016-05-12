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
  
  public :: source_deactivate

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
    

  end module source_m
 

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

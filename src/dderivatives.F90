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

module dderivatives_m

  use global_m


  implicit none

  private

  public :: dderivatives_grad, dderivatives_lapl

  contains


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

      real, intent(inout) :: f(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      real, intent(inout) :: lapl(:)
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


  end module dderivatives_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

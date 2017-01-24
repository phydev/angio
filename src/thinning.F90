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

module thinning_m

   implicit none

   private

   public :: thinning_run

 contains

   subroutine thinning(phis, direction, finish, lxyz, lxyz_inv, lsize, np)
     ! direction = 1 or -1 for x+-, 2 or -2 for y+- and 3 or -3 for z+-

      implicit none

      integer, intent(in):: np, lsize(1:3)
      integer, intent(in) :: direction
      real, allocatable, intent(inout) :: phis(:)
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      integer, intent(inout) :: finish
      integer :: x, y, z
      integer :: M(-1:1,-1:1,-1:1), temp
      real :: vijk(1:3),vijk_ref(1:3), vijk_temp(1:3)
      integer, allocatable :: remove(:)
      integer :: i, j, k, n, ip, i_ip, j_ip, k_ip, statusf(6), xpstatus(2), statusrf(6), xpstatusrf(2), rotation
      integer :: masks(1:6, -3:3,1:4,0:1)
      real :: pi, angles(-3:3), t
      real :: Rx(1:3,1:3), Ry(1:3,1:3), Rz(1:3,1:3), Rdeg(1:3,1:3) ! rotate matrix
      real :: RefSimmX(1:3,1:3), RefSimmY(1:3,1:3), RefSimmZ(1:3,1:3) ! reflection
      integer :: M1(-1:1,-1:1,-1:1), M1_rot(-1:1,-1:1,-1:1), M1_ref(-1:1,-1:1,-1:1), M1_rot_ref(-1:1,-1:1,-1:1), &
           M2(-1:1,-1:1,-1:1), M2_rot(-1:1,-1:1,-1:1), M2_ref(-1:1,-1:1,-1:1), M2_rot_ref(-1:1,-1:1,-1:1), &
           M3(-1:1,-1:1,-1:1), M3_rot(-1:1,-1:1,-1:1), M3_ref(-1:1,-1:1,-1:1), M3_rot_ref(-1:1,-1:1,-1:1), &
           M4(-1:1,-1:1,-1:1), M4_rot(-1:1,-1:1,-1:1), M4_ref(-1:1,-1:1,-1:1), M4_rot_ref(-1:1,-1:1,-1:1), &
           M5(-1:1,-1:1,-1:1), M5_rot(-1:1,-1:1,-1:1), M5_ref(-1:1,-1:1,-1:1), M5_rot_ref(-1:1,-1:1,-1:1), &
           M6(-1:1,-1:1,-1:1), M6_rot(-1:1,-1:1,-1:1), M6_ref(-1:1,-1:1,-1:1), M6_rot_ref(-1:1,-1:1,-1:1)

      ALLOCATE(remove(1:np))
      pi = 3.14159265359
      ! image point: +1
      ! background point:  -1
      ! don't care point:  3
      ! unless 1 should be +1: 2


      ! M1 mask !

      M1( -1 , -1 , -1 ) = 2
      M1( 0 , -1 , -1 ) = 2
      M1( 1 , -1 , -1 ) = 2
      M1( -1 , 0 , -1 ) = 2
      M1( 0 , 0 , -1 ) = 1
      M1( 1 , 0 , -1 ) = 2
      M1( -1 , 1 , -1 ) = 2
      M1( 0 , 1 , -1 ) = 2
      M1( 1 , 1 , -1 ) = 2
      M1( -1 , -1 , 0 ) = 2
      M1( 0 , -1 , 0 ) = 2
      M1( 1 , -1 , 0 ) = 2
      M1( -1 , 0 , 0 ) = 2
      M1( 0 , 0 , 0 ) = 1
      M1( 1 , 0 , 0 ) = 2
      M1( -1 , 1 , 0 ) = 2
      M1( 0 , 1 , 0 ) = 2
      M1( 1 , 1 , 0 ) = 2
      M1( -1 , -1 , 1 ) = -1
      M1( 0 , -1 , 1 ) = -1
      M1( 1 , -1 , 1 ) = -1
      M1( -1 , 0 , 1 ) = -1
      M1( 0 , 0 , 1 ) = -1
      M1( 1 , 0 , 1 ) = -1
      M1( -1 , 1 , 1 ) = -1
      M1( 0 , 1 , 1 ) = -1
      M1( 1 , 1 , 1 ) = -1


      ! M2 mask !


      M2( -1 , -1 , -1 ) = 3
      M2( 0 , -1 , -1 ) = 3
      M2( 1 , -1 , -1 ) = 3
      M2( -1 , 0 , -1 ) = 3
      M2( 0 , 0 , -1 ) = 1
      M2( 1 , 0 , -1 ) = 3
      M2( -1 , 1 , -1 ) = 3
      M2( 0 , 1 , -1 ) = 3
      M2( 1 , 1 , -1 ) = 3
      M2( -1 , -1 , 0 ) = 3
      M2( 0 , -1 , 0 ) = 3
      M2( 1 , -1 , 0 ) = 3
      M2( -1 , 0 , 0 ) = 1
      M2( 0 , 0 , 0 ) = 1
      M2( 1 , 0 , 0 ) = 3
      M2( -1 , 1 , 0 ) = 3
      M2( 0 , 1 , 0 ) = 3
      M2( 1 , 1 , 0 ) = 3
      M2( -1 , -1 , 1 ) = 3
      M2( 0 , -1 , 1 ) = -1
      M2( 1 , -1 , 1 ) = -1
      M2( -1 , 0 , 1 ) = 3
      M2( 0 , 0 , 1 ) = -1
      M2( 1 , 0 , 1 ) = -1
      M2( -1 , 1 , 1 ) = 3
      M2( 0 , 1 , 1 ) = -1
      M2( 1 , 1 , 1 ) = -1


      ! M3 mask !

      M3( -1 , -1 , -1 ) = 3
      M3( 0 , -1 , -1 ) = 3
      M3( 1 , -1 , -1 ) = 3
      M3( -1 , 0 , -1 ) = 3
      M3( 0 , 0 , -1 ) = 1
      M3( 1 , 0 , -1 ) = 3
      M3( -1 , 1 , -1 ) = 3
      M3( 0 , 1 , -1 ) = 3
      M3( 1 , 1 , -1 ) = 3
      M3( -1 , -1 , 0 ) = 3
      M3( 0 , -1 , 0 ) = 3
      M3( 1 , -1 , 0 ) = 3
      M3( 1 , 0 , 0 ) = 3
      M3( -1 , 1 , 0 ) = 3
      M3( 1 , 1 , 0 ) = 3
      M3( -1 , 0 , 0 ) = 1
      M3( 0 , 0 , 0 ) = 1
      M3( 0 , 1 , 0 ) = 1
      M3( 0 , -1 , 1 ) = -1
      M3( 1 , -1 , 1 ) = -1
      M3( 0 , 0 , 1 ) = -1
      M3( 1 , 0 , 1 ) = -1
      M3( -1 , -1 , 1 ) = 3
      M3( -1 , 0 , 1 ) = 3
      M3( -1 , 1 , 1 ) = 3
      M3( 0 , 1 , 1 ) = 3
      M3( 1 , 1 , 1 ) = 3


      ! M4 mask !



      M4( -1 , -1 , -1 ) = 3
      M4( 0 , -1 , -1 ) = 3
      M4( 1 , -1 , -1 ) = 3
      M4( -1 , 0 , -1 ) = 3
      M4( 0 , 0 , -1 ) = +1
      M4( 1 , 0 , -1 ) = 3
      M4( -1 , 1 , -1 ) = 3
      M4( 0 , 1 , -1 ) = 3
      M4( 1 , 1 , -1 ) = 3
      M4( -1 , -1 , 0 ) = 3
      M4( 0 , -1 , 0 ) = 3
      M4( 1 , -1 , 0 ) = 3
      M4( -1 , 0 , 0 ) = 3
      M4( 0 , 0 , 0 ) = +1
      M4( 1 , 0 , 0 ) = 3
      M4( -1 , 1 , 0 ) = +1
      M4( 0 , 1 , 0 ) = 3
      M4( 1 , 1 , 0 ) = 3
      M4( -1 , -1 , 1 ) = -1
      M4( 0 , -1 , 1 ) = -1
      M4( 1 , -1 , 1 ) = -1
      M4( -1 , 0 , 1 ) = -1
      M4( 0 , 0 , 1 ) = -1
      M4( 1 , 0 , 1 ) = -1
      M4( -1 , 1 , 1 ) = +1
      M4( 0 , 1 , 1 ) = -1
      M4( 1 , 1 , 1 ) = -1

      ! M5 mask !
      M5( 0 , 0 , 0 ) = 1
      M5( -1 , 0 , -1 ) = 1
      M5( -1 , -1 , -1 ) = 2
      M5( 0 , -1 , -1 ) = 2
      M5( -1 , 1 , -1 ) = 2
      M5( 0 , 1 , -1 ) = 2
      M5( -1 , -1 , 0 ) = 2
      M5( 0 , -1 , 0 ) = 2
      M5( -1 , 1 , 0 ) = 2
      M5( 0 , 1 , 0 ) = 2
      M5( -1 , 0 , 0 ) = 2
      M5( 1 , -1 , -1 ) = -1
      M5( 0 , 0 , -1 ) = -1
      M5( 1 , 0 , -1 ) = -1
      M5( 1 , 1 , -1 ) = -1
      M5( 1 , -1 , 0 ) = -1
      M5( 1 , 0 , 0 ) = -1
      M5( 1 , 1 , 0 ) = -1
      M5( -1 , -1 , 1 ) = -1
      M5( 0 , -1 , 1 ) = -1
      M5( 1 , -1 , 1 ) = -1
      M5( -1 , 0 , 1 ) = -1
      M5( 0 , 0 , 1 ) = -1
      M5( 1 , 0 , 1 ) = -1
      M5( -1 , 1 , 1 ) = -1
      M5( 0 , 1 , 1 ) = -1
      M5( 1 , 1 , 1 ) = -1

      ! M6 mask !

      M6( -1 , -1 , -1 ) = 3
      M6( 0 , -1 , -1 ) = -1
      M6( 1 , -1 , -1 ) = -1
      M6( -1 , 0 , -1 ) = +1
      M6( 0 , 0 , -1 ) = -1
      M6( 1 , 0 , -1 ) = -1
      M6( -1 , 1 , -1 ) = 3
      M6( 0 , 1 , -1 ) = +1
      M6( 1 , 1 , -1 ) = 3
      M6( -1 , -1 , 0 ) = 3
      M6( 0 , -1 , 0 ) = -1
      M6( 1 , -1 , 0 ) = -1
      M6( -1 , 0 , 0 ) = 3
      M6( 0 , 0 , 0 ) = +1
      M6( 1 , 0 , 0 ) = -1
      M6( -1 , 1 , 0 ) = 3
      M6( 0 , 1 , 0 ) = 3
      M6( 1 , 1 , 0 ) = 3
      M6( -1 , -1 , 1 ) = -1
      M6( 0 , -1 , 1 ) = -1
      M6( 1 , -1 , 1 ) = -1
      M6( -1 , 0 , 1 ) = -1
      M6( 0 , 0 , 1 ) = -1
      M6( 1 , 0 , 1 ) = -1
      M6( -1 , 1 , 1 ) = -1
      M6( 0 , 1 , 1 ) = -1
      M6( 1 , 1 , 1 ) = -1


      !    S, N := 1, -1
      !    E, W := 2, -2
      !    U, D := 3, -3

      !        (  D,     W  ,        N, -,       S,        E,  U )
      angles = (/ pi, pi/2.0, -pi/2.0, 0.0, pi/2.0, -pi/2.0,  0.0  /)
      t = angles(direction)

      ! arrays def. f90 : (\a11,a21,a31,a12,22,32,a13,a32,33\),  columns first:  aCL
      Rx = reshape((/ 1.0,0.0,0.0,0.0,anint(cos(t)),anint(sin(t)),0.0,anint(-sin(t)),anint(cos(t)) /), shape(Rx)) 
      Ry = reshape((/ anint(cos(t)),0.0,anint(-sin(t)),0.0,1.0,0.0,anint(sin(t)),0.0,anint(cos(t)) /), shape(Ry))
      Rz = reshape((/ anint(cos(t)),anint(sin(t)),0.0,anint(-sin(t)),anint(cos(t)),0.0,0.0,0.0,1.0 /), shape(Rz))
      RefSimmX = reshape((/ -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0/), shape(RefSimmX)) ! reflection matrix X
      RefSimmY = reshape((/ 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0/), shape(RefSimmY)) ! reflection matrix Y
      RefSimmZ = reshape((/ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0/), shape(RefSimmZ)) ! reflection matrix Z

      do i = -1, 1
         do j= -1, 1
            do k = -1, 1
               vijk = (/ i, j, k /)
               if (abs(direction).eq.1) then !
                  vijk = MATMUL(Ry,vijk)
                  vijk_ref = MATMUL(RefSimmZ,vijk)
               elseif (abs(direction).eq.2) then
                  vijk = MATMUL(Rx,vijk)
                  vijk_ref = MATMUL(RefSimmZ,vijk)
               elseif (abs(direction).eq.3) then
                  vijk = MATMUL(Ry,vijk)
                  vijk_ref = MATMUL(RefSimmX, vijk)
               end if


               M1_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M1(i,j,k)
               M2_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M2(i,j,k)
               M3_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M3(i,j,k)
               M4_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M4(i,j,k)
               M5_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M5(i,j,k)
               M6_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M6(i,j,k)


               M1_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M1(i,j,k)
               M2_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M2(i,j,k)
               M3_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M3(i,j,k)
               M4_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M4(i,j,k)
               M5_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M5(i,j,k)
               M6_ref(int(vijk_ref(1)),int(vijk_ref(2)),int(vijk_ref(3))) =  M6(i,j,k)
            end do
         end do
      end do



      M1(:,:,:) = M1_rot(:,:,:)
      M2(:,:,:) = M2_rot(:,:,:)
      M3(:,:,:) = M3_rot(:,:,:)
      M4(:,:,:) = M4_rot(:,:,:)
      M5(:,:,:) = M5_rot(:,:,:)
      M6(:,:,:) = M6_rot(:,:,:)



      t = pi/2.0



      if (abs(direction).eq.1)then
         Rdeg = reshape((/ 1.0,0.0,0.0,0.0,anint(cos(t)),anint(sin(t)),0.0,anint(-sin(t)),anint(cos(t)) /), shape(Rdeg))
      elseif(abs(direction).eq.2) then
         Rdeg = reshape((/ anint(cos(t)),0.0,anint(-sin(t)),0.0,1.0,0.0,anint(sin(t)),0.0,anint(cos(t)) /), shape(Rdeg))
      elseif(abs(direction).eq.3) then
         Rdeg = reshape((/ anint(cos(t)),anint(sin(t)),0.0,anint(-sin(t)),anint(cos(t)),0.0,0.0,0.0,1.0 /), shape(Rdeg))
      end if


      n = 0
      do rotation=1, 4

         do ip=1, np

            if (phis(ip)>0.d0) then

               statusf(1:6) = 0
               statusrf(1:6) = 0
               xpstatus(1:2) = 0
               xpstatusrf(1:2) = 0

               i_ip = Lxyz(ip,1)
               j_ip = Lxyz(ip,2)
               k_ip = Lxyz(ip,3)


               do i=-1, 1
                  do j=-1, 1
                     do k=-1, 1
                        x = i_ip+i
                        y = j_ip+j
                        z = k_ip+k

                        if ( x > Lsize(1)-1) then
                           x = -Lsize(1)
                        elseif ( x < -Lsize(1)) then
                           x = Lsize(1)-1
                        end if

                        if (y > Lsize(2)-1) then
                           y = -Lsize(2)
                        elseif ( y < -Lsize(2)) then
                           y = Lsize(2)-1
                        end if

                        if ( z > Lsize(3)-1) then
                           z = -Lsize(3)
                        elseif ( z < -Lsize(3)) then
                           z = Lsize(3) -1
                        end if


                        M(i,j,k) = int(anint(phis(Lxyz_inv(x,y,z))))

                        if(M(i,j,k).ne.M1(i,j,k)) statusf(1)=statusf(1)+1
                        if(M1(i,j,k).eq. 2) then    ! x points - unless one should be +1
                           if(M(i,j,k).eq.+1) xpstatus(1) = xpstatus(1)+1
                           statusf(1)=statusf(1)-1
                        end if


                        if(M(i,j,k).ne.M2(i,j,k)) statusf(2)=statusf(2)+1
                        if(M2(i,j,k).eq. 3) statusf(2)=statusf(2)-1

                        if(M(i,j,k).ne.M3(i,j,k)) statusf(3)=statusf(3)+1
                        if(M3(i,j,k).eq. 3) statusf(3)=statusf(3)-1

                        if(M(i,j,k).ne.M4(i,j,k)) statusf(4)=statusf(4)+1
                        if(M4(i,j,k).eq. 3) statusf(4)=statusf(4)-1

                        if(M(i,j,k).ne.M5(i,j,k)) statusf(5)=statusf(5)+1
                        if(M5(i,j,k).eq. 2) then    ! x points - unless one should be +1
                           if(M(i,j,k).eq.+1) xpstatus(2) = xpstatus(2)+1
                           statusf(5)=statusf(5)-1
                        end if

                        if(M(i,j,k).ne.M6(i,j,k)) statusf(6)=statusf(6)+1
                        if(M6(i,j,k).eq. 3) statusf(6)=statusf(6)-1

                        ! reflections masks

                        if(M(i,j,k).ne.M1_ref(i,j,k)) statusrf(1)=statusrf(1)+1
                        if(M1_ref(i,j,k).eq. 2) then    ! x points - unless one should be +1
                           if(M(i,j,k).eq.+1) xpstatusrf(1) = xpstatusrf(1)+1
                           statusrf(1)=statusrf(1)-1
                        end if

                        if(M(i,j,k).ne.M2_ref(i,j,k)) statusrf(2)=statusrf(2)+1
                        if(M2_ref(i,j,k).eq. 3) statusrf(2)=statusrf(2)-1

                        if(M(i,j,k).ne.M3_ref(i,j,k)) statusrf(3)=statusrf(3)+1
                        if(M3_ref(i,j,k).eq. 3) statusrf(3)=statusrf(3)-1

                        if(M(i,j,k).ne.M4_ref(i,j,k)) statusrf(4)=statusrf(4)+1
                        if(M4_ref(i,j,k).eq. 3) statusrf(4)=statusrf(4)-1

                        if(M(i,j,k).ne.M5_ref(i,j,k)) statusrf(5)=statusrf(5)+1
                        if(M5_ref(i,j,k).eq. 2) then    ! x points - unless one should be +1
                           if(M(i,j,k).eq.+1) xpstatusrf(2) = xpstatusrf(2)+1
                           statusrf(5)=statusrf(5)-1
                        end if

                        if(M(i,j,k).ne.M6_ref(i,j,k)) statusrf(6)=statusrf(6)+1
                        if(M6_ref(i,j,k).eq. 3) statusrf(6)=statusrf(6)-1


                     end do
                  end do
               end do



               temp = statusf(2)*statusf(3)*statusf(4)*statusf(6)*&
                    statusrf(2)*statusrf(3)*statusrf(4)*statusrf(6)
               if(temp.eq.0) then
                  n = n + 1
                  remove(n) = ip
                  finish = finish + 1
               end if


               if(statusf(1).eq.0.and.xpstatus(1)>0)  then

                  n = n + 1
                  remove(n) = ip
                  finish = finish + 1

               else if(statusf(5).eq.0.and.xpstatus(2)>0) then

                  n = n + 1
                  remove(n) = ip
                  finish = finish + 1

               else if(statusrf(1).eq.0.and.xpstatusrf(1)>0) then

                  n = n + 1
                  remove(n) = ip
                  finish = finish + 1

               else if(statusrf(5).eq.0.and.xpstatusrf(2)>0) then

                  n = n + 1
                  remove(n) = ip
                  finish = finish + 1
               end if


            end if


         end do



         do i = -1, 1
            do j= -1, 1
               do k = -1, 1
                  vijk = (/ i, j, k /)
                  vijk = MATMUL(Rdeg,vijk)

                  M1_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M1(i,j,k)
                  M2_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M2(i,j,k)
                  M3_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M3(i,j,k)
                  M4_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M4(i,j,k)
                  M5_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M5(i,j,k)
                  M6_rot(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M6(i,j,k)

                  M1_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M1_ref(i,j,k)
                  M2_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M2_ref(i,j,k)
                  M3_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M3_ref(i,j,k)
                  M4_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M4_ref(i,j,k)
                  M5_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M5_ref(i,j,k)
                  M6_rot_ref(int(vijk(1)),int(vijk(2)),int(vijk(3))) =  M6_ref(i,j,k)
               end do
            end do
         end do

         M1(:,:,:) = M1_rot(:,:,:)
         M2(:,:,:) = M2_rot(:,:,:)
         M3(:,:,:) = M3_rot(:,:,:)
         M4(:,:,:) = M4_rot(:,:,:)
         M5(:,:,:) = M5_rot(:,:,:)
         M6(:,:,:) = M6_rot(:,:,:)

         M1_ref(:,:,:) = M1_rot_ref(:,:,:)
         M2_ref(:,:,:) = M2_rot_ref(:,:,:)
         M3_ref(:,:,:) = M3_rot_ref(:,:,:)
         M4_ref(:,:,:) = M4_rot_ref(:,:,:)
         M5_ref(:,:,:) = M5_rot_ref(:,:,:)
         M6_ref(:,:,:) = M6_rot_ref(:,:,:)

      end do

      if(n>0) then
         do i=1, n-1
            phis(remove(i)) = -1.d0
         end do
      end if

      DEALLOCATE(remove)

    end subroutine thinning


    subroutine thinning_run(phis, lxyz, lxyz_inv, lsize, np)

     implicit none

     integer, intent(in):: np, lsize(1:3)
     real, allocatable, intent(inout) :: phis(:)
     integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
     integer :: finish, skel, attempt


      skel = 1
      attempt = 0

      do while(skel>0)
         finish = 0
         attempt = attempt + 1
         ! the order is important!
         call thinning(phis,3, finish, lxyz, lxyz_inv, lsize, np)
         call thinning(phis,-3,finish, lxyz, lxyz_inv, lsize, np)
         call thinning(phis,-1,finish, lxyz, lxyz_inv, lsize, np)
         call thinning(phis,1, finish, lxyz, lxyz_inv, lsize, np)
         call thinning(phis,2, finish, lxyz, lxyz_inv, lsize, np)
         call thinning(phis,-2,finish, lxyz, lxyz_inv, lsize, np)

         if(finish.eq.0.or.attempt.eq.1000) skel = -1
       end do

    end subroutine thinning_run

end module thinning_m

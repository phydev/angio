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

module misc_m  

  use global_m
  use dderivatives_m  
  
  implicit none
  
  private
  
  public :: heaviside, gen_grid_cell_domain, spherical_surface, gen_cell_points, cahn_hilliard, ran2

  contains

    function heaviside(x)
      
      implicit none
      
      real :: x, heaviside
      
      if ( x<0.d0) then
         heaviside = 0.d0
      else
         heaviside = 1.d0
      end if
    end function heaviside

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
      ! string representation  for instance, given the CHARACTER value '154',
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
     

  end module misc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

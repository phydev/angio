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

module blood_flow_m
  

  implicit none
   
  private
   
  public :: flow_calc
   
  type hydro_t
     integer :: ip
     integer :: m, n ! nodes 
     real :: flow
  end type hydro_t
   
contains
   
  subroutine flow_calc(phis, flow, lxyz, lxyz_inv, Lsize, np)    

    implicit none
    
    real, allocatable, intent(in) :: phis(:)
    real, allocatable, intent(inout) :: flow(:)
    integer, allocatable, intent(in) :: Lxyz(:,:), Lxyz_inv(:,:,:) 
    integer, intent(in) :: Lsize(3), np	 

    ! internal variables
    integer :: ip, ip2
    integer :: i, j, k, x, y, z
    real, allocatable :: neighbours(:), permittivity(:,:), B(:), midp(:), phi_copy(:)
    integer, allocatable :: nodes_matrix(:,:), nodes_matrix_count(:,:), IPIV(:)
    integer :: node_inout(1:2), paths(1:14), paths_copy(1:14), nodes(1:5000), path_length, found, attempt
    integer :: Li(2), MN_nodes(2), path_length_v(1000), connect, path_ip(1000),  connections(1000) 
    integer :: np_skel, m, nodes_connecteds, ierr, Lmn(1:1000,1:1000),  n, ip_old, np_path, np_nodes
    type(hydro_t), allocatable:: hydro(:)
    real, allocatable :: output(:)
    ! vertex variables   
    real, allocatable :: length_vertex_connections(:,:)
    integer, allocatable ::  vertex(:,:), vertex_full(:,:), vertex_vector(:), multiplier_vector(:), cluster(:,:)
    integer, allocatable :: nodes_vertex_id(:)
    integer :: n_vertex
    real :: distance(3), distance_abs
    ! new
    integer, allocatable :: path_nodes(:), temp_node_index(:)

    ALLOCATE(path_nodes(10000))
    ALLOCATE(temp_node_index(10000))
    ALLOCATE(neighbours(1:np))
    ALLOCATE(hydro(1:np)) 

    flow(:) = -1.d0
    hydro(:)%ip = 0 
    hydro(:)%flow = 0
    hydro(:)%m = 0 
    hydro(:)%n = 0 
    np_skel = 0

    np_nodes = 0
    do ip=1, np 
       neighbours(ip) = 0
       

       if(phis(ip)>0.d0) then
          np_skel = np_skel + 1
          hydro(np_skel)%ip = ip ! saving all skeleton points separated	  


          do i=-1,1
             do j=-1,1
                do k=-1,1

                   x = lxyz(ip,1) + i 
                   y = lxyz(ip,2) + j 
                   z = lxyz(ip,3) + k 

                   ip2 = lxyz_inv(x,y,z)

                   if ( phis(ip2).gt.0.and.ip2.ne.ip ) then
                      neighbours(ip) = neighbours(ip) + 1.0
                   end if



                end do
             end do
          end do

          if (neighbours(ip).gt.2 .or. neighbours(ip).eq.1) then ! real nodes or tips without connections
             np_nodes = np_nodes + 1
             nodes(np_nodes) = ip 

          end if
       end if
    end do ! nodes 


    ! search by principal nodes in/out 
    Li(1) = -Lsize(3)
    Li(2) =  Lsize(3) - 1 

    do k=1,2
       !obsolete: temp = 10000
       found = 0
       do i = -5,5!-Lsize(2), Lsize(2)-1
          do j = -5,5 !-Lsize(3), Lsize(3)-1 
             ip = lxyz_inv(i,j,Li(k))
             ! obsolete condition: .and. (j**2 + k**2.lt.temp)
             if( phis(ip).gt.0  ) then 
                do m=1,np_nodes 
                   if(ip.eq.nodes(m)) then 
                      node_inout(k) = m 
                      found = 1 
                   end if
                end do

                if (found.eq.0) then 
                   nodes(np_nodes+1) = ip 
                end if

             end if

          end do
       end do

       if(found.eq.0) then 
          np_nodes = np_nodes + 1 
          node_inout(k) = np_nodes
       end if

    end do


    ALLOCATE(nodes_matrix(1:np_nodes,1:np_nodes))  ! dynamical memory alloc
    ALLOCATE(nodes_matrix_count(1:np_nodes,1:np_nodes))
    nodes_matrix(:,:) = 0 

    ! searching by nodes neighbours
    paths(:) = 0


    
    do ip=1, np_nodes

       np_path = 0 

       ! searching by all neighbours of the node
       ! to walk by them posteriorly 

       do i=-1,1
          do j=-1,1
             do k=-1,1

                found = 0

                x = lxyz(nodes(ip),1) + i 
                y = lxyz(nodes(ip),2) + j 
                z = lxyz(nodes(ip),3) + k 


                ip2 = lxyz_inv(x,y,z)

                if( phis(ip2) > 0.d0 .and. ip2.ne.nodes(ip))  then

                   do m=1, np_nodes 
                      if(nodes(m).eq.ip2) then
                         found = m 
                      end if
                   end do

                   if(found.eq.0) then 
                      np_path = np_path + 1
                      paths(np_path) =  ip2
                      path_nodes(np_path) = 0
                   else 
                      np_path = np_path + 1
                      paths(np_path) = ip2
                      path_nodes(np_path) = 1
                      nodes_matrix(ip,found) = 1
                      temp_node_index(np_path) =  found
                   end if

                end if


             end do
          end do
       end do


       ! walking for all paths 

       connections(:) = 0
       nodes_connecteds = 0
       nodes_matrix_count(:,:) = 0

       do n=1, np_path
          ip_old = nodes(ip) 
          found = 0 
          path_length = 0 
          MN_nodes(1) = ip
          MN_nodes(2) = MN_nodes(1)	
          
          if(path_nodes(n).eq.1) then

             ip2 = paths(n)
             m = temp_node_index(n)
             path_length = 1
             nodes_connecteds = nodes_connecteds + 1

             path_length_v(n) = 1
             path_ip(path_length) = ip2
             Lmn(1:path_length, n) = path_ip(1:path_length)
             Lmn(path_length, n) = Lmn(path_length, n) + paths(n)
             
             nodes_matrix_count(ip,m) = nodes_matrix_count(ip,m) + 1
             
             if ( nodes_matrix_count(ip,m) > 1) then
                nodes_matrix(ip,m) = nodes_matrix(ip,m) + path_length 
                nodes_matrix(m,ip) = nodes_matrix(m,ip) + path_length 
             else
                nodes_matrix(ip,m) = path_length
                nodes_matrix(m,ip) = path_length
             end if
             
             path_length = path_length
             MN_nodes(2) = m
             found = 1  
          end if

          do while(found.eq.0)  
             attempt = 0 

             do i=-1,1
                do j=-1,1
                   do k=-1,1

                      x = lxyz(paths(n),1) + i  
                      y = lxyz(paths(n),2) + j 
                      z = lxyz(paths(n),3) + k 

                      ip2 = lxyz_inv(x,y,z)

                      if( (phis(ip2).gt.0).and.(ip2.ne.ip_old)&
                           .and.(ip2.ne.paths(n))) then  ! ip_old is to guarantee no return


                         do m=1, np_nodes


                            if (nodes(m).eq.ip2) then	

                               nodes_connecteds = nodes_connecteds + 1

                               path_length_v(n) = path_length+1
                               path_ip(path_length + 1) = ip2
                               Lmn(1:path_length, n) = path_ip(1:path_length)
                               Lmn(path_length+1, n) = Lmn(path_length+1, n) + paths(n)

                               nodes_matrix_count(ip,m) = nodes_matrix_count(ip,m) + 1

                               if ( nodes_matrix_count(ip,m) > 1) then
                                  nodes_matrix(ip,m) = nodes_matrix(ip,m) + path_length + 1
                                  nodes_matrix(m,ip) = nodes_matrix(m,ip) + path_length + 1
                               else
                                  nodes_matrix(ip,m) = path_length + 1
                                  nodes_matrix(m,ip) = path_length + 1
                               end if

                               path_length = path_length+1
                               MN_nodes(2) = m
                               found = 1  

                            end if
                         end do

                         if(found.eq.0) then 

                            path_ip(path_length+1) = ip2 !paths(n)  ! saving the path for later identify the flow 
                            Lmn(path_length+1,np_path) = paths(n)
                            ip_old = paths(n) 
                            path_length = path_length + 1
                            path_length_v(n) = path_length 
                            paths(n) = ip2 
                            attempt = 0 
                            hydro(ip2)%m = MN_nodes(1)
                         end if

                      else

                         attempt = attempt + 1
                      end if

                   end do
                end do
             end do


          end do ! do while found 

       
          do i=1, path_length_v(n)

             ! ip, m
             hydro(path_ip(i))%ip = path_ip(i)
             hydro(path_ip(i))%m = MN_nodes(1)
             hydro(path_ip(i))%n = MN_nodes(2)

          end do

       end do ! np_path 

    end do ! np_nodes
    
    ALLOCATE(permittivity(np_nodes,np_nodes))  ! dynamical memory alloc
    ALLOCATE(B(1:np_nodes))
    ALLOCATE(IPIV(1:np_nodes))


    permittivity(:,:) = 0.d0

    do m=1, np_nodes
       do n=1, np_nodes
          if(nodes_matrix(m,n)>0.d0 .and.m.ne.n) then 
             permittivity(m,n) = 1.d0/real(nodes_matrix(m,n))
          end if

       end do
    end do
 
    permittivity(node_inout(1),1:np_nodes) = 0.d0
    permittivity(node_inout(2),1:np_nodes) = 0.d0
    !(row,col)

    do m=1, np_nodes		
       permittivity(m,m) = permittivity(m,m) - sum(permittivity(m,1:np_nodes))
    end do

    permittivity(node_inout(1),node_inout(1)) = 1.0 
    permittivity(node_inout(2),node_inout(2)) = 1.0 		


    B(1:np_nodes) = 0.d0
    B(node_inout(1)) = 1.d0
    IPIV(1:np_nodes) = 0

    call DGESV(np_nodes, 1, permittivity, np_nodes, IPIV, B, np_nodes, ierr)

    do ip=1, np			
       if(hydro(ip)%m > 0) then             
         
          flow(ip) = abs((B(hydro(ip)%m) - B(hydro(ip)%n))/&
              (real(nodes_matrix(hydro(ip)%m, hydro(ip)%n ))) )

       end if
    end do



    DEALLOCATE(neighbours)
    DEALLOCATE(hydro)
    DEALLOCATE(nodes_matrix)  ! dynamical memory dealloc ok
    DEALLOCATE(permittivity)
    DEALLOCATE(B) 
    DEALLOCATE(IPIV)


  end subroutine flow_calc
          
end module blood_flow_m

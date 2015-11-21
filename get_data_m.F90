module get_data_m

  implicit none

  private

  public :: run_get_data,  heaviside


  type hydro_t
     integer :: ip
     integer :: m, n ! nodes 
     real :: flow
  end type hydro_t

contains

  subroutine run_get_data(xv,yv,zv,f,lx,ly,lz,nbranches)

    implicit none
    ! external
    real,  intent(in) :: f(:) 
    integer, intent(in) :: xV(:), yv(:), zv(:)
    integer, intent(in) :: lx, ly, lz
    ! GET_NBRANCHES VARIABLES
    integer, intent(out) :: nbranches

    ! intern
    real, allocatable :: phis(:)
    integer, allocatable :: lxyz(:,:), lxyz_inv(:,:,:)
    integer :: Lsize(3), boundary_points, np, ndim, np_part
    logical :: periodic, boundary
    integer :: i, j, k, l, m, n, ip, ip_part
    real :: hs(1:3)
    ! GET_NODE_MATRIX VARIABLES
    integer ::  ip2
    integer ::  x, y, z
    real, allocatable :: neighbours(:), midp(:), phi_copy(:)
    integer, allocatable :: nodes_matrix(:,:), nodes_matrix_count(:,:), IPIV(:)
    integer :: node_inout(1:2), paths(1:14), paths_copy(1:14), nodes(1:5000), path_length, found, attempt
    integer :: Li(2), MN_nodes(2), path_length_v(1000), connect, path_ip(1000),  connections(1000) 
    integer :: np_skel, nodes_connecteds, ierr, Lmn(1:1000,1:1000),  ip_old, np_path, np_nodes
    type(hydro_t), allocatable:: hydro(:)
    ! new
    integer, allocatable :: path_nodes(:), temp_node_index(:)



    Lsize(1) = lx
    Lsize(2) = ly
    Lsize(3) = lz
    

    boundary_points = 2
    np = 8 * Lsize(1) * Lsize(2) * Lsize(3)
    np_part = 8*(Lsize(1)+2*boundary_points)*(Lsize(2)+2*boundary_points)*(Lsize(3)+2*boundary_points)
    periodic = .false.

    ALLOCATE(lxyz(np_part,1:3))
    ALLOCATE(lxyz_inv(-Lsize(1)-boundary_points:Lsize(1)+boundary_points, &
         -Lsize(2)-boundary_points:Lsize(2)+boundary_points, &
         -Lsize(3)-boundary_points:Lsize(3)+boundary_points))
    ALLOCATE(phis(np))

    !######################################################################
    !################### SPACE_INIT SUBROUTINE ############################
    !######################################################################
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
    
    !#################################################################################################
    !#################################################################################################
    !#################################################################################################


    phis(:) = 0.d8
    do ip=0, size(f)
       phis(lxyz_inv(xv(ip),yv(ip),zv(ip))) = 1.0
    end do
          



          ! do ip=1, np
          !    if(phis(ip)>0.d8) then
          !       write(*,*) lxyz(ip,1:3), phis(ip)
          !    end if
          ! end do
    !#################################################################################################
    !#################################################################################################
    !#############  get_node_matrix(phis,lxyz,lxyz_inv, Lsize, np)     ###############################
    !#################################################################################################
    !#################################################################################################
    ALLOCATE(path_nodes(10000))
    ALLOCATE(temp_node_index(10000))


    ALLOCATE(neighbours(1:np))
    ALLOCATE(hydro(1:np)) 
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
      
       found = 0
       do i = -5,5
          do j = -5,5  
             ip = lxyz_inv(i,j,Li(k))
             
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

                            path_ip(path_length+1) = ip2  ! saving the path for later identify the flow 
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

             hydro(path_ip(i))%ip = path_ip(i)
             hydro(path_ip(i))%m = MN_nodes(1)
             hydro(path_ip(i))%n = MN_nodes(2)

          end do

       end do ! np_path 
       
    end do


    nbranches = 0
    do n=1, np_nodes-1
       do m=n+1, np_nodes
          if(nodes_matrix(m,n)>1) then
             nbranches = nbranches + 1
          end if
       end do

    end do

  
    DEALLOCATE(path_nodes)
    DEALLOCATE(temp_node_index)
    DEALLOCATE(neighbours)
    DEALLOCATE(hydro)  
    DEALLOCATE(nodes_matrix)
    DEALLOCATE(nodes_matrix_count)
    DEALLOCATE(lxyz)
    DEALLOCATE(lxyz_inv)
    DEALLOCATE(phis)
    

  end subroutine run_get_data



  
  function heaviside(x)
    
    implicit none
    
    real :: x, heaviside
    
    if ( x<0.d0) then
       heaviside = 0.d0
    else
       heaviside = 1.d0
    end if
  end function heaviside
  


end module get_data_m

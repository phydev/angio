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
   
  subroutine flow_calc(phis,lxyz,lxyz_inv, Lsize, np)    

    implicit none
    
    real, allocatable, intent(in) :: phis(:)
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
    integer, allocatable ::  vertex(:,:), vertex_full(:,:), vertex_vector(:), multiplier_vector(:), cluster(:,:)
    integer, allocatable :: nodes_vertex_id(:)
    integer :: n_vertex
    real :: distance(3), distance_abs

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
             write(*,*) lxyz(nodes(np_nodes),1:3), np_nodes
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
          !write(*,*) lxyz(nodes(np_nodes),1:3), np_nodes
       end if

    end do


    ALLOCATE(nodes_matrix(1:np_nodes,1:np_nodes))  ! dynamical memory alloc
    ALLOCATE(nodes_matrix_count(1:np_nodes,1:np_nodes))
    nodes_matrix(:,:) = 0 

    ! searching by nodes neighbours
    paths(:) = 0
    write(*,'(A,I10)') "np_nodes", np_nodes

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
                   else 
                      nodes_matrix(ip,found) = 1
                   end if

                end if


             end do
          end do
       end do


       ! walking for all paths 

       connections(:) = 0
       nodes_connecteds = 0
       nodes_matrix_count(:,:) = 0
       write(*,*) "node",ip, "np_path", np_path

       do n=1, np_path
          ip_old = nodes(ip) 
          found = 0 
          path_length = 0 
          MN_nodes(1) = ip
          MN_nodes(2) = MN_nodes(1)	
          write(*,*)   "path_n=",n

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
             !write(*,'(I10,I10,I10)') hydro(path_ip(i))%ip , hydro(path_ip(i))%m, hydro(path_ip(i))%n

          end do

       end do ! np_path 

    end do ! np_nodes



    ALLOCATE(vertex(1:np_nodes,1:np_nodes) )
    ALLOCATE(vertex_full(np_nodes,np_nodes))
    ALLOCATE(vertex_vector(np_nodes))
    ALLOCATE(multiplier_vector(np_nodes))

    vertex(:,:) = 0

    do i=1, np_nodes-1
       vertex(i,i) = 1
       do j=i+1, np_nodes

          distance(1:3) = lxyz(nodes(i),1:3) - lxyz(nodes(j),1:3)
          distance = distance*distance
          distance_abs = sqrt(distance(1) + distance(2) + distance(3))
          
          if(distance_abs<=3.d0) then
           
             vertex(i,j) = 1
          end if

       end do
    end do

    vertex_full(:,:) = vertex(:,:)

    do i=1, np_nodes
       do j=1, np_nodes
          if(vertex(i,j)>0) then
             vertex_full(i,:) = vertex_full(i,:) + vertex(j,:)
          end if
       end do
    end do


    do i=1, np_nodes
       do j=1, np_nodes
          if(vertex_full(i,j)>0) vertex_full(i,j) = 1
       end do
    end do

    vertex_vector(:) = 1
    multiplier_vector(:) = 0

    do i=1, np_nodes - 1
       do j = i + 1, np_nodes                         
          multiplier_vector(:) = vertex_full(i,:)*vertex_full(j,:)
          if( sum(multiplier_vector(1:5)) .ne. 0) vertex_vector(j) - 1	
       end do
    end do

    ALLOCATE(cluster(np_nodes,np_nodes))
    n_vertex = 0
    do i=1, np_nodes
       if( multiplier_vector(i) .eq. 1) then
          n_vertex = n_vertex + 1         
          cluster(n_vertex,1:np_nodes) = vertex_full(i,1:np_nodes)
       end if 
    end do

    do i=1, np_nodes
       do j=1, n_vertex
          
          if(cluster(j,i).eq.1) then
             nodes_vertex_id(i) = j
          end if

       end do
    end do
    

    
    ALLOCATE(permittivity(n_vertex,n_vertex))  ! dynamical memory alloc
    ALLOCATE(B(1:np_nodes))
    ALLOCATE(IPIV(1:np_nodes))

    permittivity(:,:) = 0.d0



    do m=1, np_nodes
       do n=1, np_nodes

          if(nodes_matrix(m,n)>0.d0 .and. nodes_vertex_id(m).ne.nodes_vertex_id(n)) then 
             permittivity(nodes_vertex_id(m),nodes_vertex_id(n)) = 1.d0/real(nodes_matrix(m,n))
          end if

       end do
    end do

    permittivity(nodes_vertex_id(node_inout(1)),1:n_vertex) = 0.d0
    permittivity(nodes_vertex_id(node_inout(2)),1:n_vertex) = 0.d0
    !(row,col)

    do m=1, n_vertex 		
       permittivity(m,m) = permittivity(m,m) - sum(permittivity(m,1:n_vertex))
    end do

    permittivity(nodes_vertex_id(node_inout(1)),nodes_vertex_id(node_inout(1))) = 1.0 
    permittivity(nodes_vertex_id(node_inout(2)),nodes_vertex_id(node_inout(2))) = 1.0 		


    B(1:n_vertex) = 0.d0
    B(nodes_vertex_id(node_inout(1))) = 1.d0
    IPIV(1:n_vertex) = 0

    write(*,*) "lel"
    call DGESV(n_vertex, 1, permittivity, n_vertex, IPIV, B, n_vertex, ierr)
    write(*,'(F10.2,F10.2,F10.2,F10.2)') B(1:n_vertex)
    write(*,'(I10)') ierr
    write(*,*) "np_skel", np_skel


    do ip=1, np			
       if(hydro(ip)%m > 0) then             
         
          hydro(ip)%flow = (B(nodes_vertex_id(hydro(ip)%m) ) - B(nodes_vertex_id(hydro(ip)%n) ))/&
               (real(nodes_matrix( nodes_vertex_id(hydro(ip)%m), nodes_vertex_id(hydro(ip)%n) ))) 
          write(*,'(I10,I10,I10,F10.5)') lxyz(ip,1:3), abs(hydro(ip)%flow)
       end if


    end do



    !write(filename,'(I5)') nstep
    !OPEN (UNIT=nstep,FILE=dir_name//'/xyzs'//trim(filename)//".out")
    !OPEN (UNIT=109929,FILE='xyzss.out')
    !finish = 0 
    !write(nstep,'(A)') "#  x       y       z       phi       flow"

    !do ip=1, np_skel
    !   i = lxyz(hydro(ip)%ip,1)
    !   j = lxyz(hydro(ip)%ip,2)
    !   k = lxyz(hydro(ip)%ip,3)
    !write(109929,'(I5,I5,I5,F10.2,F10.2)')  i,j,k, phis(hydro(ip)%ip), hydro(ip)%flow       
    !   write(*,*) i, j, k, hydro(ip)%flow
    !write(*,'(I5,I5,I5,F10.2,F10.2)')  i,j,k, phis(hydro(ip)%ip), hydro(ip)%flow
    !end do
    !write(*,*)"N_skel", np_skel 
    !close(109929)
    !CLOSE(nstep)

    DEALLOCATE(neighbours)
    DEALLOCATE(hydro)
    DEALLOCATE(nodes_matrix)  ! dynamical memory dealloc ok
    DEALLOCATE(permittivity)
    DEALLOCATE(B) 
    DEALLOCATE(IPIV)



  end subroutine flow_calc
          
end module blood_flow_m

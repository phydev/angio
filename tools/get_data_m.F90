module get_data_m

  implicit none

  private

  public :: run_get_data, space_init, nodo

  type nodo
     integer :: ip, vertex    ! ip global
  end type nodo

contains

  subroutine run_get_data(phis,Lsize,lxyz,lxyz_inv,np,periodic,output,genus, branch_length, dead_volume)

    implicit none
    ! external
    real, allocatable,  intent(in) :: phis(:)
    integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
    integer, intent(in) :: Lsize(3), np
    integer, intent(out) :: output, genus, branch_length
    real, intent(out) :: dead_volume
    logical, intent(in) :: periodic
    ! GET_NBRANCHES VARIABLES
    integer :: nbranches, carry_on
    real :: diameter
    integer :: i, j, k, l, m, n, ip, ip_part
    real :: hs(1:3), drtemp, a,b,c
    ! GET_NODE_MATRIX VARIABLES
    integer ::  ip2
    integer ::  x, y, z
    real, allocatable :: neighbours(:)
    integer, allocatable :: nodes_matrix(:,:), nodes_matrix_count(:,:), verify_nodes(:)
    integer :: node_inout(1:2), paths(1:14), nodes(1:5000), path_length, found, attempt
    integer :: Li(2), MN_nodes(2), path_length_v(1000), path_ip(1000),  connections(1000)
    integer :: np_skel, nodes_connecteds, Lmn(1:1000,1:1000),  ip_old, np_path, np_nodes
    ! new
    integer, allocatable :: path_nodes(:), temp_node_index(:)
    real :: M_Pi, dlambda
    integer :: M_ONE
    integer :: adj6(1:6), adj18(1:18), real_neighbour, real_branches, ia, real_nodes, vertex_count
    logical :: neighbour_tf
    real :: adjacent6, adjacent18
    type(nodo), allocatable :: nnodes(:)

    M_Pi = 3.14159265359
    M_ONE = 1


    np_skel = 0
    np_nodes = 0
    real_nodes = 0
   ALLOCATE(path_nodes(10000))
   ALLOCATE(temp_node_index(10000))
   ALLOCATE(neighbours(1:np))

    do ip=1, np
       neighbours(ip) = 0
       real_neighbour = 0

       if(phis(ip)>0.d0) then
          np_skel = np_skel + 1

          do i=-1,1
             do j=-1,1
                do k=-1,1

                   x = lxyz(ip,1) + i
                   y = lxyz(ip,2) + j
                   z = lxyz(ip,3) + k

                   ip2 = lxyz_inv(x,y,z)

                   dlambda = sqrt(real(i*i) + real(j*j) + real(k*k))


                   if ( phis(ip2).gt.0.and.ip2.ne.ip ) then
                     neighbours(ip) = neighbours(ip) + 1.0

                   end if



                end do
             end do
          end do

          !if(real_neighbour.gt.2 .or. real_neighbour.eq.1) then
            !print*, lxyz(ip,1:3)
          !  real_nodes = real_nodes + 1
          !  write(*,'(I10,I10,I10,I10)') lxyz(ip,1:3), real_nodes
          !  np_nodes = np_nodes + 1
          !  nodes(np_nodes) = ip
          !end if

          if (neighbours(ip).gt.2 .or. neighbours(ip).eq.1) then ! real nodes or tips without connections
             np_nodes = np_nodes + 1
             nodes(np_nodes) = ip
            ! write(*,'(I10,I10,I10,I10)') lxyz(ip,1:3), np_nodes
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
   ! write(*,'(A,I10)') "np_nodes", np_nodes


    do ip=1, np_nodes

       np_path = 0
       path_nodes(:) = 0
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

                if( phis(ip2) > 0.d0 .and. ip2.ne.nodes(ip))  then ! .and. ip2.ne.nodes(ip)

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
       !write(*,*) "node",ip, "np_path", np_path

       do n=1, np_path
          ip_old = nodes(ip)
          found = 0
          path_length = 0
          MN_nodes(1) = ip
          MN_nodes(2) = MN_nodes(1)

          ! nodes that are first neighbours
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

                               !path_length = path_length+1
                               MN_nodes(2) = m
                               found = 1
                               EXIT
                            end if
                            if(found.gt.0) exit
                         end do

                         if(found.eq.0) then

                            path_ip(path_length+1) = paths(n)  ! saving the path for later identify the flow
                            Lmn(path_length+1,np_path) = paths(n)
                            ip_old = paths(n)
                            path_length = path_length + 1
                            path_length_v(n) = path_length
                            paths(n) = ip2
                            attempt = 0

                         end if

                      else

                         attempt = attempt + 1
                      end if
                      if(found.gt.0) EXIT
                   end do
                   if(found.gt.0) EXIT
                end do
                if(found.gt.0) EXIT
             end do

             if(attempt>26) then
                found = 1
                EXIT
             end if

          end do ! do while found
       end do ! np_path
    end do ! np_nodes


    !
    ! ALLOCATE(verify_nodes(1:np_nodes))
    ! verify_nodes(:)  = 0
    ! do n=1, np_nodes
    !   do m=1, np_nodes
    !     if(nodes_matrix(n,m)>8) then
    !       !real_nodes = real_nodes + 1
    !       verify_nodes(n) = 1
    !       verify_nodes(m) = 1
    !     end if
    !   end do
    ! end do
    ! do m=1, np_nodes
    !   if(verify_nodes(m).eq.1) real_nodes = real_nodes + 1
    ! end do



    allocate(nnodes(1:np_nodes))
    nnodes(1:np_nodes)%vertex = 0
    nnodes(1:np_nodes)%ip = nodes(1:np_nodes)

    vertex_count = 0
    do n=1, np_nodes-1
       do m=n+1, np_nodes
      !   print*, nodes_matrix(m,n)
         if(nodes_matrix(m,n).eq.1 ) then

           if(nnodes(n)%vertex.eq.0.and.nnodes(m)%vertex.eq.0) then
            vertex_count = vertex_count + 1
            nnodes(n)%vertex = vertex_count
            nnodes(m)%vertex = vertex_count

          elseif(nnodes(n)%vertex.eq.0.and.nnodes(m)%vertex.ne.0) then
            nnodes(n)%vertex = nnodes(m)%vertex
          elseif(nnodes(n)%vertex.ne.0.and.nnodes(m)%vertex.eq.0) then
            nnodes(m)%vertex = nnodes(n)%vertex
          end if

         end if
       end do
    end do

    real_nodes = 0
    do n=1, np_nodes
      if(nnodes(n)%vertex.eq.0) then
        real_nodes = real_nodes + 1
      end if
    !  print*, lxyz(nnodes(n)%ip,1:3), nnodes(n)%vertex
    end do


    nbranches = 0
    real_branches = 0
    branch_length = 0
    dead_volume = 0.d0
    !print*, "distancias"
    do n=1, np_nodes-1
       do m=n+1, np_nodes
          if(nodes_matrix(m,n)>8 ) then
            !print*, n, m, nodes_matrix(m,n)
             branch_length = branch_length + nodes_matrix(m,n)

             real_branches = real_branches + 1
          elseif(nodes_matrix(m,n)<8.and.nodes_matrix(m,n)>4) then
             dead_volume = dead_volume + (4./3.)*3.14*4.0**3
          end if

          if(nodes_matrix(m,n)>1 ) then
             nbranches = nbranches + 1
          end if
       end do
    end do
!print*, 'nodes, branches, genus'
    !print*, real_nodes, vertex_count, nbranches, nbranches- (real_nodes+vertex_count) +1

    output = real_branches!nbranches
    genus = nbranches - (real_nodes+vertex_count) + 1!real_branches-np_nodes+1
  !  print*, real_branches, genus
    !print*,"genus      nbranches        rb            n        real_nodes"
    !print*, genus, nbranches, real_branches, np_nodes, real_nodes
	!write(*,*) nbranches
  !print*, "      real_nodes     np_nodes     genus         nbranches       real_branches"
  !print*, real_nodes, np_nodes, genus, nbranches, real_branches

    DEALLOCATE(nnodes)
    DEALLOCATE(path_nodes)
    DEALLOCATE(temp_node_index)
    DEALLOCATE(neighbours)
    DEALLOCATE(nodes_matrix)
    DEALLOCATE(nodes_matrix_count)


  end subroutine run_get_data

  subroutine space_init(Lsize, lxyz, lxyz_inv, boundary_points, np, ndim, periodic)

      implicit none

      ! input/output variables
      integer, intent(in) ::  Lsize(3), boundary_points
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
                    n = SIGN(1,k)*Lsize(3) - hs(3)*SIGN(1,k)
                 end if
              end if
               !
              !  if( abs(i)>Lsize(1)-hs(1)) then
              !     boundary = .true.
              !     l = i - SIGN(1,i)*(2*Lsize(1)-hs(1))
              !  end if
               !
              !  if( abs(j)>Lsize(2)-hs(2)) then
              !     boundary = .true.
              !     m = j  - SIGN(1,j)*(2*Lsize(2)-hs(2))
              !  end if
               !
              !     if(periodic) then
              !        n = k - SIGN(1,k)*(2*Lsize(3)-hs(3))
              !     else
              !        n = k - SIGN(1,k)
              !     end if


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

               ! the updating of boundaries should gra
               ! cell(ip_part) = cell( lxyz_inv( lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ) )
               ! because:
               ! lxyz(np_part,1), lxyz(np_part,2), lxyz(np_part,3) ->(l, m, n) in the bulk
               ! then lxyz_inv(m,n,l) will give the updated value in the bulk
               ! for the respective boundary position cell(ip_part)

            end do
         end do
      end do

    end subroutine space_init



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

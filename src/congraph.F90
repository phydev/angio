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

module congraph_m

   implicit none

   private

   public :: verify_graph_connection

   contains

    subroutine verify_graph_connection(phis, graph, initial_node_ip, lxyz, lxyz_inv)

      implicit none

      real, allocatable, intent(in) :: phis(:)
      real, allocatable, intent(inout) :: graph(:)
      integer, intent(in) :: initial_node_ip
      integer, allocatable, intent(in) :: lxyz(:,:), lxyz_inv(:,:,:)
      ! internal
      integer :: ip, ip_o, i, j, k
      integer :: npaths, ipaths, paths(1000), paths_temp(1000)
      logical :: loop

     ! tenho que partir do nodo principal e ir marcando com 1 cada ponto conectado com ele
     ! a cada bifurcacao devo guardar os caminhos ainda por percorrer
     ! os caminhos que ainda restarem deverão ser percorridos apenas no fim do primeiro caminho estar completo.
     ! para evitar que o loop retorne por um caminho já percorrido, deve-se realizar uma verificacao pela variavel
     ! que identifica se o ponto está ligado ao nodo principal.


      ! npaths : number of alternatives paths due to bifurcations

      graph(:) = -1.d0
      ip_o = initial_node_ip
      npaths = 0
      graph(ip_o) = 1.0
      loop = .true.

      do while(loop)

         ipaths = 0

         do i=-1,1
            do j=-1,1
               do k=-1,1
       
                  ip = lxyz_inv(lxyz(ip_o,1) +  i, lxyz(ip_o,2) +  j, lxyz(ip_o,3) +  k)

                  if(phis(ip).gt.0 .and. graph(ip).lt.0) then
                     ipaths = ipaths + 1
                     paths_temp(ipaths) = ip
                  end if

               end do
            end do
         end do

         if(ipaths.gt.0) then

            graph(paths_temp(1)) = 1
            ip_o = paths_temp(1)

            if(ipaths.gt.1) then
               do i=2, ipaths
                  npaths = npaths + 1
                  paths(npaths) = paths_temp(i)
               end do
            end if

         else ! em caso de ter chegado ao fim da linha.
            if(npaths.gt.0) then
               ip_o = paths(npaths) ! percorre os caminhos que faltam do fim da lista para o começo
               npaths = npaths - 1  !
            else
               loop = .false.
               EXIT
            end if
         end if

      end do


    end subroutine verify_graph_connection



end module congraph_m

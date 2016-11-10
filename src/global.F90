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

module global_m

  implicit none

  private

  public :: tip_cell_t, mesh_t, sec_mesh_t

  type tip_cell_t
     integer :: ip    ! ip global
     real :: x, y, z ! global mesh position
     real :: phi     ! phi_c inside tipcell
  end type tip_cell_t

  type mesh_t
     real :: phi, T, mu, lapl_phi, lapl_T, lapl_mu
     real :: alpha_p, consumption
     integer :: source
  end type mesh_t

  type sec_mesh_t
     real :: phi
  end type sec_mesh_t



      ! begin parameters
      integer, public :: i, j, k
      real, public :: cell_radius, diffusion_const, interface_width, vegf_p, vegf_c, diff_oxy_length, vegf_rate, &
           vegf_source_conc, prolif_rate, vessel_radius, dt, chi, vegf_grad_min, vegf_grad_max
      integer, public :: tstep, iseed, dr(3)
      character(len=3), public :: dir_name
      character(len=6), public :: file_name, file_id
      character(len=17), public :: file_phi, file_t
      character(len=18),public :: file_phis
      character(len=3),public :: sim_id
      real, public :: notch_distance, depletion_weight
      ! mesh variables
      integer, public, allocatable :: lxyz(:,:), lxyz_inv(:,:,:),  grid_cell_domain(:), d2sphere(:), sphere(:,:), vegf_all(:,:)
      integer, public :: Lsize(1:3), boundary_points
      integer, public :: np, np_part, ip, n_max_tipc
      real, public, allocatable :: gg(:,:), lapl(:), f(:), flow(:), flow_full(:)
      type(mesh_t), public, allocatable :: cell(:)
      type(sec_mesh_t), public, allocatable :: necrotic_tissue(:)
      type(tip_cell_t), public, allocatable :: tipc(:)
      ! vegf related variables
      integer, public :: n_source = 0, source_max, n_source_initial, np_vegf_all
      integer, public, allocatable :: vegf_xyz(:,:)
      ! misc
      integer, public :: nstep, counter = 0, ndim, output_period, extra_steps
      real, public :: hs, time_init, time_end, ctime, phi_max
      logical, public :: periodic
      ! thinning
      logical, public :: thinning
      logical, public :: calculate_flow
      integer, public :: flow_count, calc_flow_period, np_sphere
      real, public, allocatable :: phis(:), graph(:)
      ! np_vegf_s  - number of points on the spherical surface representing the limit signaling of vegf source
      ! np_tip_s   - number of points on the spherical surface representing tip cell limits
      ! np_tip_all - all points inside tip cell
      ! n_tipcell - number of activated tip cells
      ! n_source  - number of vegf sources
      integer, public :: np_vegf_s, np_tip_s, np_tip_all, n_tipcell = 0, initial_node
      integer, public, allocatable :: vegf_s(:,:), tip_s(:,:), tip_all(:,:)

      !temporary
      real, public :: phi_med, t_med


end module global_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

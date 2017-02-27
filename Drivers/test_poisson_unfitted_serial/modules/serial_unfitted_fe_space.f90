! Copyright (C) 2014 Santiago Badia, Alberto F. Martín and Javier Principe
!
! This file is part of FEMPAR (Finite Element Multiphysics PARallel library)
!
! FEMPAR is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! FEMPAR is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with FEMPAR. If not, see <http://www.gnu.org/licenses/>.
!
! Additional permission under GNU GPL version 3 section 7
!
! If you modify this Program, or any covered work, by linking or combining it 
! with the Intel Math Kernel Library and/or the Watson Sparse Matrix Package 
! and/or the HSL Mathematical Software Library (or a modified version of them), 
! containing parts covered by the terms of their respective licenses, the
! licensors of this Program grant you additional permission to convey the 
! resulting work. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module serial_unfitted_fe_space_names
  use fempar_names
  use serial_unfitted_triangulation_names

  implicit none
# include "debug.i90"
  private

  type, extends(fe_accessor_t) :: unfitted_fe_accessor_t
    private
    class(serial_unfitted_fe_space_t), pointer     :: unfitted_fe_space => NULL()
  contains
    procedure, private :: fe_accessor_create => unfitted_fe_accessor_create
    procedure          :: fe_accessor_free   => unfitted_fe_accessor_free
  end type unfitted_fe_accessor_t

  type :: unfitted_fe_iterator_t
    private
    type(unfitted_fe_accessor_t) :: unfitted_fe_accessor
  contains
    procedure, non_overridable, private :: create       => unfitted_fe_iterator_create
    procedure, non_overridable          :: free         => unfitted_fe_iterator_free
    procedure, non_overridable          :: init         => unfitted_fe_iterator_init
    procedure, non_overridable          :: next         => unfitted_fe_iterator_next
    procedure, non_overridable          :: has_finished => unfitted_fe_iterator_has_finished
    procedure, non_overridable, private :: unfitted_fe_iterator_current
    generic                             :: current      => unfitted_fe_iterator_current
  end type unfitted_fe_iterator_t

  type, extends(serial_fe_space_t) :: serial_unfitted_fe_space_t

    private

    class(serial_unfitted_triangulation_t), pointer :: unfitted_triangulation =>  NULL()

    ! All the machinery for integrating in cut elements
    type(quadrature_t)                  :: quadrature_subelem
    type(tet_lagrangian_reference_fe_t) :: reference_subelem
    type(fe_map_t)                      :: fe_map_subelem
    type(quadrature_t),        allocatable :: cut_quadratures(:)
    type(fe_map_t),            allocatable :: cut_fe_maps(:)
    type(volume_integrator_t), allocatable :: cut_vol_integrators(:)

    contains

    procedure,  private         :: serial_fe_space_create_same_reference_fes_on_all_cells => serial_unfitted_fe_space_create_same_reference_fes_on_all_cells
    procedure                   :: free                                                   => serial_unfitted_fe_space_free
    procedure, non_overridable  :: create_unfitted_fe_iterator                            => serial_unfitted_fe_space_create_unfitted_fe_iterator

    procedure, non_overridable, private :: init_reference_subelem => serial_unfitted_fe_space_init_reference_subelem
    procedure, non_overridable, private :: free_reference_subelem => serial_unfitted_fe_space_free_reference_subelem
    procedure, non_overridable, private :: init_cut_integration   => serial_unfitted_fe_space_init_cut_integration
    procedure, non_overridable, private :: free_cut_integration   => serial_unfitted_fe_space_free_cut_integration

  end type serial_unfitted_fe_space_t

  public :: unfitted_fe_accessor_t
  public :: unfitted_fe_iterator_t
  public :: serial_unfitted_fe_space_t

contains

#include "sbm_unfitted_fe_accessor.i90"
#include "sbm_unfitted_fe_iterator.i90"
#include "sbm_serial_unfitted_fe_space.i90"

end module serial_unfitted_fe_space_names

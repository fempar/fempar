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
module partitioner_input_names
  use serial_names
# include "debug.i90"
  implicit none
  private

  type partitioner_input_t 
     private 
     ! IO parameters
     character(len=:), allocatable :: default_dir_path
     character(len=:), allocatable :: default_prefix
     character(len=:), allocatable :: default_dir_path_out
     character(len=:), allocatable :: default_num_parts

     type(Command_Line_Interface)  :: cli 

     type(ParameterList_t)         :: list

     ! IO parameters
     character(len=256)            :: dir_path
     character(len=256)            :: prefix
     character(len=256)            :: dir_path_out
     integer(ip)                   :: num_parts

   contains
     procedure, non_overridable             :: create       => partitioner_input_create
     procedure, non_overridable, private    :: set_default  => partitioner_input_set_default
     procedure, non_overridable, private    :: add_to_cli   => partitioner_input_add_to_cli
     procedure, non_overridable             :: parse        => partitioner_input_parse
     procedure, non_overridable             :: store_in_fpl => partitioner_input_store_in_fpl
     procedure, non_overridable             :: get_fpl      => partitioner_input_get_fpl
     procedure, non_overridable             :: free         => partitioner_input_free
  end type partitioner_input_t

  public :: partitioner_input_t

contains

  subroutine partitioner_input_create(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    
    call this%free()
    
     ! Initialize Command Line Interface
    call this%cli%init(progname    = 'part',                                                     &
         &        version     = '',                                                                 &
         &        authors     = '',                                                                 &
         &        license     = '',                                                                 &
         &        description =  'FEMPAR driver to part a GiD mesh.', &
         &        examples    = ['part -h  ', 'part -n  ' ])
    
    call this%set_default()
    call this%add_to_cli()
  end subroutine partitioner_input_create
  
  subroutine partitioner_input_set_default(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    ! IO parameters
    this%default_dir_path       = '.'
    this%default_prefix         = 'square'
    this%default_dir_path_out   = '.'
    this%default_num_parts      = '4'
  end subroutine partitioner_input_set_default
  
  subroutine partitioner_input_add_to_cli(this)
    implicit none
    class(partitioner_input_t) , intent(inout) :: this

    ! Locals
    integer(ip) :: error

    ! IO parameters
    call this%cli%add(switch='--dir-path',switch_ab='-d',                              &
         &            help='Directory of the source files',required=.false., act='store', &
         &            def=trim(this%default_dir_path),error=error)
    check(error==0)
    call this%cli%add(switch='--prefix',switch_ab='-p',help='Name of the GiD files',  &
         &            required=.false.,act='store',def=trim(this%default_prefix),error=error) 
    check(error==0)
    call this%cli%add(switch='--dir-path-out',switch_ab='-o',help='Output Directory',&
         &            required=.false.,act='store',def=trim(this%default_dir_path_out),error=error)
    check(error==0)  
    call this%cli%add(switch='--num_parts',switch_ab='-n',help='Number of parts of the mesh',&
         &            required=.true.,act='store',def=trim(this%default_num_parts), error=error)
    check(error==0)  
    
  end subroutine partitioner_input_add_to_cli

  !==================================================================================================
  
  subroutine partitioner_input_parse(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    integer(ip) :: istat
    
    call this%cli%parse(error=istat); check(istat==0)
    
    ! IO parameters
    call this%cli%get(switch='-d',val=this%dir_path    ,error=istat); check(istat==0)
    call this%cli%get(switch='-p',val=this%prefix      ,error=istat); check(istat==0)
    call this%cli%get(switch='-o',val=this%dir_path_out,error=istat); check(istat==0)
    call this%cli%get(switch='-n',val=this%num_parts   ,error=istat); check(istat==0)
  end subroutine partitioner_input_parse  

  subroutine partitioner_input_store_in_fpl(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    integer(ip) :: error

    call this%list%init()
    error = 0
    error = error + this%list%set(key = dir_path_key, value = this%dir_path)
    error = error + this%list%set(key = prefix_key  , value = this%prefix)
    error = error + this%list%set(key = dir_path_out_key, value = this%dir_path_out)
    error = error + this%list%set(key = num_parts_key, value = this%num_parts)
    assert(error == 0)
  end subroutine partitioner_input_store_in_fpl 
  
  subroutine partitioner_input_free(this)
    implicit none
    class(partitioner_input_t), intent(inout) :: this
    if(allocated(this%default_dir_path)) deallocate(this%default_dir_path)              
    if(allocated(this%default_prefix)) deallocate(this%default_prefix)                    
    if(allocated(this%default_dir_path_out)) deallocate(this%default_dir_path_out)
    if(allocated(this%default_num_parts)) deallocate(this%default_num_parts)
    call this%cli%free()
  end subroutine partitioner_input_free

  function partitioner_input_get_fpl(this)
    implicit none
    class(partitioner_input_t), target , intent(in) :: this
    type(ParameterList_t), pointer  :: partitioner_input_get_fpl
    partitioner_input_get_fpl => this%list
  end function partitioner_input_get_fpl

end module partitioner_input_names 


program partitioner
  use serial_names
  use partitioner_input_names
  implicit none
  type(partitioner_input_t)              :: input
  type(ParameterList_t)    , pointer     :: parameters
  type(mesh_t)                           :: gmesh
  type(mesh_distribution_t), allocatable :: distr(:)
  type(mesh_t)             , allocatable :: lmesh(:)
  integer(ip) :: ipart
  integer(ip) :: error
  logical     :: is_present

  call fempar_init()
  call input%create()
  call input%parse()
  call input%store_in_fpl()
  parameters => input%get_fpl()

  is_present =  parameters%isPresent(key = num_parts_key )
  assert(is_present)
  error = parameters%get(key = num_parts_key , value = ipart)
  check(error==0)

  call gmesh%read(parameters)
  call gmesh%write_file_for_postprocess(parameters)

  ! Set partition parameters
  error = 0
  error = error + parameters%set(key = strategy_key, value = part_kway)
  error = error + parameters%set(key = debug_key   , value = 0)
  error = error + parameters%set(key = metis_option_debug_key  , value =  2)
  error = error + parameters%set(key = metis_option_ufactor_key, value = 30)
  error = error + parameters%set(key = metis_option_minconn_key, value =  0)
  error = error + parameters%set(key = metis_option_contig_key , value =  1)
  error = error + parameters%set(key = metis_option_ctype_key  , value = METIS_CTYPE_SHEM) ! METIS_CTYPE_RM
  error = error + parameters%set(key = metis_option_iptype_key , value = METIS_IPTYPE_EDGE)
  check(error==0)
  call gmesh%create_distribution (parameters, distr, lmesh)

  ! Write partition info
  call mesh_distribution_write_files           ( parameters, distr )
  call mesh_distribution_write_for_postprocess ( parameters, gmesh, distr )

  ! Write local meshes
  call mesh_write_files                 ( parameters, lmesh )
  call mesh_write_files_for_postprocess ( parameters, lmesh )

  ! Deallocate partition objects
  do ipart=1,size(distr)
     call distr(ipart)%free()
     call lmesh(ipart)%free
  end do
  deallocate (distr)
  deallocate (lmesh)

  call gmesh%free()

  call fempar_finalize()

end program partitioner

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
module par_scalar_matrix_names
  ! Serial modules
  use types_names
  use memor_names
  use serial_scalar_matrix_names
  use stdio_names
  use graph_names
#ifdef memcheck
  use iso_c_binding
#endif

  ! Parallel modules
  use par_environment_names
  use par_context_names
  use par_scalar_array_names
  use psb_penv_mod_names
  use dof_distribution_names

  ! Abstract types
  use vector_names
  use operator_names
  use matrix_names
  use vector_space_names

  implicit none
# include "debug.i90"

  private

  integer(ip), parameter :: start           = 0
  integer(ip), parameter :: properties_set  = 1
  integer(ip), parameter :: environment_set = 2
  integer(ip), parameter :: created         = 3
  
  
  type, extends(matrix_t) :: par_scalar_matrix_t
     integer(ip) :: state = start
  
     type( serial_scalar_matrix_t ) :: serial_scalar_matrix

     type(dof_distribution_t), pointer :: &
          dof_dist_domain => NULL()

     type(dof_distribution_t), pointer :: &
          dof_dist_range => NULL()

     type(par_environment_t), pointer :: &
          p_env => NULL()
   contains
   
     generic  :: create   => par_scalar_matrix_create_square, &
                             par_scalar_matrix_create_rectangular
     procedure, private :: par_scalar_matrix_create_square
     procedure, private :: par_scalar_matrix_create_rectangular
   
     generic :: set_properties => par_scalar_matrix_set_properties_square, &
                                  par_scalar_matrix_set_properties_rectangular
     procedure, private :: par_scalar_matrix_set_properties_square
     procedure, private :: par_scalar_matrix_set_properties_rectangular
                           
     procedure :: set_environment => par_scalar_matrix_set_environment
   
     generic  :: set_dof_distribution   => par_scalar_matrix_set_dof_distribution_square, &
                             par_scalar_matrix_set_dof_distribution_rectangular
     procedure, private :: par_scalar_matrix_set_dof_distribution_square
     procedure, private :: par_scalar_matrix_set_dof_distribution_rectangular
     
     procedure  :: get_graph    => par_scalar_matrix_get_graph
     procedure  :: return_graph => par_scalar_matrix_return_graph 
     
     procedure  :: allocate            => par_scalar_matrix_allocate
     procedure  :: print_matrix_market => par_scalar_matrix_print_matrix_market						 
     procedure  :: init                => par_scalar_matrix_init
     procedure  :: apply               => par_scalar_matrix_apply
     procedure  :: free_in_stages      => par_scalar_matrix_free_in_stages
     
     procedure, private :: create_vector_spaces => par_scalar_matrix_create_vector_spaces
  end type par_scalar_matrix_t

  ! Types
  public :: par_scalar_matrix_t


!***********************************************************************
! Allocatable arrays of type(par_scalar_matrix_t)
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(par_scalar_matrix_t)
# define var_size 80
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"
  
  !=============================================================================
  subroutine par_scalar_matrix_create_square(this,symmetric_storage,is_symmetric,sign,dof_dist,p_env)
    implicit none
    class(par_scalar_matrix_t)          ,intent(inout) :: this
    logical                             ,intent(in)    :: symmetric_storage
    logical                             ,intent(in)    :: is_symmetric
    integer(ip)                         ,intent(in)    :: sign
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist
    type(par_environment_t) , target    ,intent(in)    :: p_env

    assert ( this%state == start )
    this%dof_dist_domain      => dof_dist 
    this%dof_dist_range       => dof_dist
    this%p_env                => p_env
    if(this%p_env%p_context%iam>=0) then
      call this%serial_scalar_matrix%create(dof_dist%nl,symmetric_storage,is_symmetric,sign)
    end if  
    call this%create_vector_spaces()
    this%state = created
  end subroutine par_scalar_matrix_create_square
  
  !=============================================================================
  subroutine par_scalar_matrix_create_rectangular(this,dof_dist_domain,dof_dist_range,p_env)
    implicit none
    class(par_scalar_matrix_t)          ,intent(inout) :: this
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist_domain
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist_range
    type(par_environment_t) , target    ,intent(in)    :: p_env

    assert ( this%state == start )
    this%dof_dist_domain => dof_dist_domain 
    this%dof_dist_range  => dof_dist_range
    this%p_env           => p_env
    if(this%p_env%p_context%iam>=0) then
      call this%serial_scalar_matrix%create(dof_dist_domain%nl,dof_dist_range%nl)
    end if
    call this%create_vector_spaces()
    this%state = created
  end subroutine par_scalar_matrix_create_rectangular
  
  !=============================================================================
  subroutine par_scalar_matrix_set_properties_square(this,symmetric_storage,is_symmetric,sign)
    implicit none
    class(par_scalar_matrix_t)          ,intent(inout) :: this
    logical                             ,intent(in)    :: symmetric_storage
    logical                             ,intent(in)    :: is_symmetric
    integer(ip)                         ,intent(in)    :: sign
    
    assert ( this%state == start )
    if(this%p_env%p_context%iam>=0) then
      call this%serial_scalar_matrix%set_properties(symmetric_storage,is_symmetric,sign)
    end if  
    this%state = properties_set
  end subroutine par_scalar_matrix_set_properties_square
  
  !=============================================================================
  subroutine par_scalar_matrix_set_properties_rectangular(this)
    implicit none
    class(par_scalar_matrix_t), intent(inout) :: this  
    assert ( this%state == start )
    if(this%p_env%p_context%iam>=0) then
      call this%serial_scalar_matrix%set_properties()
    end if  
    this%state = properties_set
  end subroutine par_scalar_matrix_set_properties_rectangular
  
  !=============================================================================
  subroutine par_scalar_matrix_set_environment(this,p_env)
    implicit none
    class(par_scalar_matrix_t)        , intent(inout) :: this
    class(par_environment_t)  , target, intent(in)    :: p_env
    assert ( this%state == properties_set )
    this%p_env => p_env
    this%state = environment_set
  end subroutine par_scalar_matrix_set_environment
  
  !=============================================================================
  subroutine par_scalar_matrix_set_dof_distribution_square(this,dof_dist)
    implicit none
    class(par_scalar_matrix_t)          ,intent(inout) :: this
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist

    assert ( this%state == environment_set )
    this%dof_dist_domain      => dof_dist 
    this%dof_dist_range       => dof_dist
    call this%create_vector_spaces()
    this%state = created
  end subroutine par_scalar_matrix_set_dof_distribution_square
  
  !=============================================================================
  subroutine par_scalar_matrix_set_dof_distribution_rectangular(this,dof_dist_domain,dof_dist_range)
    implicit none
    class(par_scalar_matrix_t)          ,intent(inout) :: this
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist_domain
    type(dof_distribution_t), target    ,intent(in)    :: dof_dist_range

    assert ( this%state == environment_set )
    this%dof_dist_domain => dof_dist_domain 
    this%dof_dist_range  => dof_dist_range
    call this%create_vector_spaces()
    this%state = created
  end subroutine par_scalar_matrix_set_dof_distribution_rectangular
  
  subroutine par_scalar_matrix_create_vector_spaces (this)
    implicit none
    class(par_scalar_matrix_t), intent(inout) :: this
    type(par_scalar_array_t) :: range_vector
    type(par_scalar_array_t) :: domain_vector
    type(vector_space_t), pointer :: range_vector_space
    type(vector_space_t), pointer :: domain_vector_space
    
    call range_vector%create(this%dof_dist_range,this%p_env)
    call domain_vector%create(this%dof_dist_domain,this%p_env)
    range_vector_space => this%get_range_vector_space()
    call range_vector_space%create(range_vector)
    domain_vector_space => this%get_domain_vector_space()
    call domain_vector_space%create(domain_vector)
    
    call range_vector%free()
    call domain_vector%free()
  end subroutine par_scalar_matrix_create_vector_spaces
  
  subroutine par_scalar_matrix_allocate(this)
    implicit none
    class(par_scalar_matrix_t), intent(inout) :: this
    assert ( this%state == created )
    if(this%p_env%p_context%iam>=0) then
       call this%serial_scalar_matrix%allocate()
    end if
  end subroutine par_scalar_matrix_allocate
  
  !=============================================================================
  function par_scalar_matrix_get_graph (this)
    implicit none
    class(par_scalar_matrix_t), target, intent(inout) :: this
    type(graph_t), pointer                            :: par_scalar_matrix_get_graph
    assert ( this%state == created )
    if(this%p_env%p_context%iam>=0) then
      par_scalar_matrix_get_graph =>  this%serial_scalar_matrix%get_graph()
    else
      nullify(par_scalar_matrix_get_graph)
    end if  
  end function par_scalar_matrix_get_graph
  
  !=============================================================================
  subroutine par_scalar_matrix_return_graph (this,graph)
    implicit none
    class(par_scalar_matrix_t), intent(inout) :: this
    type(graph_t), pointer, intent(inout)        :: graph
    if(this%p_env%p_context%iam>=0) then
      call this%serial_scalar_matrix%return_graph(graph)
    end if  
  end subroutine par_scalar_matrix_return_graph

  !=============================================================================
  subroutine par_scalar_matrix_free_in_stages(this, action)
    implicit none
    class(par_scalar_matrix_t), intent(inout) :: this
    integer(ip)               , intent(in)    :: action

    if ( this%state == start .or. this%state == properties_set ) return

    if ( action == free_clean ) then
       if ( this%state == created ) then
          if(this%p_env%p_context%iam>=0) then
            call this%serial_scalar_matrix%free_in_stages(action)
          end if  
          nullify ( this%dof_dist_domain )
          nullify ( this%dof_dist_range )
          nullify ( this%p_env )
          call this%free_vector_spaces()
          this%state = start
       end if   
    else if ( action == free_symbolic_setup ) then
        if(this%p_env%p_context%iam>=0) then
          call this%serial_scalar_matrix%free_in_stages(action)
        end if 
    else if ( action == free_numerical_setup ) then
        if(this%p_env%p_context%iam>=0) then
          call this%serial_scalar_matrix%free_in_stages(action)
        end if     
    end if
    
  end subroutine par_scalar_matrix_free_in_stages
  

  subroutine par_scalar_matrix_print_matrix_market ( this, dir_path, prefix )
    implicit none
    class(par_scalar_matrix_t), intent(in) :: this
    character (*)            , intent(in) :: dir_path
    character (*)            , intent(in) :: prefix
    integer         :: iam, num_procs, lunou
    integer(ip)     :: j, ndigs_iam, ndigs_num_procs, id_map
    character(256)  :: name 
    character(256)  :: zeros
    character(256)  :: part_id
    
    assert ( this%state == created )
    
    if(this%p_env%p_context%iam<0) return

    name = trim(prefix) // '.par_matrix' // '.mtx'

    ! Get context info
    call par_context_info ( this%p_env%p_context, iam, num_procs )

    ! Form the file_path of the partition object to be read
    iam = iam + 1 ! Partition identifers start from 1 !!

    ndigs_num_procs = count_digits_par_matrix (num_procs)
    zeros = ' '   
    ndigs_iam = count_digits_par_matrix ( iam )

    ! write(*,*) ndgs_num_procs, ndigs_iam DBG

    do j=1,  ndigs_num_procs - ndigs_iam
       zeros (j:j) = '0'
    end do
    part_id = ch(iam)

    ! Read partition data from path_file file
    lunou =  io_open (trim(dir_path) // '/' // trim(name) // '.' // trim(zeros) // trim(part_id), 'write')

    call this%serial_scalar_matrix%print_matrix_market(lunou)

    call io_close (lunou)

  end subroutine par_scalar_matrix_print_matrix_market

  function count_digits_par_matrix ( i )
    implicit none
    ! Parameters
    integer(ip), intent(in) :: i 
    integer(ip)             :: count_digits_par_matrix
    ! Locals   
    integer(ip)             :: x 
    x = i 
    if (x < 0) x = -x;
    count_digits_par_matrix = 1;
    x = x/10;
    do while( x > 0)
       count_digits_par_matrix = count_digits_par_matrix + 1
       x = x/10;
    end do
  end function count_digits_par_matrix

  subroutine par_scalar_matrix_init (this, alpha)
    implicit none
    ! Parameters 
    class(par_scalar_matrix_t), intent(inout) :: this
    real(rp)                  , intent(in)    :: alpha
    assert ( this%state == created )
    if(this%p_env%p_context%iam<0) return
    call this%serial_scalar_matrix%init(alpha)
  end subroutine par_scalar_matrix_init

  subroutine par_scalar_matrix_apply_concrete(a,x,y)
    implicit none
    ! Parameters
    type(par_scalar_matrix_t) , intent(in)    :: a
    type(par_scalar_array_t) , intent(in)    :: x
    type(par_scalar_array_t) , intent(inout) :: y
    if(a%p_env%p_context%iam<0) return
    call a%serial_scalar_matrix%apply(x%serial_scalar_array, y%serial_scalar_array)
    call y%comm()
  end subroutine par_scalar_matrix_apply_concrete

  ! op%apply(x,y) <=> y <- op*x
  ! Implicitly assumes that y is already allocated
  subroutine par_scalar_matrix_apply(op,x,y) 
    implicit none
    class(par_scalar_matrix_t), intent(in)    :: op
    class(vector_t) , intent(in)    :: x
    class(vector_t) , intent(inout) :: y 

    assert ( op%state == created )
    
    call op%abort_if_not_in_domain(x)
    call op%abort_if_not_in_range(y)
    
    call x%GuardTemp()
    select type(x)
       class is (par_scalar_array_t)
       select type(y)
          class is(par_scalar_array_t)
          call par_scalar_matrix_apply_concrete(op, x, y)
       end select
    end select
    call x%CleanTemp()
  end subroutine par_scalar_matrix_apply

end module par_scalar_matrix_names

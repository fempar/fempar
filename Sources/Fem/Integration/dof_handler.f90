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
module dof_handler_names
  use types
  use memor

  implicit none
# include "debug.i90"
  private
  ! Data needed to go from geometrical info (import, partition, mesh) to 
  ! a dof info (dof_import, dof_partition, dof_mesh)

 !  type dof_block

!       integer(ip) :: nvars                    ! Total number of variables per block

!      integer(ip), allocatable ::     &
! !          ndofs,                     &       ! Total number of dofs (nprob)
!           nvars_prob(:)                       ! Number of vars per problem vars_prob(nprob)

!      type(array_ip2) ::              &
!           vars_prob(:)                        ! Variables per problem

!   end type dof_block

  type physical_problem
     integer(ip)        ::           &
          nvars                              ! Number of different variables
     integer(ip), allocatable ::     &
          l2g_var(:)                         ! Order chosen for variables (size nvars)
     integer(ip)        ::           &
          problem_code                       ! An internal code that defines a problem in FEMPAR
  end type physical_problem

  type dof_handler

     integer(ip) ::     &
          nblocks,                   &       ! Number of blocks
          nprobs,                     &       ! Number of problems
          nvars_global                      ! Total number of different physical variables
          !nvars_prob                         ! Number of physical variables per problem

     integer(ip), allocatable ::     &
!          l2g_vars                   &       ! Local (in block) to global Id of variable
!          phys_prob(:),              &       ! List of physical problems (nprob) 
          dof_coupl(:,:),            &        ! Dof_coupling(nvar,nvar) for avoiding allocation & assembly of zero blocks
          vars_block(:)                       ! Parameter per unknown (size nvars)

!     type(dof_handler), allocatable :: dof_blocks(:)

     type(physical_problem), allocatable :: problems(:)
!     type(fem_blocks)  ::           &
!          blocks                             ! blocks ordering

  end type dof_handler

  ! type physical_problem_pointer
  !    type(physical_problem)          , pointer :: p => NULL()
  ! end type physical_problem_pointer

  ! Types
  public :: dof_handler, physical_problem!, physical_problem_pointer

  ! Functions
  public ::  dof_handler_create, dof_handler_print !, dof_handler_fill, dof_handler_free

contains
  subroutine dof_handler_create( dhand, nblocks, nprobs, nvars_prob )
    implicit none
    ! Parameters
    type(dof_handler), intent(inout)          :: dhand
    integer(ip), intent(in)                   :: nblocks, nprobs, nvars_prob(nprobs)

    integer(ip) :: istat, i, j, nvars_global

    

    allocate( dhand%problems(nprobs), stat=istat)
    check( istat==0 )
    
    nvars_global = 0
    do i = 1, nprobs
       dhand%problems(i)%nvars = nvars_prob(i)
       call memalloc ( dhand%problems(i)%nvars, dhand%problems(i)%l2g_var, __FILE__, __LINE__ ) 
       dhand%problems(i)%l2g_var(1:nvars_prob(i)) = (/ (j, j=1,nvars_prob(i)) /) !1:dhand%problems(i)%nvars /)
       !dhand%problems(i)%problem_code = problem_code(i)
       nvars_global = nvars_global + nvars_prob(i)
    end do

    dhand%nblocks = nblocks
    dhand%nprobs = nprobs

    dhand%nvars_global = nvars_global

    call memalloc ( dhand%nvars_global, dhand%nvars_global, dhand%dof_coupl, __FILE__, __LINE__ ) 
    dhand%dof_coupl = 1
    
    call memalloc ( dhand%nvars_global, dhand%vars_block, __FILE__, __LINE__ ) 
    dhand%vars_block = 1
    
    
  end subroutine dof_handler_create


  subroutine dof_handler_print(  dhand, lunou )
    implicit none
    integer(ip)      , intent(in)           :: lunou
    type(dof_handler), intent(in)           :: dhand
    
    integer(ip) :: iprob
    
    write (lunou, '(a)')     '*** begin dof handler data structure ***'

    write (lunou,*)     'Number of problems: ',   dhand%nprobs
    write (lunou,*)     'Number of blocks: '  ,   dhand%nblocks
    write (lunou,*)     'Number of variables: ',  dhand%nvars_global

    do iprob = 1, dhand%nprobs
              
       write (lunou, '(a)')     '*** physical problem ',iprob ,'***'
       write (lunou,*)     'Number of variables of problem ',iprob, ' :' ,  dhand%problems(iprob)%nvars
       write (lunou,*)     'Local to global (of variables) for problem ',iprob, ' :' ,  dhand%problems(iprob)%l2g_var
       !write (lunou,*)     'Number of variables of problem ',iprob, ' :' ,  dhand%problems(iprob)%problem_code
       
    end do

    write (lunou,*)     'Block of every variable: ',    dhand%vars_block
    write (lunou,*)     'Coupling flag between dofs: ', dhand%dof_coupl
  
  end subroutine dof_handler_print
    
 
!  Problem codes in types.f90   
!  integer(ip), parameter :: nsi_code=1, ela_code=2, cdr_code=3, adr_code=4 
!  integer(ip), parameter :: imh_code=5, dcy_code=6, mss_code=7, lap_code=8 
 

end module dof_handler_names

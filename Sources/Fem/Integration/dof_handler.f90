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
use types_names
use memor_names
  use problem_names
  use array_names

  implicit none
# include "debug.i90"
  private



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

     !type(p_physical_problem), allocatable :: problems(:)
     type(discrete_problem_pointer), allocatable :: problems(:)


     ! Auxiliary arrays
     integer(ip), allocatable     :: g2l_vars(:,:)
     type(array_ip1), allocatable :: prob_block(:,:)

   contains
     procedure :: set_problem
     procedure :: create => dof_handler_create
  end type dof_handler

  ! Types
  public :: dof_handler !, physical_problem!, physical_problem_pointer

  ! Functions
  public ::  dof_handler_print, dof_handler_free !, dof_handler_fill

contains

  subroutine dof_handler_create( dhand, nblocks, nprobs, nvars_global, vars_block, dof_coupl )
    implicit none
    ! Parameters
    class(dof_handler), intent(inout)          :: dhand
    integer(ip), intent(in)                   :: nblocks, nprobs, nvars_global 
    integer(ip), intent(in), optional         :: vars_block(:), dof_coupl(:,:)

    integer(ip) :: istat, i, j, iprob, iblock



    allocate( dhand%problems(nprobs), stat=istat)
    check( istat==0 )

    dhand%nblocks = nblocks
    dhand%nprobs = nprobs
    dhand%nvars_global = nvars_global
    
    if (present(dof_coupl)) then
       assert ( size(dof_coupl,1) == nvars_global .and. size(dof_coupl,2) == nvars_global ) 
       dhand%dof_coupl = dof_coupl
    else
       call memalloc ( dhand%nvars_global, dhand%nvars_global, dhand%dof_coupl, __FILE__, __LINE__ ) 
       dhand%dof_coupl = 1
    end if
    
    if (present(vars_block)) then
       assert ( size(vars_block) == nvars_global )
       dhand%vars_block = vars_block
    else
       call memalloc ( dhand%nvars_global, dhand%vars_block, __FILE__, __LINE__ ) 
       dhand%vars_block = 1
    end if

    ! Global to local of variables per problem
    call memalloc( nvars_global, dhand%nprobs, dhand%g2l_vars, __FILE__, __LINE__ )
    
    call memalloc( dhand%nblocks, dhand%nprobs, dhand%prob_block, __FILE__, __LINE__ )

  end subroutine dof_handler_create


  subroutine dof_handler_print(  dhand, lunou )
    implicit none
    integer(ip)      , intent(in)           :: lunou
    type(dof_handler), intent(in)           :: dhand

    integer(ip) :: iprob, count, iblock

    write (lunou, '(a)')     '*** begin dof handler data structure ***'

    write (lunou,*)     'Number of problems: ',   dhand%nprobs
    write (lunou,*)     'Number of blocks: '  ,   dhand%nblocks
    write (lunou,*)     'Number of variables: ',  dhand%nvars_global

    do iprob = 1, dhand%nprobs

       write (lunou,*)     '*** physical problem ',iprob, ' ***'
       write (lunou,*)     'Number of variables of problem ',iprob, ' :' ,  dhand%problems(iprob)%p%nvars
       write (lunou,*)     'Local to global (of variables) for problem ',iprob, ' :' ,  dhand%problems(iprob)%p%l2g_var
       !write (lunou,*)     'Number of variables of problem ',iprob, ' :' ,  dhand%problems(iprob)%problem_code

    end do


    do iblock = 1, dhand%nblocks  
       do iprob = 1, dhand%nprobs 
          write (lunou,*) 'prob_block array iblock ',iblock,' problem ',iprob,' :', dhand%prob_block(iblock,iprob)%nd1
          write (lunou,*) 'prob_block array iblock ',iblock,' problem ',iprob,' :', dhand%prob_block(iblock,iprob)%a
         
       end do
    end do

    write (lunou,*)     'Block of every variable: ',    dhand%vars_block
    write (lunou,*)     'Coupling flag between dofs: ', dhand%dof_coupl

  end subroutine dof_handler_print

  subroutine dof_handler_free(  dhand )
    implicit none
    type(dof_handler), intent(inout)           :: dhand

    integer(ip) :: i,j

    call memfree( dhand%dof_coupl, __FILE__, __LINE__ )
    call memfree( dhand%vars_block, __FILE__, __LINE__ )

    call memfree( dhand%g2l_vars, __FILE__, __LINE__ )

    do i = 1, size( dhand%prob_block,1 )
       do j = 1, size( dhand%prob_block,2 )
          call memfree ( dhand%prob_block(i,j)%a, __FILE__, __LINE__ )
       end do
    end do
    call memfree ( dhand%prob_block, __FILE__, __LINE__ )

  end subroutine dof_handler_free
    

  subroutine set_problem ( dhand, iprob, prob )
    implicit none
    class(dof_handler), intent(inout)          :: dhand
    integer(ip), intent(in) :: iprob
    class(discrete_problem), target, intent(in) :: prob
    
    integer(ip) :: l_var, iblock, count

    dhand%problems(iprob)%p => prob
   
    do l_var = 1, dhand%problems(iprob)%p%nvars
       dhand%g2l_vars( dhand%problems(iprob)%p%l2g_var(l_var), iprob) = l_var
    end do


    do iblock = 1, dhand%nblocks  
       ! variables of a given problem for current block (prob_block)
       count = 0
       do l_var = 1, dhand%problems(iprob)%p%nvars
          if ( dhand%vars_block(dhand%problems(iprob)%p%l2g_var(l_var)) == iblock ) then
             count = count + 1 
          end if
       end do
       call array_create( count, dhand%prob_block(iblock,iprob))
       count = 0 
       do l_var = 1, dhand%problems(iprob)%p%nvars
          if ( dhand%vars_block(dhand%problems(iprob)%p%l2g_var(l_var)) == iblock ) then
             count = count + 1 
             dhand%prob_block(iblock,iprob)%a(count) = l_var !dhand%problems(iprob)%p%l2g_var(l_var)
          end if
       end do

    end do

  end subroutine set_problem

 
!  Problem codes in types.f90   
!  integer(ip), parameter :: nsi_code=1, ela_code=2, cdr_code=3, adr_code=4 
!  integer(ip), parameter :: imh_code=5, dcy_code=6, mss_code=7, lap_code=8 
 

end module dof_handler_names

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
module mem_base_names
  !-----------------------------------------------------------------------
  !
  !  This module contains memory allocation functions.
  !
  !  TODO:
  !  * This module is not thread-safe.
  !
  !  * Understand the behavior of character(*) and make a decision
  !    about it. 
  !
  !    It seems that character*(*) is the method of choice for I/O, 
  !    as explained on the following document 
  !    http://www.ibiblio.org/pub/languages/fortran/ch2-4.html
  !
  !    However, it is an obsolent feature in Fortran 95. So?
  !
  !  * Provide a function to assign units to write
  !
  !-----------------------------------------------------------------------
!$ use omp_lib
use types_names
  use hash_table_names ! Only used ifdef memcheck but this avoids passing 
                       ! the definition to configure
#ifdef memcheck
  use iso_c_binding
#endif
#ifdef exception
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
#endif
  implicit none
  private
#ifdef memcheck
  interface
     subroutine cptr2igp(cptr,longip) bind(c,name='cptr2igp')
       use iso_c_binding
       implicit none
       type(c_ptr), value, intent(in)  :: cptr 
       integer(c_long)   , intent(out) :: longip
     end subroutine cptr2igp
  end interface
#endif

  integer(ip) , save ::    &
       lumem = 0,          &  ! Logical unit to write memory evolution
       luerr = 0,          &  ! Logical unit to write errors
       meall = 0,          &  ! Number of memory allocations
       medea = 0              ! Number of memory deallocations

  integer(imp), save ::    &
       memax = 0,          &  ! Maximum memory allocation in bytes
       mecur = 0              ! Current memory allocation in bytes

  real(rp), save :: nan       ! An undefined value, initialized to nan 
                              ! if exception macro is defined

  !***********************************************************************
  ! Specific purpose hash table specification (only if memcheck is defined)
  ! using a special purpose type mem_data_t to store where allocations have
  ! been performed.
  !***********************************************************************
#ifdef memcheck
  type mem_data_t
     character(512) :: file=' '
     integer(ip)    :: line=0
  end type mem_data_t

#include "debug.i90"

#define hash_table hash_table_mem
#define hash_node  hash_node_mem
#define key_type   integer(igp)
#define key_size   igp
#define val_type   type(mem_data_t)
#define print_key_val     write(*,*) list%key, trim(list%val%file), list%val%line
#define put_hash_node     put_hash_node_mem   
#define get_hash_node     get_hash_node_mem   
#define del_hash_node     del_hash_node_mem   
#define free_hash_node    free_hash_node_mem  
#define print_hash_node   print_hash_node_mem 
#define init_hash_table   init_hash_table_mem 
#define put_hash_table    put_hash_table_mem  
#define get_hash_table    get_hash_table_mem  
#define del_hash_table    del_hash_table_mem  
#define free_hash_table   free_hash_table_mem 
#define print_hash_table  print_hash_table_mem
#define status_hash_table  status_hash_table_mem
#include "hash_table_header.i90"

  type(hash_table_mem) :: mem_db
#endif

  public :: fempar_memmax, fempar_memcur, meminit, memstatus,  &
       &    memsum, memsub, mem_alloc_error, mem_dealloc_error, mem_status_error

#ifdef memcheck
  public ::  mem_db_put, mem_db_del
#endif

contains
#ifdef memcheck
#include "hash_table_body.i90"
#endif

  !***********************************************************************
  !***********************************************************************
  ! Count, write and error functions
  !***********************************************************************
  !***********************************************************************
  subroutine memsum(lbyts)
    implicit none
    integer(imp), intent(in)    :: lbyts

!$OMP CRITICAL (memor_lock)
    mecur = mecur+lbyts
    memax = max(memax,mecur)
    meall = meall+1
!$OMP END CRITICAL (memor_lock)

  end subroutine memsum

  !-----------------------------------------------------------------------
  subroutine memsub(lbyts)
    implicit none
    integer(imp), intent(in)    :: lbyts

!$OMP CRITICAL (memor_lock)
    mecur = mecur-lbyts
    medea = medea+1
!$OMP END CRITICAL (memor_lock)

  end subroutine memsub

  !-----------------------------------------------------------------------
  subroutine mem_alloc_error(istat,lbyts,file,line)
    implicit none
    integer(ip)  , intent(in), optional :: istat
    integer(imp) , intent(in), optional :: lbyts
    character*(*), intent(in), optional :: file                     ! Calling file
    integer(ip)  , intent(in), optional :: line                     ! Calling line
    integer(ip)    :: luout
    character(20)  :: lsize,lstat
    luout = 6
    lsize = ' unknown'
    lstat = ' unknown'
    if(luerr > 0) luout = luerr
    if(present(lbyts)) write(lsize,'(i20)') lbyts
    if(present(istat)) write(lstat,'(i20)') istat
    write(luout,'(a)') '[Fempar Fatal Error] ***Memory allocation failed.'
    write(luout,'(a)') 'Error code: '// trim(lstat)
    if(present(file)) then
       if(present(line)) then
          write(luout,'(a,i10)') 'Called from file '// trim(file)//', line ',line
       else
          write(luout,'(a)') 'Called from file '// trim(file)
       end if
    end if
    write(luout,'(a)') 'Requested size: ' // trim(lsize) // ' bytes'
    ! Fatal error. Stop the program. 
    call runend
  end subroutine mem_alloc_error

  subroutine mem_dealloc_error(istat,file,line)
    implicit none
    integer(ip)  , intent(in), optional :: istat
    character*(*), intent(in), optional :: file                     ! Calling file
    integer(ip)  , intent(in), optional :: line                     ! Calling line
    integer(ip)    :: luout
    character(20)  :: lsize,lstat
    luout = 6
    lstat = ' unknown'
    if(luerr > 0) luout = luerr
    if(present(istat)) write(lstat,'(i20)') istat
    write(luout,'(a)') '[Fempar Fatal Error] ***Memory deallocation failed.'
    write(luout,'(a)') 'Error code: '// trim(lstat)
    if(present(file)) then
       if(present(line)) then
          write(luout,'(a,i10)') 'Called from file '// trim(file)//', line ',line
       else
          write(luout,'(a)') 'Called from file '// trim(file)
       end if
    end if
    ! Fatal error. Stop the program. 
    call runend
  end subroutine mem_dealloc_error

  subroutine mem_status_error(file,line)
    implicit none
    character*(*), intent(in), optional :: file                     ! Calling file
    integer(ip)  , intent(in), optional :: line                     ! Calling line
    integer(ip)    :: luout
    character(20)  :: lsize,lstat
    luout = 6
    if(luerr > 0) luout = luerr
    write(luout,'(a)') '[Fempar Fatal Error] ***Attempting to (re)allocate an (un)allocated variable.'
    if(present(file)) then
       if(present(line)) then
          write(luout,'(a,i10)') 'Called from file '// trim(file)//', line ',line
       else
          write(luout,'(a)') 'Called from file '// trim(file)
       end if
    end if
    ! Fatal error. Stop the program. 
    call runend
  end subroutine mem_status_error

  !-----------------------------------------------------------------------
  subroutine fempar_memmax (maxmem)
    implicit none
    integer(imp), intent(out)    :: maxmem
    maxmem = memax
  end subroutine fempar_memmax

  !-----------------------------------------------------------------------
  subroutine fempar_memcur (curmem)
    implicit none
    integer(imp), intent(out)    :: curmem
    curmem = mecur
  end subroutine fempar_memcur

  !-----------------------------------------------------------------------
  ! An initialization routine, ieee stuff can be handled here.
  ! Error handling behavior (stop or return) defined here
  ! Logical unit to write messages defined here
  subroutine meminit
    implicit none
#ifdef exception
    nan = ieee_value( nan, ieee_signaling_nan)
    call ieee_set_halting_mode(ieee_invalid, .true.)   ! Enable trapping of invalid.
#endif
#ifdef memcheck
    call mem_db%init
#endif
  end subroutine meminit

#ifdef memcheck
  !-----------------------------------------------------------------------
  ! Memory allocation/deallocation database
  subroutine mem_db_put(varptr,file,line)
    implicit none
    type(c_ptr)  , intent(in)           :: varptr
    character*(*), intent(in), optional :: file              ! Calling file
    integer(ip)  , intent(in), optional :: line              ! Calling line
    integer(igp)   :: key
    type(mem_data_t) :: val, old_val
    integer(ip)    :: istat
    !key = transfer(varptr, 0_igp)
    call cptr2igp(varptr,key)

    if(present(file)) val%file = file
    if(present(line)) val%line = line
    !write(*,*) 'MEM_DB Storing ', key, varptr, file, line
    call mem_db%put(key, val, istat)
    if(istat/=now_stored) then
       write(*,*) 'An error ocurred storing allocation in mem_db: ', stat(istat)
       write(*,*) 'Called from '//file , line
       !write(*,*) 'key: ', key, 'varptr: ', varptr ! DBG
       if(istat==was_stored) then
          call mem_db%get(key, old_val, istat)
          if(istat==key_found) then
             write(*,*) 'Stored allocation was performed from:'//old_val%file, old_val%line
          else
             write(*,*) 'Error getting stored allocation', istat
          end if
       end if
       !call memstatus
       !assert(.false.) ! istat = 1/0
       !stop
    end if
  end subroutine mem_db_put

  subroutine mem_db_del(varptr,file,line)
    implicit none
    type(c_ptr)  , intent(in)           :: varptr
    character*(*), intent(in), optional :: file              ! Calling file
    integer(ip)  , intent(in), optional :: line              ! Calling line
    integer(igp)   :: key
    integer(ip)    :: istat
    !key = transfer(varptr, 0_igp)
    call cptr2igp(varptr,key)

    call mem_db%del(key, istat)
    if(istat/=deleted) then
       write(*,*) 'An error ocurred deleting allocation from mem_db: ',stat(istat)
       if(present(file)) then
          if(present(line)) then
             write(*,*) 'Called from '//file , line
          else
             write(*,*) 'Called from '//file 
          end if
       end if
       !write(*,*) 'varptr: ', varptr
       !write(*,*) 'key: ', key
       !call mem_db%status(0,istat)
    end if
  end subroutine mem_db_del
#endif

  subroutine memstatus
    implicit none
    write(*,*) '====================================================='
    write(*,*) 'Current memory usage:',mecur
    write(*,*) 'Maximum memory usage:',memax
#ifdef memcheck
    write(*,*) 'Currently allocated variables are:'
    call mem_db%print
#endif
    write(*,*) '====================================================='
  end subroutine memstatus

end module mem_base_names

!***********************************************************************
!***********************************************************************
! Allocatable routines
!***********************************************************************
!***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc
!***********************************************************************
! integer(ip)
!***********************************************************************
module mem_ip_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type integer(ip)
# define var_size ip
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module mem_ip_allocatable_names
!***********************************************************************
! integer(ieep)
!***********************************************************************
module mem_ieep_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type integer(ieep)
# define var_size ieep
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module mem_ieep_allocatable_names
!***********************************************************************
! integer(igp)
!***********************************************************************
module mem_igp_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   integer(igp)
# define var_size   igp
# define bound_kind igp
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module mem_igp_allocatable_names
!***********************************************************************
! real(rp)
!***********************************************************************
module mem_rp_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   real(rp)
# define var_size   rp
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
#ifdef exception
#define exception_real
#endif
# include "mem_body.i90"
#ifdef exception
#undef exception_real
#endif
end module mem_rp_allocatable_names
!***********************************************************************
! logical
!***********************************************************************
module mem_lg_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type logical
# define var_size 1 
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc,  memfree, memmovealloc
contains
# include "mem_body.i90"
end module mem_lg_allocatable_names
!***********************************************************************
! mixed integer(ip)+integer(igp) specialization
! free and movealloc subroutines have exactly the same interface as in
! the case integer(ip) and therefore need not (cannot) be included here.
!***********************************************************************
#undef generic_memfree_interface
#undef generic_memmovealloc_interface
module mem_ip_igp_allocatable_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   integer(igp)
# define var_size   igp
# define bound_kind ip
# include "mem_header.i90"
  public :: memalloc,  memrealloc
contains
# include "mem_body.i90"
end module mem_ip_igp_allocatable_names
!***********************************************************************
!***********************************************************************
! Pointer routines
!***********************************************************************
!***********************************************************************
# define var_attr pointer
# define point(a,b) b => a
# define generic_status_test             associated
# define generic_memalloc_interface      memallocp
# define generic_memrealloc_interface    memreallocp
# define generic_memfree_interface       memfreep
!***********************************************************************
! integer(ip)
!***********************************************************************
module mem_ip_pointer_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type integer(ip)
# define var_size ip
# define bound_kind ip
# include "mem_header.i90"
  public :: memallocp,  memreallocp,  memfreep
contains
# include "mem_body.i90"
end module mem_ip_pointer_names
!***********************************************************************
! integer(igp)
!***********************************************************************
module mem_igp_pointer_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   integer(igp)
# define var_size   igp
# define bound_kind igp
# include "mem_header.i90"
  public :: memallocp,  memreallocp,  memfreep
contains
# include "mem_body.i90"
end module mem_igp_pointer_names
!***********************************************************************
! real(rp)
!***********************************************************************
module mem_rp_pointer_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   real(rp)
# define var_size   rp
# define bound_kind ip
# include "mem_header.i90"
  public :: memallocp,  memreallocp,  memfreep
contains
#ifdef exception
#define exception_real
#endif
# include "mem_body.i90"
#ifdef exception
#undef exception_real
#endif
end module mem_rp_pointer_names
!***********************************************************************
! mixed integer(ip)+integer(igp) specialization
! free and movealloc subroutines have exactly the same interface as in
! the case integer(ip) and therefore need not (cannot) be included here.
!***********************************************************************
#undef generic_memfree_interface
module mem_ip_igp_pointer_names
use types_names
use mem_base_names
#ifdef memcheck
use iso_c_binding
#endif
  implicit none
  private
# define var_type   integer(igp)
# define var_size   igp
# define bound_kind ip
# include "mem_header.i90"
  public :: memallocp,  memreallocp
contains
# include "mem_body.i90"
end module mem_ip_igp_pointer_names

module memor_names
use mem_base_names
use mem_ip_allocatable_names
use mem_ieep_allocatable_names
use mem_igp_allocatable_names
use mem_rp_allocatable_names
use mem_lg_allocatable_names
use mem_ip_igp_allocatable_names
use mem_ip_pointer_names
use mem_igp_pointer_names
use mem_rp_pointer_names
use mem_ip_igp_pointer_names
end module memor_names

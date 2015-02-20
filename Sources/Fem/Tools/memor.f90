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
module memor
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
  use types
  use hash_table_class ! Only used ifdef memcheck but this avoids passing the definition to configure
#ifdef memcheck
  use iso_c_binding
#endif
#ifdef exception
  use, intrinsic :: ieee_arithmetic
  use, intrinsic :: ieee_exceptions
#endif
  implicit none
  private

  integer(ip), save        ::     &
       lumem = 0,                 &              ! Logical unit to write memory evolution
       luerr = 0,                 &              ! Logical unit to write errors
       meall = 0,                 &              ! Number of memory allocations
       medea = 0                                 ! Number of memory deallocations

  integer(imp) , save      ::     &
       memax = 0,                 &              ! Maximum memory allocation in bytes
       mecur = 0                                 ! Current memory allocation in bytes

  real(rp), save :: nan ! An undefined value, initialized to nan if exception macro is defined

#ifdef memcheck
  type mem_data
     character(30) :: file=' '
     integer(ip)   :: line=0
  end type mem_data
# include "debug.i90"

#define hash_table hash_table_mem
#define hash_node  hash_node_mem
#define key_type   integer(igp)
#define val_type   type(mem_data)
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


  ! Allocatables
  interface memalloc
     module procedure memalloca_ip1, memalloca_ip2, memalloca_ip3, memalloca_ip4, &
          &           memalloca_igp1, memalloca_igp2, memalloca_igp3, memalloca_igp4, & 
          &           memalloca_ip_igp1, memalloca_ip_igp2, memalloca_ip_igp3, memalloca_ip_igp4, &
          &           memalloca_rp1, memalloca_rp2, memalloca_rp3, memalloca_rp4, &
          &           memalloca_lg1, memalloca_lg2, memalloca_lg3, memalloca_lg4
  end interface memalloc
  interface memrealloc
     module procedure memrealloca_ip1, memrealloca_ip2, memrealloca_ip3, memrealloca_ip4, &
          &           memrealloca_igp1, memrealloca_igp2, memrealloca_igp3, memrealloca_igp4, &
          &           memrealloca_ip_igp1, memrealloca_ip_igp2, memrealloca_ip_igp3, memrealloca_ip_igp4, &
          &           memrealloca_rp1, memrealloca_rp2, memrealloca_rp3, memalloca_rp4, &
          &           memrealloca_lg1, memrealloca_lg2, memrealloca_lg3, memalloca_lg4
  end interface memrealloc
  interface memfree
     module procedure memfreea_ip1, memfreea_ip2, memfreea_ip3, memfreea_ip4, &
          &           memfreea_igp1, memfreea_igp2, memfreea_igp3, memfreea_igp4, &
          ! AFM: this set of four subroutines has exactly the same interface as the ones
          !      in the previous line (i.e., memfreea_ipgp*), then it is a mistake to consider
          !      them in the same generic interface
          ! &           memfreea_ip_igp1, memfreea_ip_igp2, memfreea_ip_igp3, memfreea_ip_igp4, &
          &           memfreea_rp1, memfreea_rp2, memfreea_rp3, memfreea_rp4, &
          &           memfreea_lg1, memfreea_lg2, memfreea_lg3, memfreea_lg4
  end interface memfree
  interface memmovealloc
     module procedure memmovealloc_ip1, memmovealloc_ip2, memmovealloc_ip3, memmovealloc_ip4, &
          &           memmovealloc_igp1, memmovealloc_igp2, memmovealloc_igp3, memmovealloc_igp4, &
          ! AFM: this set of four subroutines has exactly the same interface as the ones
          !      in the previous line (i.e., memmovealloc_ipgp*), then it is a mistake to consider
          !      them in the same generic interface
          ! &         memmovealloc_ip_igp1, memmovealloc_ip_igp2, memmovealloc_ip_igp3, memmovealloc_ip_igp4, & 
          &           memmovealloc_rp1, memmovealloc_rp2, memmovealloc_rp3, memmovealloc_rp4, &
          &           memmovealloc_lg1, memmovealloc_lg2, memmovealloc_lg3, memmovealloc_lg4
  end interface memmovealloc

  ! Pointers
  interface memallocp
     module procedure memallocp_ip1, memallocp_ip2, memallocp_ip3, memallocp_ip4, &
          &           memallocp_igp1, memallocp_igp2, memallocp_igp3, memallocp_igp4, &
          &           memallocp_ip_igp1, memallocp_ip_igp2, memallocp_ip_igp3, memallocp_ip_igp4, &
          &           memallocp_rp1, memallocp_rp2, memallocp_rp3, memallocp_rp4
  end interface memallocp
  interface memreallocp
     module procedure memreallocp_ip1, memreallocp_ip2, memreallocp_ip3, memreallocp_ip4, &
          &           memreallocp_igp1, memreallocp_igp2, memreallocp_igp3, memreallocp_igp4, &
          &           memreallocp_ip_igp1, memreallocp_ip_igp2, memreallocp_ip_igp3, memreallocp_ip_igp4, &
          &           memreallocp_rp1, memreallocp_rp2, memreallocp_rp3, memallocp_rp4
  end interface memreallocp
  interface memfreep
     module procedure memfreep_ip1, memfreep_ip2, memfreep_ip3, memfreep_ip4, &
          &           memfreep_igp1, memfreep_igp2, memfreep_igp3, memfreep_igp4, &
          ! AFM: this set of four subroutines has exactly the same interface as the ones
          !      in the previous line (i.e., memfreep_ipgp*), then it is a mistake to consider
          !      them in the same generic interface
          ! &           memfreep_ip_igp1, memfreep_ip_igp2, memfreep_ip_igp3, memfreep_ip_igp4, &
          &           memfreep_rp1, memfreep_rp2, memfreep_rp3, memfreep_rp4
  end interface memfreep

  public :: memalloc,  memrealloc,  memfree, memmovealloc, &
       &    memallocp, memreallocp, memfreep, & 
       &    fempar_memmax, fempar_memcur, meminit, memstatus,  &
       &    memsum, memsub, mem_alloc_error, mem_dealloc_error, mem_status_error

#ifdef memcheck
  public ::  mem_db_put, mem_db_del
#endif

contains

  !***********************************************************************
  !***********************************************************************
  ! Specific purpose hash table specification (only if memcheck is defined)
  !***********************************************************************
  !***********************************************************************
#ifdef memcheck
#define hash_table hash_table_mem
#define hash_node  hash_node_mem
#define key_type   integer(igp)
#define key_size   igp
#define val_type   type(mem_data)
!#define print_key_val     write(*,'(i10,2x,a10,2x,i10)') list%key, list%val%file, list%val%line
#define print_key_val     write(*,*) list%key, trim(list%val%file), list%val%line
#define put_hash_node     put_hash_node_mem   
#define get_hash_node     get_hash_node_mem   
#define del_hash_node     del_hash_node_mem   
#define free_hash_node    free_hash_node_mem  
#define print_hash_node   print_hash_node_mem 
#define init_hash_table   init_hash_table_mem 
#define put_hash_table    put_hash_table_mem  
#define del_hash_table    del_hash_table_mem  
#define get_hash_table    get_hash_table_mem  
#define free_hash_table   free_hash_table_mem 
#define print_hash_table  print_hash_table_mem
#define status_hash_table  status_hash_table_mem
#include "hash_table_body.i90"
#endif

  !***********************************************************************
  !***********************************************************************
  ! integer(ip), allocatable routines
  !***********************************************************************
  !***********************************************************************
# define var_type integer(ip)
# define var_size ip
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloca_ip1    
# define generic_memalloc_2      memalloca_ip2    
# define generic_memalloc_3      memalloca_ip3    
# define generic_memalloc_4      memalloca_ip4    
# define generic_memrealloc_1    memrealloca_ip1  
# define generic_memrealloc_2    memrealloca_ip2  
# define generic_memrealloc_3    memrealloca_ip3  
# define generic_memrealloc_4    memrealloca_ip4  
# define generic_memfree_1       memfreea_ip1     
# define generic_memfree_2       memfreea_ip2     
# define generic_memfree_3       memfreea_ip3     
# define generic_memfree_4       memfreea_ip4     
# define generic_memmovealloc_1  memmovealloc_ip1
# define generic_memmovealloc_2  memmovealloc_ip2
# define generic_memmovealloc_3  memmovealloc_ip3
# define generic_memmovealloc_4  memmovealloc_ip4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! integer(igp), allocatable routines
  !***********************************************************************
  !***********************************************************************
# define var_type   integer(igp)
# define var_size   igp
# define var_attr allocatable, target
# define generic_status_test allocated
# define point(a,b) call move_alloc(a,b)
# define bound_kind igp

# define generic_memalloc_1      memalloca_igp1    
# define generic_memalloc_2      memalloca_igp2    
# define generic_memalloc_3      memalloca_igp3    
# define generic_memalloc_4      memalloca_igp4    
# define generic_memrealloc_1    memrealloca_igp1  
# define generic_memrealloc_2    memrealloca_igp2  
# define generic_memrealloc_3    memrealloca_igp3  
# define generic_memrealloc_4    memrealloca_igp4  
# define generic_memfree_1       memfreea_igp1     
# define generic_memfree_2       memfreea_igp2     
# define generic_memfree_3       memfreea_igp3     
# define generic_memfree_4       memfreea_igp4     
# define generic_memmovealloc_1  memmovealloc_igp1
# define generic_memmovealloc_2  memmovealloc_igp2
# define generic_memmovealloc_3  memmovealloc_igp3
# define generic_memmovealloc_4  memmovealloc_igp4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! mixed integer(ip)+integer(igp) specialization, allocatable routines
  !***********************************************************************
  !***********************************************************************
# define var_type integer(igp)
# define var_size igp
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_status_test     allocated
# define generic_memalloc_1      memalloca_ip_igp1    
# define generic_memalloc_2      memalloca_ip_igp2    
# define generic_memalloc_3      memalloca_ip_igp3    
# define generic_memalloc_4      memalloca_ip_igp4    
# define generic_memrealloc_1    memrealloca_ip_igp1  
# define generic_memrealloc_2    memrealloca_ip_igp2  
# define generic_memrealloc_3    memrealloca_ip_igp3  
# define generic_memrealloc_4    memrealloca_ip_igp4  
# define generic_memfree_1       memfreea_ip_igp1     
# define generic_memfree_2       memfreea_ip_igp2     
# define generic_memfree_3       memfreea_ip_igp3     
# define generic_memfree_4       memfreea_ip_igp4     
# define generic_memmovealloc_1  memmovealloc_ip_igp1
# define generic_memmovealloc_2  memmovealloc_ip_igp2
# define generic_memmovealloc_3  memmovealloc_ip_igp3
# define generic_memmovealloc_4  memmovealloc_ip_igp4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! real(rp), allocatable routines
  !***********************************************************************
  !***********************************************************************
# define var_type real(rp)
# define var_size rp
# define var_attr allocatable, target
# define generic_status_test allocated
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_memalloc_1      memalloca_rp1    
# define generic_memalloc_2      memalloca_rp2    
# define generic_memalloc_3      memalloca_rp3    
# define generic_memalloc_4      memalloca_rp4    
# define generic_memrealloc_1    memrealloca_rp1  
# define generic_memrealloc_2    memrealloca_rp2  
# define generic_memrealloc_3    memrealloca_rp3  
# define generic_memrealloc_4    memrealloca_rp4  
# define generic_memfree_1       memfreea_rp1     
# define generic_memfree_2       memfreea_rp2     
# define generic_memfree_3       memfreea_rp3     
# define generic_memfree_4       memfreea_rp4     
# define generic_memmovealloc_1  memmovealloc_rp1
# define generic_memmovealloc_2  memmovealloc_rp2
# define generic_memmovealloc_3  memmovealloc_rp3
# define generic_memmovealloc_4  memmovealloc_rp4
#ifdef exception
#define exception_real
#endif
# include "memor.i90"
#ifdef exception
#undef exception_real
#endif
  !***********************************************************************
  !***********************************************************************
  ! logical(lg), allocatable routines
  !***********************************************************************
  !***********************************************************************
# define var_type logical(lg)
# define var_size lg
# define var_attr allocatable, target
# define generic_status_test allocated
# define point(a,b) call move_alloc(a,b)
# define bound_kind ip

# define generic_memalloc_1      memalloca_lg1    
# define generic_memalloc_2      memalloca_lg2    
# define generic_memalloc_3      memalloca_lg3    
# define generic_memalloc_4      memalloca_lg4    
# define generic_memrealloc_1    memrealloca_lg1  
# define generic_memrealloc_2    memrealloca_lg2  
# define generic_memrealloc_3    memrealloca_lg3  
# define generic_memrealloc_4    memrealloca_lg4  
# define generic_memfree_1       memfreea_lg1     
# define generic_memfree_2       memfreea_lg2     
# define generic_memfree_3       memfreea_lg3     
# define generic_memfree_4       memfreea_lg4     
# define generic_memmovealloc_1  memmovealloc_lg1
# define generic_memmovealloc_2  memmovealloc_lg2
# define generic_memmovealloc_3  memmovealloc_lg3
# define generic_memmovealloc_4  memmovealloc_lg4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! integer(ip), pointer routines
  !***********************************************************************
  !***********************************************************************
# define var_type integer(ip)
# define var_size ip
# define var_attr pointer
# define point(a,b) b => a
# define bound_kind ip

# define generic_status_test     associated
# define generic_memalloc_1      memallocp_ip1    
# define generic_memalloc_2      memallocp_ip2    
# define generic_memalloc_3      memallocp_ip3    
# define generic_memalloc_4      memallocp_ip4    
# define generic_memrealloc_1    memreallocp_ip1  
# define generic_memrealloc_2    memreallocp_ip2  
# define generic_memrealloc_3    memreallocp_ip3  
# define generic_memrealloc_4    memreallocp_ip4  
# define generic_memfree_1       memfreep_ip1     
# define generic_memfree_2       memfreep_ip2     
# define generic_memfree_3       memfreep_ip3     
# define generic_memfree_4       memfreep_ip4     
# undef generic_memmovealloc_1
# undef generic_memmovealloc_2
# undef generic_memmovealloc_3
# undef generic_memmovealloc_4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! mixed integer(ip)+integer(igp) specialization, pointer routines
  !***********************************************************************
  !***********************************************************************
# define var_type integer(igp)
# define var_size ip
# define var_attr pointer
# define point(a,b) b => a
# define bound_kind ip

# define generic_status_test     associated
# define generic_memalloc_1      memallocp_ip_igp1    
# define generic_memalloc_2      memallocp_ip_igp2    
# define generic_memalloc_3      memallocp_ip_igp3    
# define generic_memalloc_4      memallocp_ip_igp4    
# define generic_memrealloc_1    memreallocp_ip_igp1  
# define generic_memrealloc_2    memreallocp_ip_igp2  
# define generic_memrealloc_3    memreallocp_ip_igp3  
# define generic_memrealloc_4    memreallocp_ip_igp4  
# define generic_memfree_1       memfreep_ip_igp1     
# define generic_memfree_2       memfreep_ip_igp2     
# define generic_memfree_3       memfreep_ip_igp3     
# define generic_memfree_4       memfreep_ip_igp4     
# undef generic_memmovealloc_1
# undef generic_memmovealloc_2
# undef generic_memmovealloc_3
# undef generic_memmovealloc_4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! integer(igp), pointer routines
  !***********************************************************************
  !***********************************************************************
# define var_type integer(igp)
# define var_size igp
# define var_attr pointer
# define point(a,b) b => a
# define bound_kind igp

# define generic_status_test     associated
# define generic_memalloc_1      memallocp_igp1    
# define generic_memalloc_2      memallocp_igp2    
# define generic_memalloc_3      memallocp_igp3    
# define generic_memalloc_4      memallocp_igp4    
# define generic_memrealloc_1    memreallocp_igp1  
# define generic_memrealloc_2    memreallocp_igp2  
# define generic_memrealloc_3    memreallocp_igp3  
# define generic_memrealloc_4    memreallocp_igp4  
# define generic_memfree_1       memfreep_igp1     
# define generic_memfree_2       memfreep_igp2     
# define generic_memfree_3       memfreep_igp3     
# define generic_memfree_4       memfreep_igp4     
# undef generic_memmovealloc_1
# undef generic_memmovealloc_2
# undef generic_memmovealloc_3
# undef generic_memmovealloc_4
# include "memor.i90"

  !***********************************************************************
  !***********************************************************************
  ! real(rp), pointer routines
  !***********************************************************************
  !***********************************************************************
# define var_type real(rp)
# define var_size rp
# define var_attr pointer
# define point(a,b) b => a
# define bound_kind ip

# define generic_status_test     associated
# define generic_memalloc_1      memallocp_rp1    
# define generic_memalloc_2      memallocp_rp2    
# define generic_memalloc_3      memallocp_rp3    
# define generic_memalloc_4      memallocp_rp4    
# define generic_memrealloc_1    memreallocp_rp1  
# define generic_memrealloc_2    memreallocp_rp2  
# define generic_memrealloc_3    memreallocp_rp3  
# define generic_memrealloc_4    memreallocp_rp4  
# define generic_memfree_1       memfreep_rp1     
# define generic_memfree_2       memfreep_rp2     
# define generic_memfree_3       memfreep_rp3     
# define generic_memfree_4       memfreep_rp4     
# undef generic_memmovealloc_1
# undef generic_memmovealloc_2
# undef generic_memmovealloc_3
# undef generic_memmovealloc_4

#ifdef exception
#define exception_real
#endif

# include "memor.i90"

#ifdef exception
#undef exception_real
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
    stop
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
    stop
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
    stop
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
    type(mem_data) :: val, old_val
    integer(ip)    :: istat
    key = transfer(varptr, 0_igp)

    if(present(file)) val%file = file
    if(present(line)) val%line = line
    !write(*,*) 'MEM_DB Storing ', key, varptr, file, line
    call mem_db%put(key, val, istat)
    if(istat/=now_stored) then
       write(*,*) 'An error ocurred storing allocation in mem_db: ', stat(istat)
       write(*,*) 'Called from '//file , line
       write(*,*) 'key: ', key, 'varptr: ', varptr ! DBG
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
    key = transfer(varptr, 0_igp)
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

end module memor

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
# include "debug.i90"
module fem_space_types
  use types
  use memor
  use sort_names
#ifdef memcheck
  use iso_c_binding
#endif

  implicit none
  private

  integer(ip), parameter       :: max_nnode  = 512  ! Maximum amount of nodes in an element (512=2**9)
  integer(ip), parameter       :: max_nobje  = 28   ! Maximum amount of objects in an element
  integer(ip), parameter       :: ht_length  = 10   ! hash tables length
  integer(ip), parameter       :: max_eltype = 12   ! Total number of types
  integer(ip), parameter       :: max_ndime  = 3    ! Maximum amount of space dimension
  integer(ip), parameter       :: max_elinf  = 8    ! Amount of different interpolations in a mesh
  integer(ip), parameter       :: max_order  = 7    ! Maximum interpolation order admited

  integer(ip), parameter       :: scond_off = 0
  integer(ip), parameter       :: scond_on  = 1

  ! IMPORTANT: All element identifiers MUST BE POSITIVE!!
  integer(ip), parameter       :: NULL_type_id = 0
  integer(ip), parameter       :: P_type_id    = 1
  integer(ip), parameter       :: Q_type_id    = 2
  integer(ip), parameter       :: max_FE_types = 2

  ! Types
  type fem_fixed_info_pointer
     type(fem_fixed_info)     , pointer :: p => NULL()  
  end type fem_fixed_info_pointer

  type fem_fixed_info

     integer(ip)              ::    &
          ftype,                    &        ! type of fem, e.g. 'P' 'Q' 'prism'...
          ndime,                    &        ! ndime
          order,                    &        ! FE order
          nobje,                    &        ! Number of objects
          nnode                              ! Number of nodes

     integer(ip) ::                 &
	  nobje_dim(5),             &        ! Pointer to object for each dimension
	  nodes_obj(4)                       ! Stores the number of nodes in each object
     ! SB: nodes_obj SHOULD NOT BE used because it cannot accomodate PRISMS, PYRAMIDS, etc.

     integer(ip), allocatable  :: o(:)       ! Orientation of the objects

     type(list)   :: ndxob       !array of interior nodes per object
     type(list)   :: ntxob       !array of all nodes per object
     type(list)   :: crxob       !array of corners per object
     type(list)   :: ndxob_int   !array of interior nodes per object when all nodes belong to interior
     type(list)   :: obxob       !array that lists all the objects in an object (idem ntxob for p = 2)

  end type fem_fixed_info

  ! Parameters 
  public :: max_nobje, ht_length, max_eltype, max_ndime, max_nnode, max_order
  public :: max_FE_types, max_elinf
  public :: scond_off, scond_on ! Static condensation flags (not active yet)
  public :: P_type_id, Q_type_id, NULL_type_id

  ! Types
  public :: fem_fixed_info, fem_fixed_info_pointer

  ! Functions
  public :: fem_element_fixed_info_create, fem_element_fixed_info_free, &
       fem_element_fixed_info_write, permute_nodes_object, get_order

  ! Functions
  public :: Q_ijkg, Q_gijk, Q_nnods, Q_set_integ
  public :: Q_coord_1d, Q_refcoord, Q_face_outno, Q_fixed_info_fill

  ! Functions
  public :: P_set_integ, P_refcoord

  !***********************************************************************
  ! Allocatable arrays of type(fem_fixed_info_pointer)
  !***********************************************************************
# define var_attr allocatable, target
# define point(a,b) call move_alloc(a,b)
# define generic_status_test             allocated
# define generic_memalloc_interface      memalloc
# define generic_memrealloc_interface    memrealloc
# define generic_memfree_interface       memfree
# define generic_memmovealloc_interface  memmovealloc

# define var_type type(fem_fixed_info_pointer)
# define var_size 8
# define bound_kind ip
# include "mem_header.i90"

  public :: memalloc,  memrealloc,  memfree, memmovealloc

contains

# include "mem_body.i90"

  !==================================================================================================
  subroutine fem_element_fixed_info_create ( f_info, f_type, f_order, dim_space,created)
    implicit none
    ! Parameters
    type(fem_fixed_info),  intent(inout) :: f_info 
    integer(ip)          ,  intent(in)    :: f_type, f_order, dim_space
    logical(lg)             ,  intent(inout) :: created

    created = .false.
    if (f_type == P_type_id) then
       call P_fixed_info_fill( f_info, dim_space,f_order)
       created = .true.
    elseif (f_type == Q_type_id) then
       call Q_fixed_info_fill( f_info, dim_space,f_order)
       created = .true.
    end if
  end subroutine fem_element_fixed_info_create

  !==================================================================================================
  subroutine fem_element_fixed_info_write ( f )
    implicit none
    ! Parameters
    type(fem_fixed_info),  intent(inout) :: f

    integer(ip) :: i

    write(*,*) 'ftype', f%ftype
    write(*,*) 'order', f%order
    write(*,*) 'nobje', f%nobje
    write(*,*) 'nnode', f%nnode
    write(*,*) 'nobje_dim', f%nobje_dim
    write(*,*) 'nodes_obj', f%nodes_obj

    write(*,*) 'ndxob'
    do i=1,f%nobje+1
       write(*,*) f%ndxob%l(f%ndxob%p(i):f%ndxob%p(i+1)-1)
    end do

    write(*,*) 'ntxob'
    do i=1,f%nobje+1
       write(*,*) f%ntxob%l(f%ntxob%p(i):f%ntxob%p(i+1)-1)
    end do
    
    write(*,*) 'obxob'
    do i=1,f%nobje+1
       write(*,*) f%obxob%l(f%obxob%p(i):f%obxob%p(i+1)-1)
    end do

    write(*,*) 'crxob'
    do i=1,f%nobje+1
       write(*,*) f%crxob%l(f%crxob%p(i):f%crxob%p(i+1)-1)
    end do
  end subroutine fem_element_fixed_info_write

  !==================================================================================================
  subroutine fem_element_fixed_info_free ( f_info)
    implicit none
    ! Parameters
    type(fem_fixed_info),  intent(inout) :: f_info 

    !Deallocate nobje_dim, nodes_obj
    !call memfree(f_info%nobje_dim,__FILE__,__LINE__)
    !call memfree(f_info%nodes_obj,__FILE__,__LINE__)
    call memfree(        f_info%o,__FILE__,__LINE__)   

    !Deallocate arrays
    call memfree(f_info%ndxob%p,__FILE__,__LINE__)   !Pointer to f_info%ndxob%l for each object
    call memfree(f_info%ndxob%l,__FILE__,__LINE__)   !Array of interior nodes of each object
    call memfree(f_info%ndxob_int%p,__FILE__,__LINE__)   !Pointer to f_info%ndxob_int%l for each object
    call memfree(f_info%ndxob_int%l,__FILE__,__LINE__)   !Array of nodes of each object when all interior
    call memfree(f_info%ntxob%p,__FILE__,__LINE__)   !Pointer to ntxob%l for each object
    call memfree(f_info%ntxob%l,__FILE__,__LINE__)   !Array of all nodes of each object
    call memfree(f_info%obxob%p,__FILE__,__LINE__)   !Pointer to obxob%l for each object
    call memfree(f_info%obxob%l,__FILE__,__LINE__)   !Array of all objects of each object
    call memfree(f_info%crxob%p,__FILE__,__LINE__)   !Pointer to crxob%l for each object
    call memfree(f_info%crxob%l,__FILE__,__LINE__)   !Array of corners for each object
  end subroutine fem_element_fixed_info_free

  !==================================================================================================
  ! This routine gives a permutation vector 'permu' that gives the relative position of nodes in
  ! object o2 of element e2 wrt the nodes in object o1 in element e1
  subroutine permute_nodes_object(e1,e2,permu,o1,o2,ln1,ln2,od,q,subface1,subface2)
    implicit none
    ! Parameters
    type(fem_fixed_info), intent(in)   :: e1, e2   ! Info of the elements
    integer(ip)         , intent(out)  :: permu(:) ! Permutation vector
    integer(ip)         , intent(in)   :: o1,o2    ! Local identifier of the object in each element
    integer(ip)         , intent(in)   :: ln1(e1%nobje), ln2(e2%nobje) ! lnods of each object
    integer(ip)         , intent(in)   :: od       ! Dimension of the object
    integer(ip)         , intent(in)   :: q  
    integer(ip), optional, intent(in)  :: subface1,subface2

    ! Local variables
    integer(ip) :: i,c1,r, o=0, r0, num_corners

    ! TODO: CHECK THE R0 implementation, it is probably worng


    if (present(subface1)) then
       if (subface1 == 0) then
          r0=0
       else
          r0 = subface1 -1
       end if
       assert (subface2 == 0) 
    else
       r0 = 0
    end if

    permu = 1
    !c1 = ln1(e1%crxob%l(e1%crxob%p(o1)))  ! Global identifier of the object of the first corner
    c1 = ln1(e1%crxob%l(e1%crxob%p(o1)+r0))  ! Global identifier of the object of the first corner
    r = 1
    do i = e2%crxob%p(o2),e2%crxob%p(o2+1)-1
       if ( ln2(e2%crxob%l(i)) == c1 ) exit
       r = r+1
    end do
    check ( ln2(e2%crxob%l(i)) == c1 )

    if (r0>0) then
       r = r-r0
       if (r < 1) then
          num_corners = e2%crxob%p(o2+1)- e2%crxob%p(o2)
          r = r + num_corners 
       end if
    end if

    if (od == 2) then
       o = modulo(e1%o(o1)+e1%o(o2)+1,2)
    else
       o = 0
    end if

    if (e2%ftype == P_type_id) then
       call P_permute_or( permu,q,o,r,od )
    elseif (e2%ftype == Q_type_id) then
       call Q_permute_or( permu,q,o,r,od )
    else
       write(*,*) __FILE__,__LINE__,'WARNING! elem type not identified!'
    end if


  end subroutine permute_nodes_object

  !==================================================================================================
  !==================================================================================================
  !
  ! P_fem_fixed_info
  !
  !==================================================================================================
  !==================================================================================================
  subroutine  P_fixed_info_fill(fefi,nd,p)
    implicit none
    ! Parameters
    type(fem_fixed_info) ,  intent(inout) :: fefi
    integer(ip)          ,  intent(in)    :: nd, p ! nd = dimension. p = order.

    ! Local variables
    integer(ip)               :: i,j,k
    integer(ip)               :: co,c2,c3,c4          ! Counters
    integer(ip)               :: aux1,aux2,aux3
    integer(ip)               :: idm(nd+1)         ! Stores the corners f each element
    integer(ip)               :: nodes(nd,nd+1)    ! Coordinates of the corners of the element
    integer(ip)               :: objec(nd,nd+1)    ! Coordinates that define an object
    integer(ip)               :: no                ! #objects in the element
    integer(ip)               :: nn                ! #nodes 
    integer(ip)               :: nt                ! Length of ntxob
    integer(ip)               :: nc                ! #corners x #{objects delimited by each corner}
    integer(ip), allocatable  :: idcro(:,:)        ! Array of corners belonging to the nth object

    nt = 0
    nc = 0
    nn = 0

    ! Initialize nobje_dim, nodes_obj
    ! call memalloc(nd+2,fefi%nobje_dim,__FILE__,__LINE__)
    ! call memalloc(nd+1,fefi%nodes_obj,__FILE__,__LINE__)
    fefi%nodes_obj = 0
    fefi%nobje_dim = 0
    fefi%nobje_dim(1) = 1

    ! Fill nobje_dim, nodes_obj and compute nt, nc, nn
    do k = 0,nd
       i = bnm(nd+1,k+1)          ! #objects of dimension k
       j = P_nnods(k,p)           ! #nodes object dim k order p
       nt = nt + i*j              ! nodes in the clousure of the object
       nc = nc + i*(k+1)          ! corners delimiting objects of dimension k
       ! Pointer to obj id by dim. Local obj of dim k are fefi%nobje_dim(k):fefi%nobje_dim(k+1)
       fefi%nobje_dim(k+2) = fefi%nobje_dim(k+1) + i  
       ! #nodes in objects of dimension k
       fefi%nodes_obj(k+1) = inods(k,p)
       nn = nn + i*fefi%nodes_obj(k+1)      
    end do

    ! Set no
    no = fefi%nobje_dim(nd+2)-1    ! #objects

    ! Set constant values of fem_fixed_info
    fefi%ftype = P_type_id
    fefi%order = p
    fefi%nobje = no-1
    fefi%nnode = nn

    ! Allocate arrays
    call memalloc(no,        fefi%o,__FILE__,__LINE__)   ! Array of orientation of each object
    call memalloc(no+1,fefi%ndxob%p,__FILE__,__LINE__)   !Pointer to fefi%ndxob%l for each object
    call memalloc(nn,  fefi%ndxob%l,__FILE__,__LINE__)   !Array of interior nodes of each object
    call memalloc(no+1,fefi%ndxob_int%p,__FILE__,__LINE__)   !Pointer to fefi%ndxob%l for each object
    call memalloc(nn,  fefi%ndxob_int%l,__FILE__,__LINE__)   !Array of interior nodes of each object
    call memalloc(no+1,fefi%ntxob%p,__FILE__,__LINE__)   !Pointer to ntxob%l for each object
    call memalloc(nt,  fefi%ntxob%l,__FILE__,__LINE__)   !Array of all nodes of each object
    call memalloc(no+1,fefi%obxob%p,__FILE__,__LINE__)   !Pointer to obxob%l for each object
    call memalloc(nt,  fefi%obxob%l,__FILE__,__LINE__)   !Array of all objects of each object
    call memalloc(no+1,fefi%crxob%p,__FILE__,__LINE__)   !Pointer to crxob%l for each object
    call memalloc(nc,  fefi%crxob%l,__FILE__,__LINE__)   !Array of corners for each object
    call memalloc(nd+2,no,idcro,__FILE__,__LINE__) !Array of dim and corners belonging to each object

    fefi%ndxob%p=0   !Pointer to fefi%ndxob%l for each object
    fefi%ndxob%l=0     !Array of interior nodes of each object
    fefi%ndxob_int%p=0   !Pointer to fefi%ndxob_int%l for each object
    fefi%ndxob_int%l=0     !Array of interior nodes of each object when all nodes belong to volume
    fefi%ntxob%p=0   !Pointer to ntxob%l for each object
    fefi%ntxob%l=0     !Array of all nodes of each object
    fefi%obxob%p=0   !Pointer to ntxob%l for each object
    fefi%obxob%l=0     !Array of all objects of each object
    fefi%crxob%p=0   !Pointer to crxob%l for each object
    fefi%crxob%l=0     !Array of corners for each object

    ! Create auxiliar matrix nodes with the coordinates of the corners
    nodes = 0
    do k=1,nd
       nodes(k,k+1) = 1
    end do

    ! Initialize pointers
    fefi%ndxob%p(1) = 1
    fefi%ntxob%p(1) = 1
    fefi%crxob%p(1) = 1

    ! Loop over dimensions
    do k = 0,nd
       aux1 = inods(k,p)            ! interior nodes for an object of dim k
       aux3 = P_nnods(k,p)            ! Total nodes for an object of dim k
       aux2 = k+1                   ! Corners for an object of dim k

       !Loop over objects of dimension k
       do i = fefi%nobje_dim(k+1),fefi%nobje_dim(k+2)-1 
          fefi%ndxob%p(i+1) = fefi%ndxob%p(i) + aux1 !assign pointers
          fefi%ntxob%p(i+1) = fefi%ntxob%p(i) + aux3 !assign pointers
          fefi%crxob%p(i+1) = fefi%crxob%p(i) + aux2 !assign pointers 
       end do
    end do

    fefi%ndxob_int%p = 1
    fefi%ndxob_int%p(no+1) = nn+1

    ! Each object of dimension k is defined by a set of k+1 corners (idcro stores this info)
    i = 1
    j = 1
    k = 1 
    idcro = -1
    do k=0,nd
       call crxob(idcro,k,i,idm,nd,1,no,j)
    end do

    ! Initialize counters
    co = 0 ! Counter of object
    c2 = 1 ! ndxob%p counter
    c3 = 0 ! crxob%p counter
    c4 = 1 ! ntxob%p counter

    ! Loop over objects dimensions
    do k = 0,nd

       ! Loop over the number of objects of dimension k
       do i=1, bnm(nd+1,k+1)
          ! Fill crxob%l for object co
          co = co+1 
          call P_orientation_object(fefi%o(co),k,nd,i)
          fefi%crxob%l(c3+1:c3+k+1) = idcro(2:k+2,co)
          c3 = c3 + k +1

          ! Objec stores the coordinates of the corners defining object co
          objec(:,1) = nodes(:,idcro(2,co))
          do j=1,k
             objec(:,j+1) = nodes(:,idcro(j+2,co)) - nodes(:,idcro(2,co))
          end do

          ! Fill ntxobj and ndxobj for object co
          call ntxob_fill(fefi%ntxob%l,c4,0,p,1,idm,nd,k,p,objec,nt)
          call ntxob_fill(fefi%ndxob%l,c2,1,p-1,1,idm,nd,k,p,objec,nt)
       end do
    end do

    do i=1,nn
       fefi%ndxob_int%l(i) = i
    end do

    ! Deallocation of variable
    call memfree(idcro,__FILE__,__LINE__)

    ! write(*,*) 'orientation objects'
    ! do k = 1,nd
    !    write(*,*) 'dime', k, '--------------------------'
    !    write(*,*) fefi%o(fefi%nobje_dim(k):fefi%nobje_dim(k+1)-1)
    ! end do
    ! write(*,*) 'no+1', no+1, 'ndxob%p'
    ! do k = 1,no+1
    !    write(*,*) fefi%ndxob%p(k), ', &'
    ! end do
    ! write(*,*) 'nn', nn, 'ndxob%l'
    ! do k = 1,nn
    !    write(*,*) fefi%ndxob%l(k), ', &'
    ! end do

    ! write(*,*) 'no+1', no+1, 'ntxob%p'
    ! do k = 1,no+1
    !    write(*,*) fefi%ntxob%p(k), ', &'
    ! end do
    ! write(*,*) 'nt', nt, 'ntxob%l'
    ! do k = 1,nt
    !    write(*,*) fefi%ntxob%l(k), ', &'
    ! end do

    ! write(*,*) 'no+1', no+1, 'crxob%p'
    ! do k = 1,no+1
    !    write(*,*) fefi%crxob%p(k), ', &'
    ! end do
    ! write(*,*) 'nc', nc, 'crxob%l'
    ! do k = 1,nc
    !    write(*,*) fefi%crxob%l(k), ', &'
    ! end do
  end subroutine P_fixed_info_fill

  !==================================================================================================
  subroutine P_orientation_object(o,od,nd,io)
    implicit none
    ! Parameters
    integer(ip), intent(out) :: o
    integer(ip), intent(in)  :: od,nd,io  ! io=numbering of the object in the od dimension

    if (nd == 3 .and. od == 2) then
       o = modulo(io+1,2)
    elseif (nd>3) then
       write(*,*) __FILE__,__LINE__,'WARNING!! the orientation is not defined for dimension >3'
    else
       o = 0
    end if
  end subroutine P_orientation_object

  !==================================================================================================
  subroutine P_o2n_2d_create(o2n,p)
    implicit none
    ! Parameters
    integer(ip), allocatable, intent(inout) :: o2n(:,:)
    integer(ip)             , intent(in)    :: p

    integer(ip) :: o,r,o_r

    call memalloc(6,int(((p+1)**2+p+1)/2),o2n,__FILE__,__LINE__)

    ! Loop over possible orientations of a face
    do o = 0,1
       ! Loop over possible rotations of a face
       do r= 1,3  
          o_r = 3*o+r
          call P_permute_or_2d(o2n(:,o_r),p,o,r)
       end do
    end do
  end subroutine P_o2n_2d_create

  !=================================================================================================
  ! This subroutine gives the reodering (o2n) of the nodes of an object given an orientation 'o'
  ! and a delay 'r' wrt to a refence element sharing the same object.
  subroutine  P_permute_or( o2n,p,o,r,nd )
    implicit none
    integer(ip), intent(in)    :: p,o,r,nd
    integer(ip), intent(inout) :: o2n(:)

    if     (nd == 0) then
       o2n(1) = 1
    elseif (nd == 1) then
       call P_permute_or_1d(o2n(1:p+1),p,r)
    elseif (nd == 2) then
       call P_permute_or_2d(o2n(1:int(((p+1)**2+p+1)/2)),p,o,r)
    else
       o2n(1) = 0
       o2n(2) = 1/o2n(1)
       write(*,*) __FILE__,__LINE__,'P_permute_or:: WARNING! Permutations not given for nd>3'
    end if
  end subroutine P_permute_or

  !=================================================================================================
  subroutine  P_permute_or_1d( o2n,p,r )
    implicit none
    integer(ip), intent(in)    :: p,r
    integer(ip), intent(inout) :: o2n(p+1)

    ! Local variables
    integer(ip) :: i

    ! Generic loop+rotation identifier  
    if (r==1) then
       o2n = (/(i,i=1,p+1)/)
    elseif (r==2) then
       o2n = (/(p+1-i,i=0,p)/)
    else
       write(*,*) __FILE__,__LINE__,'P_permute_or_1d:: ERROR! Delay cannot be >1 for edges'
    end if
  end subroutine P_permute_or_1d

  !=================================================================================================
  subroutine  P_permute_or_2d( o2n,p,o,r )
    implicit none
    integer(ip), intent(in)    :: p,o,r
    integer(ip), intent(inout) :: o2n(int(((p+1)**2+p+1)/2))

    ! Local variables
    integer(ip) :: o_r,i,j,ij_t(3)     ! ij_t = (i,j,p-i-j)
    integer(ip) :: ij_n(2),go,gn
    integer(ip) :: faces_perm_tet(4,4) = reshape((/ 1, 0, 1, 0, &
         &                                         0, 1, 0, 1, &
         &                                         1, 0, 1, 0, &
         &                                         0, 1, 0, 1/), (/4,4/) )
    integer(ip) :: ij_perm_tet(2,6) = reshape((/ 1, 2, 2, 3, 3, 1, 2, 1, 3, 2, 1, 3/), (/2,6/) )

    ! Generic loop+rotation identifier  
    o_r = 3*o+r
    do j = 0,p
       ij_t(2) = j
       do i = 0,p-j
          ij_t(1) = i
          ij_t(3) = p-i-j
          ! Get the global numbering of node (i,j)
          go  = gijk(ij_t(1:2),2,p)
          ! i,j coordinates for the o_r permutation
          ij_n(1:2) = ij_t(ij_perm_tet(1:2,o_r)) 
          ! Store the global numbering of node ij_n 
          o2n(go)    = gijk(ij_n,2,p)
       end do
    end do
  end subroutine P_permute_or_2d

  ! !=================================================================================================
  ! ! FC(k)=k! computes the factorial of k 
  ! integer(ip) function fc(i)
  !   implicit none
  !   integer(ip) :: i, k
  !   fc = 1
  !   do k=2,i
  !      fc = fc*k
  !   end do
  ! end function fc

  ! !=================================================================================================
  ! ! BNM(A,B)=A!/((A-B)!B!) computes the binomial coefficient of (A,B), A>B
  ! integer (ip) function bnm(a,b)
  !   implicit none
  !   integer(ip) :: a,b
  !   if (a >= b) then
  !      bnm = int(fc(a)/(fc(b)*fc(a-b)))
  !   else
  !      write(*,*) 'ERROR: no binomial coef for b > a'
  !      check(.false.)
  !   end if
  ! end function bnm

  !=================================================================================================
  ! P_NNODS(k,p) computes the #nodes in a simplex of dim k and order p
  integer(ip) recursive function P_nnods(k,p) RESULT(nnods)
    implicit none
    !integer(ip)                :: nnods
    integer(ip), intent(in)    :: k, p

    integer(ip) :: q
    if (k == 0) then
       nnods = 1
    elseif (p == 0) then
       nnods = 1
    elseif (k == 1) then
       nnods = p+1
    elseif (k == 2) then
       nnods = int((p+1)*(p+2)/2)
    else
       nnods = P_nnods(k-1,0)
       do q=1,p
          nnods = nnods + P_nnods(k-1,q)
       end do
    end if
  end function P_nnods

  !=================================================================================================
  ! PNODS(k,pi,pf) computes the sum of nodes in simplices of dimension k of order pi,...,pf
  integer(ip) function pnods(k,pi,pf)

    implicit none
    integer(ip), intent(in)    :: k, pi,pf

    integer(ip) :: q
    if (k == 0) then
       pnods = 1
    elseif (pi > pf) then
       pnods = 0
    else
       pnods = P_nnods(k-1,pi)
       do q=pi+1,pf
          pnods = pnods + P_nnods(k-1,q)
       end do
    end if
  end function pnods

  !=================================================================================================
  ! INODS(k,p) computes the #{interior nodes} in a simplex of dim k and order p
  integer(ip) recursive function inods(k,p) RESULT(inodes)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: k, p

    ! Variables
    integer(ip) :: q
    if (k == 0) then
       inodes = 1
    elseif (p == 1) then
       inodes = 0
    elseif (k == 1) then
       inodes = p-1
    elseif (k == 2) then
       inodes = int((p-2)*(p-1)/2)
    else
       inodes = inods(k-1,1)
       do q=2,p-1
          inodes = inodes+ inods(k-1,q)
       end do
    end if
  end function inods

  !==================================================================================================
  ! CRXOB constructs idcro and idm
  ! idcro(:,cc) = (k,idm(1:k+1))
  ! k   :: object dimension
  ! i   :: position of idm we are modifying
  ! idm :: array of dimension k+1 that will give the corners in the element
  ! n   :: global dimension
  ! in  :: value of the corner we can count from (knowing that idm(i)<idm(i+1))
  ! no  :: total amount of objects
  ! cc  :: counter of the objects
  recursive subroutine crxob(idcro,k,i,idm,n,in,no,cc)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: k,i,n,in,no
    integer(ip), intent(inout) :: idcro(n+2,no),idm(n+1),cc

    ! Local variables
    integer(ip)  :: j,m

    do j=in,n+1
       idm(i) = j
       if(i+1<k+2) then
          call crxob(idcro,k,i+1,idm,n,idm(i)+1,no,cc)
       else
          idcro(1,cc) = k
          idcro(2:k+2,cc) = idm(1:k+1)
          cc = cc + 1
       end if
    end do

  end subroutine crxob

  !==================================================================================================
  ! NTXOB_FILL constructs ntxob%l
  recursive subroutine ntxob_fill(ntxob,c3,ini,end,i,idm,nd,k,p,objec,nt)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: ini,end,nd,k,p,objec(nd,nd+1),nt,i
    integer(ip), intent(inout) :: ntxob(nt),c3,idm(nd)

    ! Local variables
    integer(ip) :: j,ijk(nd),m

    if (k==0) then
       ! Dimension 0: objec gives the coordinates of the corner
       ijk = p*objec(:,1)
       ntxob(c3) = gijk(ijk,nd,p)
       c3 = c3 + 1
    else
       ! Loop over the possible values of the i-th factor
       do j=ini,end
          ! Assign value of the i-th factor
          idm(i) = j
          if (i<k) then
             ! Assign values of the i+1-th factor
             call ntxob_fill(ntxob,c3,ini,end-j,i+1,idm,nd,k,p,objec,nt)
          else
             ! Compute the coordinate of the c3-th node: objec(:,1)+idm*objec(:,2:k)
             ijk = p*objec(:,1)
             do m=1,k
                ijk = ijk+idm(k-m+1)*objec(:,m+1)
             end do
             ! Store in ntxob the corresponding identifier of the node
             ntxob(c3) = gijk(ijk,nd,p)
             c3 = c3 + 1
          end if
       end do
    end if
  end subroutine ntxob_fill

  !==================================================================================================
  ! GIJK(i,nd,p) returns the generic identifier of a node with coordinates i in an elem(nd,p)
  integer(ip) function gijk(i,nd,p)
    implicit none
    integer(ip) :: nd, i(nd),p,k,q

    if (sum(i)>p) then
       write(*,*) 'ijk', i
       write(*,*) 'tets_setting:: gijk:: ERROR i+j+k<=p'
       stop
    end if

    gijk = 1
    do k =1,nd-1
       q = p-i(k+1)
       gijk = gijk + pnods(k,q-i(k)+1,q)
    end do
    gijk = gijk + pnods(k,p-i(nd)+1,p)
  end function gijk

  !==================================================================================================
  ! P_IJKG(i,g,nd,p) returns coordinates of the g-th node in an elem(nd,p)
  subroutine P_ijkg(i,g,nd,p)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: nd,g,p
    integer(ip), intent(out) :: i(nd)

    ! Local variables
    integer(ip)              :: k,g0,g1,g2,j,q

    g0 = g-1
    j  = 0
    do while (pnods(nd,p-j+1,p)<=g0)
       j=j+1
    end do
    i(nd) = j-1
    g0 = g0 - pnods(nd,p-i(nd)+1,p)

    do k=1,nd-1
       g2 = 0
       j  = 0
       q = p-i(nd-k+1)
       do while (pnods(nd-k,q-j+1,q)<=g0)
          j=j+1
       end do
       i(nd-k) = j-1
       g0 = g0 - pnods(nd-k,q-i(nd-k)+1,q)
    end do
  end subroutine P_ijkg

  !==================================================================================================
  subroutine P_set_integ(ndime,order,nnode,ngaus,nlocs,lrule,llapl,mnode)
    implicit none
    ! Parameters
    integer(ip),           intent(in)  :: ndime,order
    integer(ip),           intent(out) :: nnode,ngaus,lrule,llapl,nlocs
    integer(ip), optional, intent(in)  :: mnode

    ! Assign nnode and ngaus
    nnode = P_nnods(ndime,order)
    if (present(mnode)) then
       ngaus = mnode
       if (ndime == 1) then
          nlocs = mnode
       elseif (ndime == 2) then
          nlocs = int(-0.5_rp+0.5*sqrt(1.0_rp+8*real(mnode)))
       else          
          write(*,*) __FILE__,__LINE__, 'WARNING! nlocs not computed'
       end if
       if(ndime==2) then
          if(mnode == 10) then ! P3
             ngaus = 13
             write(*,*) 'new_integration:: WARNING! nnode = 10, ngaus = 13' 
          end if
       elseif(ndime==3) then
          if( mnode == 10) then ! P2
             ngaus = 11
             write(*,*) 'new_integration:: WARNING! nnode = 10, ngaus = 11'
          else if( mnode > 13 ) then ! P3
             ngaus = 14
             write(*,*) 'new_integration:: ERROR! quadrature for 3D tetrahedra up to 14 gauss points' 
             stop
          end if
       end if
    else
       ngaus = nnode
       nlocs = order + 1
       ! Assign some exceptions of ngaus
       if(ndime==2) then
          if(order == 3) then ! P3
             ngaus = 13
             write(*,*) 'new_integration:: WARNING! nnode = 10, ngaus = 13' 
          end if
       elseif(ndime==3) then
          if( order == 2) then ! P2
             ngaus = 11
             write(*,*) 'new_integration:: WARNING! nnode = 10, ngaus = 11'
          else if( order == 3 ) then ! P3
             ngaus = 14
             write(*,*) 'new_integration:: ERROR! quadrature for 3D tetrahedra up to 14 gauss points' 
             stop
          end if
       end if
    end if

    ! Assign lrule
    if (ndime == 1) then
       lrule = ruqope_id
    else
       lrule = rutope_id
    end if

    ! Assign llapl
    if (order > 1 ) then
       llapl = 1
    else
       llapl = 0
    end if

  end subroutine P_set_integ

  !==================================================================================================
  ! Coord(:,g) = coordinates of the g-th node in the reference element
  subroutine P_refcoord (coord,nd,p,nn)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nd,p,nn
    real(rp),    intent(inout) :: coord(nd,nn)

    ! Local variables
    integer(ip) :: ijk(nd),i,d
    real(rp)    :: c1d(p+1)

    call P_coord_1d(c1d,p+1)

    do i=1,nn

       call P_ijkg(ijk,i,nd,p)
       do d=1,nd
          coord(d,i) = c1d(ijk(d)+1)
       end do
    end do
  end subroutine P_refcoord

  ! =================================================================================================
  ! Set n equidistant points in the segment [0,1]
  subroutine P_coord_1d (x,n)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: n
    real(rp)   , intent(out) :: x(n)

    ! Local variables
    integer(ip)              :: i

    do i = 0,n-1
       x(i+1) = real(i)/(real(n)-1)
    end do
  end subroutine P_coord_1d

  !==================================================================================================
  !==================================================================================================
  !
  ! Q_fem_fixed_info
  !
  !==================================================================================================
  !==================================================================================================
  !==================================================================================================
  subroutine  Q_fixed_info_fill(fefi,nd,p)
    implicit none
    ! Parameters
    type(fem_fixed_info), intent(inout) :: fefi
    integer(ip)          , intent(in)    :: nd,p ! nd = dimension. p = order.

    ! Local variables
    integer(ip)               :: i,j,k,l,m
    integer(ip)               :: aux1,aux2,aux3,aux4,aux(nd)
    integer(ip)               :: idm(nd),fdm(nd)
    integer(ip)               :: od,cd,kk,c
    integer(ip)               :: ijk(nd),ijk_g(nd)
    integer(ip)               :: c2,c3,c4,c5,c6,co
    integer(ip)               :: no  ! #objects in nd-quad
    integer(ip)               :: nn  ! #nodes in the nd-quad of order p
    integer(ip)               :: nt  ! sum of the #nodes (all, not only interior) of the objects
    integer(ip)               :: nt2 ! sum of the #objects (all, not only interior) of the objects
    integer(ip)               :: nc  ! #corners x #{objects delimited by each corner}
    integer(ip)               :: nod ! #corners in the quad. nod = 2^(nd)
    integer(ip), allocatable  :: auxt1(:,:),auxt2(:,:),auxt3(:,:),auxt4(:,:),auxt5(:,:), auxt6(:,:)
    integer(ip), allocatable  :: obdla(:,:),node2ob(:),ob2node(:)

    ! Initilize values
    no = 0
    nn = 0
    nt = 0
    nc = 0
    nt2 = 0
    nod = 0

    ! Initialize nobje_dim, nodes_obj
    !call memalloc(nd+2,fefi%nobje_dim,__FILE__,__LINE__)
    !call memalloc(nd+1,fefi%nodes_obj,__FILE__,__LINE__)
    fefi%nodes_obj = 0
    fefi%nobje_dim = 0
    fefi%nobje_dim(1) = 1

    do k = 0,nd
       i = int(2**(nd-k)*bnm(nd,k)) ! #objects of dimension k
       no = no + i                  ! compute #objects
       nn = nn + int(i*((p-1)**k))  ! #nodes inside objects of dimension k
       nt = nt + int(i*((p+1)**k))  ! nodes in the clousure of the object
       nt2 = nt2 + int(i*((3)**k))! objects in the clousure of the object
       nc = nc + int(i*2**k)        ! corners delimiting objects of dimension k
       nod = nod + bnm(nd,k)        ! #nodes/{displacement} (2^n = sum(bnm(n,k)), k=0,..,n)
       ! Pointer to obj id by dim. Local obj of dim k are fefi%nobje_dim(k):fefi%nobje_dim(k+1)
       fefi%nobje_dim(k+2) = fefi%nobje_dim(k+1) + i 
       ! #nodes in objects of dimension k 
       fefi%nodes_obj(k+1) = int((p-1)**k)
    end do

    ! Set constant values of fem_fixed_info
    fefi%ftype = Q_type_id
    fefi%order = p
    fefi%nobje = no-1
    fefi%nnode = int((p+1)**nd) 

    ! Allocate arrays
    call memalloc(no,        fefi%o,__FILE__,__LINE__)  ! Array of orientation of each object
    call memalloc(no+1,fefi%ndxob%p,__FILE__,__LINE__)  !Pointer to fefi%ndxob%l for each object
    call memalloc(nn,  fefi%ndxob%l,__FILE__,__LINE__)  !Array of interior nodes of each object
    call memalloc(no+1,fefi%ndxob_int%p,__FILE__,__LINE__)  !Pointer to fefi%ndxob%l for each object
    call memalloc(nn,  fefi%ndxob_int%l,__FILE__,__LINE__)  !Array of interior nodes of each object
    call memalloc(no+1,fefi%ntxob%p,__FILE__,__LINE__)  !Pointer to ntxob%l for each object
    call memalloc(nt,  fefi%ntxob%l,__FILE__,__LINE__)  !Array of all nodes of each object
    call memalloc(no+1,fefi%obxob%p,__FILE__,__LINE__)  !Pointer to obxob%l for each object
    call memalloc(nt2,  fefi%obxob%l,__FILE__,__LINE__)  !Array of all objects of each object
    call memalloc(no+1,fefi%crxob%p,__FILE__,__LINE__)  !Pointer to crxob%l for each object
    call memalloc(nc,  fefi%crxob%l,__FILE__,__LINE__)  !Array of corners for each object
    call memalloc(nod,nd+1,obdla,__FILE__,__LINE__)
    call memalloc(no, node2ob,__FILE__,__LINE__)        ! Auxiliar array
    call memalloc(no, ob2node,__FILE__,__LINE__)        ! Auxiliar array

    fefi%ndxob%p=0   !Pointer to fefi%ndxob%l for each object
    fefi%ndxob%l=0   !Array of interior nodes of each object
    fefi%ndxob_int%p=0   !Pointer to fefi%ndxob%l for each object
    fefi%ndxob_int%l=0   !Array of interior nodes of each object
    fefi%ntxob%p=0   !Pointer to ntxob%l for each object
    fefi%obxob%l=0   !Array of all nodes of each object
    fefi%obxob%p=0   !Pointer to obxob%l for each object
    fefi%ntxob%l=0   !Array of all nodes of each object
    fefi%crxob%p=0   !Pointer to crxob%l for each object
    fefi%crxob%l=0   !Array of corners for each object


    !Initialize pointers
    fefi%ndxob%p(1) = 1
    fefi%ntxob%p(1) = 1
    fefi%obxob%p(1) = 1
    fefi%crxob%p(1) = 1

    !Loop over dimensions
    do k = 0,nd
       aux1 = int(((p-1)**k)) ! interior nodes for an object of dim k
       aux3 = int(((p+1)**k)) ! Total nodes for an object of dim k
       aux2 = int(2**k)       ! Corners for an object of dim k
       aux4 = int((3**k)) ! Corners for an object of dim k (idem p=2)

       ! Loop over objects of dimension k
       do i = fefi%nobje_dim(k+1),fefi%nobje_dim(k+2)-1 
          fefi%ndxob%p(i+1) = fefi%ndxob%p(i) + aux1 ! assign pointers
          fefi%ntxob%p(i+1) = fefi%ntxob%p(i) + aux3 ! assign pointers
          fefi%crxob%p(i+1) = fefi%crxob%p(i) + aux2 ! assign pointers
          fefi%obxob%p(i+1) = fefi%obxob%p(i) + aux4 ! assign pointers 
       end do
    end do

    fefi%ndxob_int%p = 1
    fefi%ndxob_int%p(no+1) = nn+1

    ! Initialize auxiliar values
    k = 0
    i = 0
    idm = 0
    j = 2

    ! Construction of obdla matrix
    ! For each object, up to a displacement, we have an identifier id={1..nod}
    ! obdla(id,1) = dimension of the object
    ! obdla(id,2:obdla(id,1)+1) = gives the directions that define the object
    obdla = -1
    obdla(1,1) = 0
    do od = 0,nd
       if (od > 0) then
          call r_dim(j,idm(1:od),k,i,nd,od,obdla,nod)
       end if
    end do

    ! Initialize auxiliar values
    idm = 0
    fdm = 0
    cd  = 0
    c2  = 0 ! fefi%ndxob%p counter
    c3  = 0 ! crxob%p counter
    c4  = 0 ! ntxob%p counter
    c5  = 0 ! obxob%p counter
    c6  = 0 ! ob2node   counter
    co  = 0 ! counter of objects

    ! Loop over objects dimensions
    do od = 0,nd

       ! Create ijk tables (od)
       ! Compute auxt1 the local numbering of the corners of an object of dimension nd-od
       ! It allows to know how many translations for each paralel set of objects
       if (od < nd) then
          call memalloc(nd-od,2**(nd-od),auxt1,__FILE__,__LINE__)
          auxt1 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,nd-od,auxt1,0,1)
       end if
       ! Compute auxt2 the local numbering of the corners in an object of dim od
       if (od >0) then
          call memalloc(od,2**(od),auxt2,__FILE__,__LINE__)
          auxt2 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,od,auxt2,0,1)

          if (p > 1) then
             ! Compute auxt3 the local numbering of the interior nodes in an object of dim od
             call memalloc(od,(p-1)**(od),auxt3,__FILE__,__LINE__)
             auxt3 = 0
             kk    = 0
             aux1  = 1
             ijk   = 0
             call Q_r_ijk(kk,aux1,ijk,od,auxt3,1,p-1)
          end if

          ! Compute auxt4 the local numbering of all nodes in an object of dim od
          call memalloc(od,(p+1)**(od),auxt4,__FILE__,__LINE__)
          auxt4 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,od,auxt4,0,p)

          ! Compute auxt5 the local numbering of all nodes in an object of dim od
          call memalloc(od,(2+1)**(od),auxt5,__FILE__,__LINE__)
          auxt5 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,od,auxt5,0,2)

          ! Compute auxt6 the local numbering of the interior objects in an object of dim od
          call memalloc(od,(2-1)**(od),auxt6,__FILE__,__LINE__)
          auxt6 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,od,auxt6,1,2-1)

       end if


       ! For each dimension, there are bnm(nd,od) objects up to translation
       do j = 1,bnm(nd,od)
          idm = -1 ! positions in which the nodes variates
          fdm = -1 ! Positions corresponding to the translation between paralel objects
          aux = -1 ! auxiliar vector to construct fdm from idm
          cd = cd+1

          ! Take the position that will vary inside the object
          do k = 1,od
             idm(k) = obdla(cd,k+1) 
          end do

          !Mark the positions already taken by idm
          aux = 0
          do k=1,od
             aux(idm(k)+1) = 1
          end do

          !Construct the array of orthogonal space wrt idm. 
          !It gives the translations for each paralel object
          c = 1
          do k=1,nd
             if (aux(k) == 0) then
                fdm(c) = k-1
                c = c+1
             end if
          end do

          ! Corner numbering
          ! Loop over the translations
          do l = 1,2**(nd-od) 

             ! Set orientation of the object
             co = co +1
             call Q_orientation_object(fefi%o(co),fdm(1:nd-od),od,nd,l)

             !ijk_g(jdm) will contain the translations from one object to another
             !ijk_g(idm) will contain the variations inside the object
             do k = 1,nd-od
                ijk_g(fdm(k)+1) = auxt1(k,l)
             end do

             ! Loop over the corners inside the object
             do m = 1,2**od
                do k = 1,od
                   ijk_g(idm(k)+1) = auxt2(k,m)
                end do
                c2 = c2+1
                fefi%crxob%l(c2) = Q_gijk(ijk_g,nd,1) !store the object numbering of the corner 
             end do
          end do

          ! Interior node numbering
          ! Loop over the translations
          do l = 1,2**(nd-od)
             !ijk_g(jdm) will contain the translations from 1 object to another; must be scaled by p
             !ijk_g(idm) will contain the variations inside the object
             do k = 1,nd-od
                ijk_g(fdm(k)+1) = auxt1(k,l)*p
             end do

             ! Loop over the interior nodes of the object
             do m = 1,(p-1)**(od)
                do k = 1,od
                   ijk_g(idm(k)+1) = auxt3(k,m)
                end do
                c3 = c3+1
                fefi%ndxob%l(c3) = Q_gijk(ijk_g,nd,p) ! Store the local numbering in fefi%ndxob%l
             end do
          end do

          ! All node numbering
          !Loop over the translations
          do l = 1,2**(nd-od)
             !ijk_g(jdm) will contain the translations from 1 object to another; must be scaled by p
             !ijk_g(idm) will contain the variations inside the object
             do k = 1,nd-od
                ijk_g(fdm(k)+1) = auxt1(k,l)*p
             end do

             ! Loop over the interior nodes of the object
             do m = 1,(p+1)**(od)
                do k = 1,od
                   ijk_g(idm(k)+1) = auxt4(k,m)
                end do
                c4 = c4+1
                fefi%ntxob%l(c4) = Q_gijk(ijk_g,nd,p) ! Store the local numbering in ntxob%l
             end do
          end do

          ! obxob array and auxiliar ob2node array 

          ! Interior node numbering
          ! Loop over the translations
          do l = 1,2**(nd-od)
             !ijk_g(jdm) will contain the translations from 1 object to another; must be scaled by p
             !ijk_g(idm) will contain the variations inside the object
             do k = 1,nd-od
                ijk_g(fdm(k)+1) = auxt1(k,l)*2
             end do

             ! Loop over the interior nodes of the object
             do m = 1,(2-1)**(od)
                do k = 1,od
                   ijk_g(idm(k)+1) = auxt6(k,m)
                end do
                c6 = c6+1
                ob2node(c6) = Q_gijk(ijk_g,nd,2) ! Store the local numbering in fefi%ndxob%l
             end do
          end do


          do l = 1,2**(nd-od)

             ! Fixed values for the object
             do k = 1,nd-od
                ijk_g(fdm(k)+1) = auxt1(k,l)*2
             end do

             ! Fill obxob (equivalent to ntxob for p=2)
             do m = 1,3**od
                do k = 1,od
                   ijk_g(idm(k)+1) = auxt5(k,m)
                end do
                c5 = c5 +1
                fefi%obxob%l(c5) = Q_gijk(ijk_g,nd,2)
             end do

             ! Define ijk_g for the node in the center of the object
             do k = 1,od
                ijk_g(idm(k)+1) = 1
             end do

             ! Find the generic node numbering
             m = Q_gijk(ijk_g,nd,2)

             ! Fill ob2node array
             !ob2node(m) = co
          end do

       end do

       !Deallocate auxiliar arrays
       if (od < nd) call memfree(auxt1,__FILE__,__LINE__)
       if (od >0) then
          call memfree(auxt2,__FILE__,__LINE__)
          if (p > 1) call memfree(auxt3,__FILE__,__LINE__)
          call memfree(auxt4,__FILE__,__LINE__)
          call memfree(auxt5,__FILE__,__LINE__)
          call memfree(auxt6,__FILE__,__LINE__)
       end if
    end do

    do i=1,no
       node2ob(ob2node(i)) = i
    end do

    ! Modify the identifiers of the nodes by the ids of the object in obxob
    do c5 = 1, nt2
       fefi%obxob%l(c5) = node2ob(fefi%obxob%l(c5))
    end do

    ! Sort the array 
    !do co = 1, no
    !   call sort(fefi%obxob%p(co+1)-fefi%obxob%p(co),fefi%obxob%l(fefi%obxob%p(co):fefi%obxob%p(co+1)))
    !end do

    do i=1,nn
       fefi%ndxob_int%l(i) = i
    end do

    ! Deallocate OBDLA
    call memfree(obdla,__FILE__,__LINE__)
    call memfree(ob2node,__FILE__,__LINE__)
    call memfree(node2ob,__FILE__,__LINE__)

    ! ! Create the face permutation of nodes
    ! if (nd>2) then call memalloc(2*2**2,fefi%nodes_obj(3),fefi%o2n,__FILE__,__LINE__)

    ! write(*,*) 'orientation objects'
    ! do od = 1,nd
    !    write(*,*) 'dime', od, '--------------------------'
    !    write(*,*) fefi%o(fefi%nobje_dim(od):fefi%nobje_dim(od+1)-1)
    ! end do
    ! write(*,*) 'no+1', no+1, 'ndxob%p'
    ! do od = 1,no+1
    !    write(*,*) fefi%ndxob%p(od), ', &'
    ! end do
    ! write(*,*) 'nn', nn, 'ndxob%l'
    ! do od = 1,nn
    !    write(*,*) fefi%ndxob%l(od), ', &'
    ! end do

    ! write(*,*) 'no+1', no+1, 'ntxob%p'
    ! do od = 1,no+1
    !    write(*,*) fefi%ntxob%p(od), ', &'
    ! end do
    ! write(*,*) 'nt', nt, 'ntxob%l'
    ! do od = 1,nt
    !    write(*,*) fefi%ntxob%l(od), ', &'
    ! end do

    ! write(*,*) 'no+1', no+1, 'crxob%p'
    ! do od = 1,no+1
    !    write(*,*) fefi%crxob%p(od), ', &'
    ! end do
    ! write(*,*) 'nc', nc, 'crxob%l'
    ! do od = 1,nc
    !    write(*,*) fefi%crxob%l(od), ', &'
    ! end do
  end subroutine Q_fixed_info_fill

  !==================================================================================================
  subroutine Q_orientation_object(o,fdm,od,nd,l)
    implicit none
    ! Parameters
    integer(ip), intent(out) :: o
    integer(ip), intent(in)  :: fdm(:)  ! fdm gives the orthogonal directions to the object
    integer(ip), intent(in)  :: od,nd,l ! l=translation ordering

    if (nd == 2 .and. od == 1) then
       o = modulo(l+fdm(1),2)
    elseif (nd == 3 .and. od == 2) then
       o = modulo(l+fdm(1),2)
    elseif (nd>3) then
       write(*,*) __FILE__,__LINE__,'WARNING!! the orientation is not defined for dimension >3'
    else
       o = 0
    end if
  end subroutine Q_orientation_object

  !==================================================================================================
  ! Given an orientation 'o' and a rotation 'r' with respect to a reference element,
  ! o2n(o_r,:) gives the local nodes in the refence element corresponding to the 
  subroutine Q_o2n_2d_create(o2n,p)
    implicit none
    ! Parameters
    integer(ip), allocatable, intent(inout) :: o2n(:,:)
    integer(ip)             , intent(in)    :: p

    ! Local variables
    integer(ip) :: o,r,o_r

    call memalloc(8,int((p+1)**2),o2n,__FILE__,__LINE__)
    o2n = 1

    ! Loop over possible orientations of a face
    do o = 0,1
       ! Loop over possible rotations of a face
       do r= 1,4  
          o_r = 4*o+r
          call Q_permute_or_2d(o2n(:,o_r),p,o,r)
       end do
    end do
  end subroutine Q_o2n_2d_create

  !=================================================================================================
  ! This subroutine gives the reodering (o2n) of the nodes of an object given an orientation 'o'
  ! and a delay 'r' wrt to a refence element sharing the same object.
  subroutine  Q_permute_or( o2n,p,o,r,nd )
    implicit none
    integer(ip), intent(in)    :: p,o,r,nd
    integer(ip), intent(inout) :: o2n(:)

    if     (nd == 0) then
       o2n(1) = 1
    elseif (nd == 1) then
       call Q_permute_or_1d(o2n(1:p+1),p,r)
    elseif (nd == 2) then
       call Q_permute_or_2d(o2n(1:int((p+1)**2)),p,o,r)
    else
       write(*,*) __FILE__,__LINE__,'WARNING! Permutations not given for nd>3'
    end if
  end subroutine Q_permute_or

  !=================================================================================================
  subroutine  Q_permute_or_1d( o2n,p,r )
    implicit none
    integer(ip), intent(in)    :: p,r
    integer(ip), intent(inout) :: o2n(p+1)

    ! Local variables
    integer(ip) :: i

    ! Generic loop+rotation identifier  
    if (r==1) then
       o2n = (/(i,i=1,p+1)/)
    elseif (r==2) then
       o2n = (/(p+1-i,i=0,p)/)
    else
       write(*,*) __FILE__,__LINE__,'Q_permute_or_1d:: ERROR! Delay cannot be >1 for edges'
    end if
  end subroutine Q_permute_or_1d

  !==================================================================================================
  subroutine Q_permute_or_2d( o2n,p,o,r )
    implicit none
    integer(ip), intent(in)    :: p,o,r
    integer(ip), intent(inout) :: o2n((p+1)**2)

    ! Local variables
    integer(ip) :: o_r,i,j,ij_q(4) ! ij_q = (i,j,p-i,p-j)
    integer(ip) :: ij_n(2),go
    integer(ip) :: ij_perm_quad(2,8) = reshape((/ 1, 2, 2, 3, 4, 1, 3, 4, 2, 1, 3, 2, 1, 4, 4, 3/), &
         &                                     (/2,8/))

    ! Generic loop+rotation identifier
    o_r = 4*o+r
    do j = 0,p
       ij_q(2) = j
       ij_q(4) = p-j
       do i = 0,p
          ij_q(1) = i
          ij_q(3) = p-i
          ! Get the global numbering of node (i,j)
          go = Q_gijk(ij_q(1:2),2,p)
          ! i,j coordinates for the o_r permutation
          ij_n(1:2) = ij_q(ij_perm_quad(1:2,o_r)) 
          ! Get the global numbering of node ij_n
          o2n(go) = Q_gijk(ij_n,2,p)
       end do
    end do
  end subroutine Q_permute_or_2d

  !==================================================================================================
  ! Given the dimension nd of the quad and the positions p0,...,p1 we want to consider is going
  ! to take values. The routine Q_R_IJK returns the number of nodes (co) and a matrix auxt1 that 
  ! for the nth node of the quad it gives the ijk position that it corresponds to
  !         | p0 p0   ... p0 p0   ... p1 |
  ! auxt1 = | p0 p0   ... p0 p0+1 ... p1 |
  !         | p0 p0+1 ... p1 p0   ... p1 |
  recursive subroutine Q_r_ijk(co,d,ijk,nd,auxt1,p0,p1)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: d,nd,p0,p1
    integer(ip), intent(inout) :: ijk(nd),co              !#nodes
    integer(ip), intent(inout) :: auxt1(nd,(p1+1-p0)**nd) ! ode ijk position matrix

    ! Local variables
    integer(ip)                :: dp,i,ip

    do ip = p0,p1
       if (d>nd) exit
       ijk(nd-d+1) = ip 
       if (d < nd) call Q_r_ijk(co,d+1,ijk,nd,auxt1,p0,p1)
       if (d == nd) then
          co = co + 1
          do i = 1,nd
             auxt1(i,co) = ijk(i)
          end do
       end if
    end do
  end subroutine Q_r_ijk

  !==================================================================================================
  ! R_DIM construct OBDLA matrix
  ! For each object, up to a displacement, we have an identifier id={1..nod}
  ! obdla(id,1) = dimension of the object
  ! obdla(id,2:obdla(id,1)+1) = gives the directions that define the object
  recursive subroutine r_dim(co,idm,ko,i,nd,od,obdla,nod)
    implicit none
    integer(ip), intent(in)    :: ko !space position we begin to count from (ij=>i<j)
    integer(ip), intent(in)    :: i  !i=local space position we are currently modifying (i=0..od)
    integer(ip), intent(in)    :: nd,od,nod
    integer(ip), intent(inout) :: obdla(nod,nd+1)
    integer(ip), intent(inout) :: idm(od) 
    integer(ip), intent(inout) :: co      ! Pointer to the position of the object
    integer(ip)                :: aux,ijk_c(nd),j,cd,kn,s

    !Given dimension od of the object
    do kn = ko,nd-od+i
       idm(i+1) = kn 
       if (i+1 < od) then
          call r_dim(co,idm,kn+1,i+1,nd,od,obdla,nod)
       else
          obdla(co,1) = od
          do j = 1,od
             obdla(co,j+1) = idm(j)
          end do
          co = co + 1
       end if
    end do
  end subroutine r_dim

  !=================================================================================================
  ! Q_NNODS(k,p) computes the #nodes in a simplex of dim k and order p
  integer(ip) function Q_nnods(k,p)
    implicit none
    integer(ip), intent(in)    :: k, p

    Q_nnods = int((p+1)**k)
  end function Q_nnods

  !==================================================================================================
  ! FC(k)=k! computes the factorial of k 
  integer(ip) function fc(i)
    implicit none
    integer(ip) :: i, k
    fc = 1
    do k=2,i
       fc = fc*k
    end do
  end function fc

  !=================================================================================================
  ! BNM(A,B)=A!/((A-B)!B!) computes the binomial coefficient of (A,B), A>B
  integer (ip) function bnm(a,b)
    implicit none
    integer(ip) :: a,b
    if (a >= b) then
       bnm = int(fc(a)/(fc(b)*fc(a-b)))
    else
       write(*,*) 'ERROR: no binomial coef for b > a'
       check(.false.)
    end if
  end function bnm

  !==================================================================================================
  ! Q_GIJK(i,nd,p) returns the generic identifier of a node with coordinates i in an elem(nd,p)
  !Given the coordinates ijk (in i) and the dimension nd and the order p of the 
  !quad, it returns the local numbering of the node: Q_gijk=i+j*(p+1)+k*(p+1)^2
  integer(ip) function Q_gijk(i,nd,p)
    implicit none
    integer(ip) :: nd,i(nd),p,k
    Q_gijk = 1
    do k = 1,nd
       Q_gijk = Q_gijk + i(k)*((p+1)**(k-1))
    end do
  end function  Q_gijk

  !==================================================================================================
  ! Q_IJKG(i,g,nd,p) returns coordinates of the g-th node in an elem(nd,p)
  subroutine Q_ijkg(i,g,nd,p)
    implicit none
    integer(ip), intent(in)  :: nd,g,p
    integer(ip), intent(out) :: i(nd)

    integer(ip)              :: k,g2

    g2 = g-1
    do k=1,nd
       i(k) = int(mod(g2,(p+1)**k)/(p+1)**(k-1))
       g2 = g2 - i(k)*(p+1)**(k-1)
    end do
  end subroutine Q_ijkg

  !==================================================================================================
  subroutine Q_set_integ(ndime,order,nnode,ngaus,nlocs,lrule,llapl,mnode)
    implicit none
    ! Parameters
    integer(ip)          , intent(in)  :: ndime,order
    integer(ip)          , intent(out) :: nnode,ngaus,nlocs,lrule,llapl
    integer(ip), optional, intent(in)  :: mnode

    ! Second order derivative storage
    if (order > 1) then
       llapl = 1
    else
       llapl = 0
    end if

    ! Quadrilateral elements ---> tensorial product
    !ngaus = order + 1
    if (present(mnode)) then
       ngaus = mnode
       nlocs = int(real(mnode)**(1.0_rp/real(ndime)))
    else
       nlocs = order + 1
       ngaus = (order + 1)**ndime
    end if
    nnode = (order + 1)**ndime
    lrule = ruqope_id
    !lrule = ruqope_tp_id
  end subroutine Q_set_integ

  !==================================================================================================
  ! Coord(:,g) = coordinates of the g-th node in the reference element
  subroutine Q_refcoord (coord,nd,p,nn)
    implicit none
    ! Parameters
    integer(ip), intent(in)    :: nd,p,nn
    real(rp),    intent(inout) :: coord(nd,nn)

    ! Local variables
    integer(ip) :: ijk(nd),i,d
    real(rp)    :: c1d(p+1)

    call Q_coord_1d(c1d,p+1)

    do i=1,nn

       call Q_ijkg(ijk,i,nd,p)
       do d=1,nd
          coord(d,i) = c1d(ijk(d)+1)
       end do
    end do
  end subroutine Q_refcoord

  ! =================================================================================================
  ! Set n equidistant points in the segment [-1,1]
  subroutine Q_coord_1d (x,n)
    implicit none
    ! Parameters
    integer(ip), intent(in)  :: n
    real(rp)   , intent(out) :: x(n)

    ! Local variables
    integer(ip)              :: i

    do i = 0,n-1
       x(i+1) = 2*real(i)/(real(n)-1)-1
    end do
  end subroutine Q_coord_1d

  ! =================================================================================================
  subroutine Q_face_outno(ndime,nlocf,outno)
    implicit none
    ! Parameters
    real(rp),    intent(out) :: outno(ndime,2*ndime,nlocf)
    integer(ip), intent(in)  :: ndime,nlocf

    ! Local variables
    logical(lg)                  :: assig = .false.
    integer(ip)              :: ilocs,i,j,iface 

    ! Initialize
    outno = 0

    iface = 1
    do i=0,ndime-1
       do j = 1,2
          do ilocs = 1, nlocf
             outno(ndime-i,iface,ilocs) = int((-1)**j)
          end do
          iface = iface + 1
       end do
    end do
  end subroutine Q_face_outno

  integer(ip) function get_order(type, ndime, nnodes)
    implicit none
    integer(ip) :: type, ndime, nnodes

    if ( type == P_type_id ) then
       do get_order=1,max_order
          if ( P_nnods(ndime,get_order) == nnodes ) exit                
       end do
    else if ( type == Q_type_id ) then
       do get_order=1,max_order
          if ( int((get_order+1)**ndime) == nnodes ) exit                
       end do
    end if
    assert ( get_order == max_order + 1 )

  end function get_order

end module fem_space_types



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
module fem_block_matrix_class
  use types
  use memor
  use block_elmat_class
  use fem_matrix_class
  use fem_space_class
  use elmat_class
  use dof_handler_class
  use fem_blocks_class
  implicit none
# include "debug.i90"

  private

  ! Pointer to matrix
  type p_fem_matrix
    type(fem_matrix), pointer :: p_f_matrix
  end type p_fem_matrix

  ! Block Matrix
  type fem_block_matrix
    integer(ip)                     :: nblocks
    type(p_fem_matrix), allocatable :: blocks(:,:)
  end type fem_block_matrix

  interface fem_block_matrix_assembly
     module procedure fem_block_matrix_assembly_standard
     module procedure fem_block_matrix_assembly_w_dof_handler
     module procedure fem_block_matrix_assembly_nosq_w_dof_handler
  end interface fem_block_matrix_assembly

  ! Types
  public :: fem_block_matrix

  ! Functions
  public :: fem_block_matrix_alloc, fem_block_matrix_alloc_block,       & 
            fem_block_matrix_set_block_to_zero, fem_block_matrix_print, & 
            fem_block_matrix_free, fem_block_matrix_assembly,           & 
            fem_block_matrix_info, fem_block_matrix_zero

contains

  !=============================================================================
  subroutine fem_block_matrix_alloc(nblocks, bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: nblocks
    type(fem_block_matrix), intent(out) :: bmat
    integer(ip) :: ib,jb

    bmat%nblocks = nblocks
    allocate ( bmat%blocks(nblocks,nblocks) )
    do ib=1, nblocks 
      do jb=1, nblocks
           allocate ( bmat%blocks(ib,jb)%p_f_matrix )
      end do
    end do
  end subroutine fem_block_matrix_alloc

  subroutine fem_block_matrix_alloc_block (ib,jb,bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: ib,jb
    type(fem_block_matrix), intent(inout) :: bmat
    if ( .not. associated( bmat%blocks(ib,jb)%p_f_matrix)) then
       allocate ( bmat%blocks(ib,jb)%p_f_matrix )
    end if
  end subroutine

  subroutine fem_block_matrix_set_block_to_zero (ib,jb,bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)  :: ib,jb
    type(fem_block_matrix), intent(inout) :: bmat

    if ( associated(bmat%blocks(ib,jb)%p_f_matrix) ) then
       deallocate (bmat%blocks(ib,jb)%p_f_matrix)
       nullify ( bmat%blocks(ib,jb)%p_f_matrix )
    end if
  end subroutine

  subroutine fem_block_matrix_zero(bmat)
    implicit none
    ! Parameters
    type(fem_block_matrix), intent(inout) :: bmat
    integer(ip)                           :: ib,jb

    do ib=1,bmat%nblocks
       do jb=1,bmat%nblocks
          call fem_matrix_zero (bmat%blocks(ib,jb)%p_f_matrix)
       end do
    end do

  end subroutine fem_block_matrix_zero


  subroutine fem_block_matrix_print (lunou, f_b_matrix)
    implicit none
    type(fem_block_matrix), intent(in)    :: f_b_matrix
    integer(ip)           , intent(in)    :: lunou
    integer(ip)                           :: i

  end subroutine fem_block_matrix_print

  subroutine fem_block_matrix_free(f_b_matrix)
    implicit none
    type(fem_block_matrix), intent(inout) :: f_b_matrix
    integer(ip) :: ib,jb

    do ib=1, f_b_matrix%nblocks 
       do jb=1, f_b_matrix%nblocks
          if ( associated(f_b_matrix%blocks(ib,jb)%p_f_matrix) ) then
             deallocate (f_b_matrix%blocks(ib,jb)%p_f_matrix) 
          end if
       end do
    end do

    f_b_matrix%nblocks = 0
    deallocate ( f_b_matrix%blocks ) 
  end subroutine fem_block_matrix_free

  !=============================================================================
  subroutine fem_block_matrix_assembly_standard(nn,ln,bea,bmat)
    implicit none
    ! Parameters
    integer(ip)           , intent(in)    :: nn
    integer(ip)           , intent(in)    :: ln(:)
    type(block_elmat)     , intent(in)    :: bea
    type(fem_block_matrix), intent(inout) :: bmat

    ! Locals
    integer(ip) :: ib, jb
    ! call block_elmat_print (6, bea) ! DBG:
    do ib=1, bmat%nblocks
      do jb=1, bmat%nblocks
         if ( associated(bmat%blocks(ib,jb)%p_f_matrix) ) then
            call fem_matrix_assembly (nn,ln,bea%blocks(ib,jb),bmat%blocks(ib,jb)%p_f_matrix)
         end if
      end do
    end do

  end subroutine fem_block_matrix_assembly_standard

!=============================================================================
  subroutine fem_block_matrix_assembly_w_dof_handler(nn,dofh,el,l2g1,l2g2,l2g3,l2g4,ea,mat)
    implicit none
    integer(ip) ,           intent(in)    :: nn(:)
    type(dof_handler),      intent(in)    :: dofh
    type(fem_element),      intent(in)    :: el
    type(elmat) ,           intent(in)    :: ea
    type(fem_block_matrix), intent(inout) :: mat
    integer(ip),            intent(inout) :: l2g1(:),l2g2(:),l2g3(:),l2g4(:)
    ! Locals
    integer(ip)              :: nd1,nd2,inode,i,j,k,ibloc,jbloc,ib(3),idofs,jdofs,maxnode
    integer(ip), allocatable :: jb(:) 

    ! Asserts
    nd1=size(ea%data,dim=1)
    nd2=size(ea%data,dim=2)
    assert(nd1 == 1)
    assert(nd2 == 1)
    assert(size(nn,1)==el%nint) 

    maxnode=size(el%jvars,1)       

    do ibloc = 1,dofh%blocks%nb

       ! Construct block structures
       ib(1) = 1
       ib(2) = ib(1) + dofh%blknvarsxprob(ibloc)%a(el%prob) ! Rows

       ! ROW DOFS
       idofs = 0; l2g1 = 0; l2g3 = 0
       do j = dofh%i_varsxprob(ibloc)%a(el%prob),dofh%i_varsxprob(ibloc)%a(el%prob+1)-1
          do inode = 1,nn(el%iv(j))
             idofs = idofs+1
             l2g1(idofs) = el%elem2dof(inode,dofh%j_varsxprob(ibloc)%a(j)) ! idof
             if(l2g1(idofs)/=0) l2g3(idofs) = dofh%j_varsxprob(ibloc)%a(j) 
          end do
       end do

       do jbloc = 1,dofh%blocks%nb

          if ( associated(mat%blocks(ibloc,jbloc)%p_f_matrix)) then

          ! Assert
          assert(nd1==mat%blocks(ibloc,jbloc)%p_f_matrix%nd1)
          assert(nd2==mat%blocks(ibloc,jbloc)%p_f_matrix%nd2)

          ! Construct block structures
          ib(3) = ib(2) + dofh%blknvarsxprob(jbloc)%a(el%prob) ! Columns
          call memalloc(ib(3)-1,jb,__FILE__,__LINE__)
          jb(1:(ib(2)-1)) = dofh%j_varsxprob(ibloc)%a(dofh%i_varsxprob(ibloc)%a(el%prob) : &
               &                                     (dofh%i_varsxprob(ibloc)%a(el%prob+1)-1))
          jb(ib(2):(ib(3)-1)) = dofh%j_varsxprob(jbloc)%a(dofh%i_varsxprob(jbloc)%a(el%prob): &
               &                                         (dofh%i_varsxprob(jbloc)%a(el%prob+1)-1))

          ! COLUMN DOFS
          jdofs = 0; l2g2 = 0; l2g4 = 0
          do j = dofh%i_varsxprob(jbloc)%a(el%prob),dofh%i_varsxprob(jbloc)%a(el%prob+1)-1
             do inode = 1,nn(el%iv(j))
                jdofs = jdofs+1
                l2g2(jdofs) = el%elem2dof(inode,dofh%j_varsxprob(jbloc)%a(j)) ! jdof
                if(l2g2(jdofs)/=0) l2g4(jdofs) = dofh%j_varsxprob(jbloc)%a(j)
             end do
          end do

          if (mat%blocks(ibloc,jbloc)%p_f_matrix%storage == blk) then ! gr must describe block sparsity pattern of mat
             if(mat%blocks(ibloc,jbloc)%p_f_matrix%type==css_mat) then

             else if(mat%blocks(ibloc,jbloc)%p_f_matrix%type==csr_mat) then

             else if(mat%blocks(ibloc,jbloc)%p_f_matrix%type==csc_mat) then

             end if
          else if (mat%blocks(ibloc,jbloc)%p_f_matrix%storage == scal) then ! gr must describe scalar sparsity pattern of mat
             if( mat%blocks(ibloc,jbloc)%p_f_matrix%type==css_mat ) then

             else if(mat%blocks(ibloc,jbloc)%p_f_matrix%type==csr_mat) then
                call ass_csr_blkmat_scal_w_dof_handler(mat%blocks(ibloc,jbloc)%p_f_matrix%symm,el%nint,  &
                     &                                 nn,size(el%p_nod,1),idofs,jdofs,el%ldof,ib,jb,    &
                     &                                 size(el%iv,1),el%iv,el%p_nod,l2g1,l2g2,l2g3,l2g4, &
                     &                                 ea%data,mat%blocks(ibloc,jbloc)%p_f_matrix%gr%nv, &
                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%gr%ia,         &
                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%gr%ja,         &
                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%a,             &
                     &                                 dofh%nvars,dofh%dof_coupl,maxnode,el%jvars)
!!$                call ass_csr_blkmat_scal_w_dof_handler(mat%blocks(ibloc,jbloc)%p_f_matrix%symm,nn(1),    &
!!$                     &                                 dofh%nvarsxprob(el%prob),ib,jb,l2g1,l2g2,ea%data, &
!!$                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%gr%nv,         &
!!$                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%gr%ia,         &
!!$                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%gr%ja,         &
!!$                     &                                 mat%blocks(ibloc,jbloc)%p_f_matrix%a)
             else if(mat%blocks(ibloc,jbloc)%p_f_matrix%type==csc_mat) then

             end if
          end if

          call memfree(jb)

       end if
    end do
 end do

  end subroutine fem_block_matrix_assembly_w_dof_handler

  !=============================================================================
  subroutine fem_block_matrix_assembly_nosq_w_dof_handler(nn,dofh,el,l2g1,l2g2,l2g3,l2g4,bl,ea,mat)
    implicit none
    integer(ip) ,           intent(in)    :: nn(:)
    type(dof_handler),      intent(in)    :: dofh
    type(fem_element),      intent(in)    :: el
    type(fem_blocks),       intent(in)    :: bl
    type(elmat) ,           intent(in)    :: ea
    type(fem_block_matrix), intent(inout) :: mat
    integer(ip),            intent(inout) :: l2g1(:),l2g2(:),l2g3(:),l2g4(:)
    ! Locals
    integer(ip)              :: nd1,nd2,inode,i,j,k,ibloc,jbloc,ib(3)
    integer(ip), allocatable :: jb(:) ! nn*dofh%nvarsxprob(el%prob))

    write (*,*) 'fem_block_matrix_assembly_nosq_w_dof_handler TO BE DONE'

  end subroutine fem_block_matrix_assembly_nosq_w_dof_handler

  subroutine fem_block_matrix_info ( f_blk_mat, me, np )
    implicit none

    ! Parameters 
    type(fem_block_matrix) , intent(in)    :: f_blk_mat
    integer                , intent(out)   :: me
    integer                , intent(out)   :: np
    
    me = 0
    np = 1 
  end subroutine fem_block_matrix_info

  !============================================================================
  subroutine ass_csr_blkmat_scal_w_dof_handler(symm,nint,nn,nd,id,jd,ld,ib,jb,nva,iv,pn, &
       &                                       l2g1,l2g2,l2g3,l2g4,ea,nv,ia,ja,la,nvar,  &
       &                                       dofc,mn,jbn)

    implicit none
    integer(ip) , intent(in)      :: symm
    integer(ip) , intent(in)      :: nint, nv, nd, nva, id, jd, ld, nvar, mn
    integer(ip) , intent(in)      :: nn(nint),pn(nd),iv(nva)
    integer(ip) , intent(in)      :: ib(3),jb(ib(3)-1)
    integer(ip) , intent(in)      :: l2g1(id), l2g2(jd)
    integer(ip) , intent(in)      :: l2g3(id), l2g4(jd)
    integer(ip) , intent(in)      :: ia(nv+1), ja(ia(nv+1)-1)
    integer(ip) , intent(in)      :: dofc(nvar,nvar)
    integer(ip) , intent(in)      :: jbn(mn,nva)
    real(rp)    , intent(in)      :: ea(ld,ld)
    real(rp)    , intent(inout)   :: la(ia(nv+1)-1)

    ! local variables
    integer(ip)                   :: in, jn, ip, jp, il, ig, jl, jg, k, ie, je, ivar, jvar

    ! CSR
    il = 0
    do ip = ib(1),ib(2)-1
       do in = 1,nn(iv(jb(ip)))
          il = il + 1
          ig = l2g1(il)
          ivar = l2g3(il)
          if (ig /= 0) then
             ie = pn(in)-1 + jbn(in,jb(ip))
             jl = 0
             do jp = ib(2),ib(3)-1
                do jn = 1,nn(iv(jb(jp)))
                   jl = jl + 1
                   jg = l2g2(jl)
                   jvar = l2g4(jl)
                   if (jg /= 0) then
                      if(dofc(ivar,jvar)==1) then
                         je = pn(jn)-1 + jbn(jn,jb(jp))
                         do k = ia(ig),ia(ig+1)-1
                            if(ja(k) == jg) then
                               la(k) = la(k) + ea(ie,je)
                               exit
                            end if
                         end do
                      end if
                   end if
                end do
             end do
          end if
       end do
    end do

    ! CSS
    ! idem for 
    ! jl=il,nl
    ! k = ia(min(ig,jg)) 
    ! do while(ja(k)/= max(ig,jg))

  end subroutine ass_csr_blkmat_scal_w_dof_handler

end module fem_block_matrix_class

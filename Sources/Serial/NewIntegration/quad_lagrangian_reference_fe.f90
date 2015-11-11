module quad_lagrangian_reference_fe_names
  use reference_fe_names
  use quadrature_names
  use allocatable_array_ip1_names
  use types_names
  use memor_names
  implicit none
# include "debug.i90"

  private

  type, extends(reference_fe_t) :: quad_lagrangian_reference_fe_t
  private
contains 
  procedure :: create
  procedure :: create_interpolation
  procedure :: set_integration_rule
  procedure :: create_quadrature
  procedure :: fill
  !procedure :: local_to_ijk_node     
  !procedure :: ijk_to_local_node     
end type quad_lagrangian_reference_fe_t

public :: quad_lagrangian_reference_fe_t

contains 

subroutine create ( this, number_dimensions, order, continuity )
 implicit none 
 class(quad_lagrangian_reference_fe_t), intent(out) :: this 
 integer(ip), intent(in)  :: number_dimensions, order
 logical, optional, intent(in) :: continuity

 call this%set_common_data ( number_dimensions, order, continuity )
 call this%set_topology( "tet" )
 call this%set_fe_type( "Lagrangian" )
 call this%fill( )

end subroutine create

subroutine set_integration_rule ( this , quadrature, interpolation )
 implicit none 
 class(quad_lagrangian_reference_fe_t), intent(in) :: this
 class(quadrature_t), intent(in) :: quadrature
 class(interpolation_t), intent(out) :: interpolation

 ! Here we should put all the things in interpolation.f90
end subroutine set_integration_rule

subroutine create_interpolation ( this, quadrature, interpolation )
 implicit none 
 class(quad_lagrangian_reference_fe_t), intent(in) :: this 
 class(quadrature_t), intent(in) :: quadrature
 type(interpolation_t), intent(out) :: interpolation
end subroutine create_interpolation

subroutine create_quadrature ( this, quadrature, max_order )
 implicit none 
 class(quad_lagrangian_reference_fe_t), intent(in) :: this        
 integer(ip), optional, intent(in) :: max_order
 class(quadrature_t), intent(out) :: quadrature
end subroutine create_quadrature

subroutine fill (this)
 implicit none
 ! Parameters
 class(quad_lagrangian_reference_fe_t), intent(inout) :: this

 ! Local variables
 integer(ip)               :: nd,p
 integer(ip)               :: i,j,k,l,m
 integer(ip)               :: aux1,aux2,aux3,aux4
 integer(ip)               :: od,cd,kk,c
 integer(ip)               :: c2,c3,c4,c5,c6,co
 integer(ip)               :: no  ! #vefs in nd-quad
 integer(ip)               :: nn  ! #nodes in the nd-quad of order p
 integer(ip)               :: nt  ! sum of the #nodes (all, not only interior) of the vefs
 integer(ip)               :: nt2 ! sum of the #vefs (all, not only interior) of the vefs
 integer(ip)               :: nc  ! #corners x #{vefs delimited by each corner}
 integer(ip)               :: nod ! #corners in the quad. nod = 2^(nd)
 integer(ip), allocatable  :: auxt1(:,:),auxt2(:,:),auxt3(:,:),auxt4(:,:),auxt5(:,:), auxt6(:,:)
 integer(ip), allocatable  :: obdla(:,:),node2ob(:),ob2node(:)

 integer(ip), allocatable  :: aux(:),idm(:),fdm(:),ijk(:),ijk_g(:)
 
 integer(ip), pointer :: number_vefs, number_nodes, number_vefs_dimension(:)
 type(allocatable_array_ip1_t), pointer ::  orientation
 type(list_t), pointer :: interior_nodes_vef, nodes_vef, corners_vef, vefs_vef

 number_vefs => this%get_pointer_number_vefs()
 number_nodes => this%get_pointer_number_nodes()
 number_vefs_dimension => this%get_pointer_number_vefs_dimension()
 orientation => this%get_pointer_orientation()
 interior_nodes_vef => this%get_pointer_interior_nodes_vef()
 interior_nodes_vef => this%get_pointer_interior_nodes_vef()
 nodes_vef => this%get_pointer_nodes_vef()
 corners_vef => this%get_pointer_corners_vef()
 vefs_vef => this%get_pointer_vefs_vef()

 !  ! Initilize values
 nd = this%get_number_dimensions()
 p  = this%get_order()

 call memalloc( nd, aux, __FILE__, __LINE__ )
 call memalloc( nd, idm, __FILE__, __LINE__ )
 call memalloc( nd, fdm, __FILE__, __LINE__ )
 call memalloc( nd, ijk, __FILE__, __LINE__ )
 call memalloc( nd, ijk_g, __FILE__, __LINE__ )

 no = 0
 nn = 0
 nt = 0
 nc = 0
 nt2 = 0
 nod = 0

 ! Initialize nvef_dim, nodes_vef
 !call memalloc(nd+2,nvef_dim,__FILE__,__LINE__)
 !call memalloc(nd+1,nodes_vef,__FILE__,__LINE__)
 number_vefs_dimension = 0
 number_vefs_dimension(1) = 1

 do k = 0,nd
    i = int(2**(nd-k)*bnm(nd,k)) ! #vefs of dimension k
    no = no + i                  ! compute #vefs
    nn = nn + int(i*((p-1)**k))  ! #nodes inside vefs of dimension k
    nt = nt + int(i*((p+1)**k))  ! nodes in the clousure of the vef
    nt2 = nt2 + int(i*((3)**k))! vefs in the clousure of the vef
    nc = nc + int(i*2**k)        ! corners delimiting vefs of dimension k
    nod = nod + bnm(nd,k)        ! #nodes/{displacement} (2^n = sum(bnm(n,k)), k=0,..,n)
    ! Pointer to obj id by dim. Local obj of dim k are nvef_dim(k):nvef_dim(k+1)
    number_vefs_dimension(k+2) = number_vefs_dimension(k+1) + i 
    ! #nodes in vefs of dimension k 
 end do

 ! Set constant values of reference_element
 number_vefs = no-1
 number_nodes = int((p+1)**nd) 

 ! Allocate arrays
 call orientation%create(no)
 call memalloc(no+1,interior_nodes_vef%p,__FILE__,__LINE__)  !Pointer to interior_nodes_vef%l for each vef
 call memalloc(nn,  interior_nodes_vef%l,__FILE__,__LINE__)  !Array of interior nodes of each vef
 call memalloc(no+1,nodes_vef%p,__FILE__,__LINE__)  !Pointer to nodes_vef%l for each vef
 call memalloc(nt,  nodes_vef%l,__FILE__,__LINE__)  !Array of all nodes of each vef
 call memalloc(no+1,vefs_vef%p,__FILE__,__LINE__)  !Pointer to vefs_vef%l for each vef
 call memalloc(nt2,  vefs_vef%l,__FILE__,__LINE__)  !Array of all vefs of each vef
 call memalloc(no+1,corners_vef%p,__FILE__,__LINE__)  !Pointer to corners_vef%l for each vef
 call memalloc(nc,  corners_vef%l,__FILE__,__LINE__)  !Array of corners for each vef
 call memalloc(nod,nd+1,obdla,__FILE__,__LINE__)
 call memalloc(no, node2ob,__FILE__,__LINE__)        ! Auxiliar array
 call memalloc(no, ob2node,__FILE__,__LINE__)        ! Auxiliar array

 interior_nodes_vef%p=0   !Pointer to ndxob%l for each vef
 interior_nodes_vef%l=0   !Array of interior nodes of each vef
 nodes_vef%p=0   !Pointer to ntxob%l for each vef
 vefs_vef%l=0   !Array of all nodes of each vef
 vefs_vef%p=0   !Pointer to obxob%l for each vef
 nodes_vef%l=0   !Array of all nodes of each vef
 corners_vef%p=0   !Pointer to crxob%l for each vef
 corners_vef%l=0   !Array of corners for each vef


 !Initialize pointers
 interior_nodes_vef%p(1) = 1
 nodes_vef%p(1) = 1
 vefs_vef%p(1) = 1
 corners_vef%p(1) = 1

 !Loop over dimensions
 do k = 0,nd
    aux1 = int(((p-1)**k)) ! interior nodes for an vef of dim k
    aux3 = int(((p+1)**k)) ! Total nodes for an vef of dim k
    aux2 = int(2**k)       ! Corners for an vef of dim k
    aux4 = int((3**k)) ! Corners for an vef of dim k (idem p=2)

    ! Loop over vefs of dimension k
    do i = number_vefs_dimension(k+1),number_vefs_dimension(k+2)-1 
       interior_nodes_vef%p(i+1) = interior_nodes_vef%p(i) + aux1 ! assign pointers
       nodes_vef%p(i+1) = nodes_vef%p(i) + aux3 ! assign pointers
       corners_vef%p(i+1) = corners_vef%p(i) + aux2 ! assign pointers
       vefs_vef%p(i+1) = vefs_vef%p(i) + aux4 ! assign pointers 
    end do
 end do

 ! Initialize auxiliar values
 k = 0
 i = 0
 idm = 0
 j = 2

 ! Construction of obdla matrix
 ! For each vef, up to a displacement, we have an identifier id={1..nod}
 ! obdla(id,1) = dimension of the vef
 ! obdla(id,2:obdla(id,1)+1) = gives the directions that define the vef
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
 c2  = 0 ! ndxob%p counter
 c3  = 0 ! crxob%p counter
 c4  = 0 ! ntxob%p counter
 c5  = 0 ! obxob%p counter
 c6  = 0 ! ob2node   counter
 co  = 0 ! counter of vefs

 ! Loop over vefs dimensions
 do od = 0,nd
    ! Create ijk tables (od)
    ! Compute auxt1 the local numbering of the corners of an vef of dimension nd-od
    ! It allows to know how many translations for each paralel set of vefs
    if (od < nd) then
       call memalloc(nd-od,2**(nd-od),auxt1,__FILE__,__LINE__)
       auxt1 = 0
       kk    = 0
       aux1  = 1
       ijk   = 0
       call Q_r_ijk(kk,aux1,ijk,nd-od,auxt1,0,1)
    end if
    ! Compute auxt2 the local numbering of the corners in an vef of dim od
    if (od >0) then
       call memalloc(od,2**(od),auxt2,__FILE__,__LINE__)
       auxt2 = 0
       kk    = 0
       aux1  = 1
       ijk   = 0
       call Q_r_ijk(kk,aux1,ijk,od,auxt2,0,1)

       if (p > 1) then
          ! Compute auxt3 the local numbering of the interior nodes in an vef of dim od
          call memalloc(od,(p-1)**(od),auxt3,__FILE__,__LINE__)
          auxt3 = 0
          kk    = 0
          aux1  = 1
          ijk   = 0
          call Q_r_ijk(kk,aux1,ijk,od,auxt3,1,p-1)
       end if

       ! Compute auxt4 the local numbering of all nodes in an vef of dim od
       call memalloc(od,(p+1)**(od),auxt4,__FILE__,__LINE__)
       auxt4 = 0
       kk    = 0
       aux1  = 1
       ijk   = 0
       call Q_r_ijk(kk,aux1,ijk,od,auxt4,0,p)

       ! Compute auxt5 the local numbering of all nodes in an vef of dim od
       call memalloc(od,(2+1)**(od),auxt5,__FILE__,__LINE__)
       auxt5 = 0
       kk    = 0
       aux1  = 1
       ijk   = 0
       call Q_r_ijk(kk,aux1,ijk,od,auxt5,0,2)

       ! Compute auxt6 the local numbering of the interior vefs in an vef of dim od
       call memalloc(od,(2-1)**(od),auxt6,__FILE__,__LINE__)
       auxt6 = 0
       kk    = 0
       aux1  = 1
       ijk   = 0
       call Q_r_ijk(kk,aux1,ijk,od,auxt6,1,2-1)

    end if


    ! For each dimension, there are bnm(nd,od) vefs up to translation
    do j = 1,bnm(nd,od)
       idm = -1 ! positions in which the nodes variates
       fdm = -1 ! Positions corresponding to the translation between paralel vefs
       aux = -1 ! auxiliar vector to construct fdm from idm
       cd = cd+1

       ! Take the position that will vary inside the vef
       do k = 1,od
          idm(k) = obdla(cd,k+1) 
       end do

       !Mark the positions already taken by idm

       aux = 0
       do k=1,od
          aux(idm(k)+1) = 1
       end do

       !Construct the array of orthogonal space wrt idm. 
       !It gives the translations for each paralel vef
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

          ! Set orientation of the vef
          co = co +1
          call Q_orientation_vef(orientation%a(co),fdm(1:nd-od),od,nd,l)

          !ijk_g(jdm) will contain the translations from one vef to another
          !ijk_g(idm) will contain the variations inside the vef
          do k = 1,nd-od
             ijk_g(fdm(k)+1) = auxt1(k,l)
          end do

          ! Loop over the corners inside the vef
          do m = 1,2**od
             do k = 1,od
                ijk_g(idm(k)+1) = auxt2(k,m)
             end do
             c2 = c2+1
             corners_vef%l(c2) = Q_gijk(ijk_g,nd,1) !store the vef numbering of the corner 
          end do
       end do

       ! Interior node numbering
       ! Loop over the translations
       do l = 1,2**(nd-od)
          !ijk_g(jdm) will contain the translations from 1 vef to another; must be scaled by p
          !ijk_g(idm) will contain the variations inside the vef
          do k = 1,nd-od
             ijk_g(fdm(k)+1) = auxt1(k,l)*p
          end do

          ! Loop over the interior nodes of the vef
          do m = 1,(p-1)**(od)
             do k = 1,od
                ijk_g(idm(k)+1) = auxt3(k,m)
             end do
             c3 = c3+1
             interior_nodes_vef%l(c3) = Q_gijk(ijk_g,nd,p) ! Store the local numbering in ndxob%l
          end do
       end do

       ! All node numbering
       !Loop over the translations
       do l = 1,2**(nd-od)
          !ijk_g(jdm) will contain the translations from 1 vef to another; must be scaled by p
          !ijk_g(idm) will contain the variations inside the vef
          do k = 1,nd-od
             ijk_g(fdm(k)+1) = auxt1(k,l)*p
          end do

          ! Loop over the interior nodes of the vef
          do m = 1,(p+1)**(od)
             do k = 1,od
                ijk_g(idm(k)+1) = auxt4(k,m)
             end do
             c4 = c4+1
             nodes_vef%l(c4) = Q_gijk(ijk_g,nd,p) ! Store the local numbering in ntxob%l
          end do
       end do

       ! obxob array and auxiliar ob2node array 

       ! Interior node numbering
       ! Loop over the translations
       do l = 1,2**(nd-od)
          !ijk_g(jdm) will contain the translations from 1 vef to another; must be scaled by p
          !ijk_g(idm) will contain the variations inside the vef
          do k = 1,nd-od
             ijk_g(fdm(k)+1) = auxt1(k,l)*2
          end do

          ! Loop over the interior nodes of the vef
          do m = 1,(2-1)**(od)
             do k = 1,od
                ijk_g(idm(k)+1) = auxt6(k,m)
             end do
             c6 = c6+1
             ob2node(c6) = Q_gijk(ijk_g,nd,2) ! Store the local numbering in ndxob%l
          end do
       end do


       do l = 1,2**(nd-od)

          ! Fixed values for the vef
          do k = 1,nd-od
             ijk_g(fdm(k)+1) = auxt1(k,l)*2
          end do

          ! Fill obxob (equivalent to ntxob for p=2)
          do m = 1,3**od
             do k = 1,od
                ijk_g(idm(k)+1) = auxt5(k,m)
             end do
             c5 = c5 +1
             vefs_vef%l(c5) = Q_gijk(ijk_g,nd,2)
          end do

          ! Define ijk_g for the node in the center of the vef
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

 ! Modify the identifiers of the nodes by the ids of the vef in obxob
 do c5 = 1, nt2
    vefs_vef%l(c5) = node2ob(vefs_vef%l(c5))
 end do

 ! Sort the array 
 !do co = 1, no
 !   call sort(obxob%p(co+1)-obxob%p(co),obxob%l(obxob%p(co):obxob%p(co+1)))
 !end do

 ! Deallocate OBDLA
 call memfree(obdla,__FILE__,__LINE__)
 call memfree(ob2node,__FILE__,__LINE__)
 call memfree(node2ob,__FILE__,__LINE__)

 ! ! Create the face permutation of nodes
 ! if (nd>2) then call memalloc(2*2**2,nodes_vef(3),o2n,__FILE__,__LINE__)

 ! write(*,*) 'orientation vefs'
 ! do od = 1,nd
 !    write(*,*) 'dime', od, '--------------------------'
 !    write(*,*) orientation(nvef_dim(od):nvef_dim(od+1)-1)
 ! end do
 ! write(*,*) 'no+1', no+1, 'ndxob%p'
 ! do od = 1,no+1
 !    write(*,*) ndxob%p(od), ', &'
 ! end do
 ! write(*,*) 'nn', nn, 'ndxob%l'
 ! do od = 1,nn
 !    write(*,*) ndxob%l(od), ', &'
 ! end do

 ! write(*,*) 'no+1', no+1, 'ntxob%p'
 ! do od = 1,no+1
 !    write(*,*) ntxob%p(od), ', &'
 ! end do
 ! write(*,*) 'nt', nt, 'ntxob%l'
 ! do od = 1,nt
 !    write(*,*) ntxob%l(od), ', &'
 ! end do

 ! write(*,*) 'no+1', no+1, 'crxob%p'
 ! do od = 1,no+1
 !    write(*,*) crxob%p(od), ', &'
 ! end do
 ! write(*,*) 'nc', nc, 'crxob%l'
 ! do od = 1,nc
 !    write(*,*) crxob%l(od), ', &'
 ! end do
end subroutine fill





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
 integer(ip)                :: dp,i,ipp

 do ipp = p0,p1
    if (d>nd) exit
    ijk(nd-d+1) = ipp 
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
subroutine Q_orientation_vef(o,fdm,od,nd,l)
 implicit none
 ! Parameters
 integer(ip), intent(out) :: o
 integer(ip), intent(in)  :: fdm(:)  ! fdm gives the orthogonal directions to the vef
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
end subroutine Q_orientation_vef

!==================================================================================================
! R_DIM construct OBDLA matrix
! For each vef, up to a displacement, we have an identifier id={1..nod}
! obdla(id,1) = dimension of the vef
! obdla(id,2:obdla(id,1)+1) = gives the directions that define the vef
recursive subroutine r_dim(co,idm,ko,i,nd,od,obdla,nod)
 implicit none
 integer(ip), intent(in)    :: ko !space position we begin to count from (ij=>i<j)
 integer(ip), intent(in)    :: i  !i=local space position we are currently modifying (i=0..od)
 integer(ip), intent(in)    :: nd,od,nod
 integer(ip), intent(inout) :: obdla(nod,nd+1)
 integer(ip), intent(inout) :: idm(od) 
 integer(ip), intent(inout) :: co      ! Pointer to the position of the vef
 integer(ip)                :: aux,ijk_c(nd),j,cd,kn,s

 !Given dimension od of the vef
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

end module quad_lagrangian_reference_fe_names


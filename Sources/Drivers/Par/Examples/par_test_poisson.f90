program par_test_poisson
  use serial_names
  use par_test_poisson_driver_names
  implicit none
  type(par_test_poisson_fe_driver_t)     :: test_driver  
  call fempar_init()  
  call test_driver%run_simulation()
  call fempar_finalize()
  !! Our data
  !type(par_context_t)                                     :: w_context
  !type(par_environment_t)                                 :: par_env
  !type(par_mesh_t)                                        :: par_mesh
  !type(par_conditions_t)                                  :: par_conditions
  !type(par_triangulation_t)                               :: par_triangulation
  !type(par_fe_space_t)                                    :: par_fe_space
  !type(p_reference_fe_t)                                  :: reference_fe_array_one(1)
  !type(fe_affine_operator_t)                              :: fe_affine_operator
  !type(poisson_discrete_integration_t)                    :: poisson_integration
  !type(iterative_linear_solver_t)                         :: iterative_linear_solver
  !type(fe_function_t)                                     :: fe_function
  !class(vector_t), pointer                                :: dof_values
  !type(mlbddc_t)                                          :: mlbddc
  
  !integer(ip)              :: num_levels
  !integer(ip), allocatable :: parts_mapping(:), num_parts_per_level(:)
  !integer(ip), allocatable :: coarsening_ratios(:,:)
  !integer(ip)              :: npsoc(3), npdir(3)
  
  !type(par_test_reference_fe_parameters_t) :: test_params
  
  !! Uniform mesh generation related data
  !type(uniform_mesh_descriptor_t)       :: gdata
  !type(uniform_conditions_descriptor_t) :: bdata
  
  !call fempar_init()

  !! Start parallel execution
  !call w_context%create()

  !call test_params%set_default_params()
  !call test_params%set_cli()
  !call test_params%parse()

  ! !num_levels = 3
  ! !call memalloc(num_levels, parts_mapping , __FILE__, __LINE__)
  ! !call memalloc(num_levels, num_parts_per_level, __FILE__, __LINE__)
  ! 
  ! !num_parts_per_level = (/test_params%nparts, 2, 1/)
  ! !if ( w_context%get_rank() < test_params%nparts ) then
  ! !  parts_mapping       = (/w_context%get_rank()+1, w_context%get_rank()/2+1, 1/)
  ! !else if ( w_context%get_rank() >= test_params%nparts) then
  ! !  parts_mapping       = (/w_context%get_rank()+1, w_context%get_rank()+1-test_params%nparts, 1/)
  ! !end if
  
  ! num_levels = 2
  ! call memalloc(num_levels, parts_mapping , __FILE__, __LINE__)
  ! call memalloc(num_levels, num_parts_per_level, __FILE__, __LINE__)
  ! num_parts_per_level = (/test_params%nparts, 1/)
  ! parts_mapping       = (/w_context%get_rank()+1, 1/)
  !  
  !!num_levels = 3
  !!call memalloc(num_levels, parts_mapping , __FILE__, __LINE__)
  !!call memalloc(num_levels, num_parts_per_level, __FILE__, __LINE__)
  !!call memalloc(3, num_levels, coarsening_ratios, __FILE__, __LINE__ )
  !!coarsening_ratios(1,1) = 1; coarsening_ratios(1,2) = 2 ; coarsening_ratios(1,3) = 2 ;
  !!coarsening_ratios(2,1) = 1; coarsening_ratios(2,2) = 2 ; coarsening_ratios(2,3) = 2 ;
  !!coarsening_ratios(3,1) = 1; coarsening_ratios(3,2) = 1 ; coarsening_ratios(3,3) = 1 ;
  !!npsoc(1) = 1; npsoc(2) = 1; npsoc(3)=1;
  !!npdir(1) = 4; npdir(2) = 4; npdir(3)=1;
  !!call mlevel_coarsening ( w_context%get_rank()+1, 2, & 
  !!                         num_levels, coarsening_ratios, & 
  !!                         npdir, npsoc, parts_mapping, num_parts_per_level)
  
  !! Create parallel environment
  !call par_env%create (w_context,&
  !                     num_levels,&
  !                     num_parts_per_level,&
  !                     parts_mapping )
  
  !!call par_env%print()
  
  !! Read mesh
  !!call par_mesh_read ( test_params%dir_path, test_params%prefix, par_env, par_mesh )
  !!Read boundary conditions
  !!call par_conditions_read(test_params%dir_path, test_params%prefix, par_mesh%f_mesh%npoin, par_env, par_conditions)
  !! Generate triangulation
  !!call par_mesh_to_triangulation (par_mesh, par_triangulation, par_conditions)
  
  !! 4x4 quadrilateral mesh distributed among 2x2 subdomains
  !call uniform_mesh_descriptor_create(gdata, &
  !                                    nex=8, ney=8, nez=0, &
  !                                    npx=2, npy=2, npz=0, &
  !                                    nsx=1, nsy=1, nsz=0, &
  !                                    lx=1.0, ly=1.0, lz=0.0 )
  
  !! Create uniform conditions descriptor for a  
  !! scalar problem (i.e., ncode=nvalu=1).
  !! The uniform conditions descriptor is created within
  !! this subroutine s.t. all DoFs on top of vefs at the
  !! boundary of the domain are subject to strong Dirichlet
  !! boundary conditions. Recall that the function to be
  !! imposed is not extracted from the one prescribed within
  !! bdata, but from the call to par_fe_space%update_bc_value
  !! below.
  !call uniform_conditions_descriptor_create(bdata, &
  !                                          ncode = 1,&
  !                                          nvalu = 1,&
  !                                          ndime=gdata%ndime) 
  
  !!! Actually generate type(par_triangulation_t) and 
  !!! type(par_conditions_t).
  !call par_generate_uniform_triangulation(par_env, &
  !                                        gdata, &
  !                                        bdata, &
  !                                        par_triangulation, &
  !                                        par_conditions)
  
  !! Simple case
  !reference_fe_array_one(1) =  make_reference_fe ( topology = topology_quad, &
  !                                                 fe_type = fe_type_lagrangian, &
  !                                                 number_dimensions = 2, &
  !                                                 order = 1, &
  !                                                 field_type = field_type_scalar, &
  !                                                 continuity = .true. )
  
  !! call reference_fe_array_one(1)%p%print() 
  !  
  !call par_fe_space%create( par_triangulation = par_triangulation, &
  !                          par_boundary_conditions = par_conditions, &
  !                          reference_fe_phy = reference_fe_array_one ) 

  !call par_fe_space%update_bc_value (scalar_function=constant_scalar_function_t(1.0_rp), &
  !                                   bc_code = 1, &
  !                                   fe_space_component = 1 )
  
  !call par_fe_space%fill_dof_info()
  !! call par_fe_space%print()
  !call par_fe_space%renumber_dofs_first_interior_then_interface()
  !! call par_fe_space%print()

  !  
  !call fe_affine_operator%create (sparse_matrix_storage_format=csr_format, &
  !                                diagonal_blocks_symmetric_storage=(/.true./), &
  !                                diagonal_blocks_symmetric=(/.true./), &
  !                                diagonal_blocks_sign=(/SPARSE_MATRIX_SIGN_POSITIVE_DEFINITE/), &
  !                                environment=par_env, &
  !                                fe_space=par_fe_space, &
  !                                discrete_integration=poisson_integration )

  !call fe_affine_operator%symbolic_setup()
  !call fe_affine_operator%numerical_setup()
  
  !call par_fe_space%create_global_fe_function(fe_function)
  
  !call mlbddc%create(fe_affine_operator)
  !call mlbddc%symbolic_setup()
  !call mlbddc%numerical_setup()
  
  !! Create iterative linear solver, set operators and solve linear system
  !call iterative_linear_solver%create(par_env)
  !call iterative_linear_solver%set_type_from_string(cg_name)
  !call iterative_linear_solver%set_operators(fe_affine_operator, mlbddc)
  !dof_values => fe_function%get_dof_values()
  !call iterative_linear_solver%solve(fe_affine_operator%get_translation(),dof_values)
  !call iterative_linear_solver%free()
  !call mlbddc%free()
  
  !!select type(dof_values)
  !! type is (par_scalar_array_t)
  !!   call dof_values%print(6)
  !!end select
  
  
  !call fe_function%free()
  !call fe_affine_operator%free()
  !call par_fe_space%free()
  !call reference_fe_array_one(1)%free()

  !call par_triangulation_free(par_triangulation)
  !call par_conditions_free (par_conditions)
  !!call par_mesh_free (par_mesh)
  
  !call uniform_conditions_descriptor_free(bdata)

  !call memfree(parts_mapping , __FILE__, __LINE__)
  !call memfree(num_parts_per_level, __FILE__, __LINE__)
  !!call memfree(coarsening_ratios, __FILE__, __LINE__)
  
  !call par_env%free()
  !call w_context%free(finalize=.true.)
  !call fempar_finalize()

contains

  !!================================================================================================
  !subroutine global_to_ijk(lpart,nsckt,npsoc,ndime,ijk)
  !  implicit none
  !  integer(ip), intent(in)  :: lpart,ndime,nsckt(ndime),npsoc(ndime)
  !  integer(ip), intent(out) :: ijk(3)

  !  integer(ip) :: isckt,jsckt,ksckt,lsckt
  !  integer(ip) :: ipart_aux,jpart_aux,kpart_aux,lpart_aux
  !  integer(ip) :: ipart,jpart,kpart
  !  real(rp) :: aux1,aux2,aux3,aux4,aux5,aux6,aux7

  !  if(ndime==2) then
  !     ! First level (Socket level)
  !     aux1 = (lpart-1)/(npsoc(1)*npsoc(2))
  !     lsckt = floor(aux1) + 1
  !     aux2 = (lsckt-1)/nsckt(2)
  !     jsckt = lsckt - floor(aux2)*nsckt(2)
  !     isckt = floor(aux2) + 1

  !     ! Second level (Part inside the socket)
  !     lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)
  !     aux3 = (lpart_aux-1)/npsoc(2)
  !     jpart_aux = lpart_aux - floor(aux3)*npsoc(2)
  !     ipart_aux = floor(aux3) + 1

  !     ! ijk numeration
  !     ipart = ipart_aux + (isckt-1)*npsoc(1)
  !     jpart = jpart_aux + (jsckt-1)*npsoc(2)
  !     kpart = 1
  !     ijk = (/ipart,jpart,kpart/)
  !  else if(ndime==3) then
  !     ! First level (Socket level)
  !     aux1 = (lpart-1)/(npsoc(1)*npsoc(2)*npsoc(3))
  !     lsckt = floor(aux1) + 1
  !     aux2 = (lsckt-1)/(nsckt(2)*nsckt(3))
  !     isckt = floor(aux2) + 1
  !     aux3 = (lsckt - (isckt-1)*nsckt(2)*nsckt(3) - 1)/nsckt(3)
  !     jsckt = floor(aux3) + 1
  !     aux4 = lsckt - (isckt-1)*nsckt(2)*nsckt(3) - (jsckt-1)*nsckt(3) - 1
  !     ksckt = floor(aux4) + 1

  !     ! Second level (Part inside the socket)
  !     lpart_aux = lpart - (lsckt-1)*npsoc(1)*npsoc(2)*npsoc(3)
  !     aux5 = (lpart_aux-1)/(npsoc(2)*npsoc(3))
  !     ipart_aux = floor(aux5) + 1
  !     aux6 = (lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - 1)/npsoc(3)
  !     jpart_aux = floor(aux6) + 1
  !     aux7 = lpart_aux - (ipart_aux-1)*npsoc(2)*npsoc(3) - (jpart_aux-1)*npsoc(3) - 1
  !     kpart_aux = floor(aux7) + 1

  !     ! ijk numeration
  !     ipart = ipart_aux + (isckt-1)*npsoc(1)
  !     jpart = jpart_aux + (jsckt-1)*npsoc(2)
  !     kpart = kpart_aux + (ksckt-1)*npsoc(3)
  !     ijk = (/ipart,jpart,kpart/)
  !  end if

  !end subroutine global_to_ijk

  !! AFM: We should certainly move the following three subroutines:
  !!       * mlevel_coarsening; 
  !!       * mlevel_ijk; 
  !!       * gijk 
  !!     within Fempar (Any ideas with respect to which module 
  !!     could they fit?)
  !! **** VERY IMPORTANT NOTE** : These subroutines only work for the very 
  !! particular case of NPSOCX==NPSOCY==NPSOCZ==1. If this is not true, then
  !! the coarse-grid tasks may get coarse elements which are not contiguous in
  !! the coarse-grid mesh, rendering the recursive application of the BDDC possibly 
  !! leading to solvability issues in the local/global problems.
  !recursive subroutine mlevel_coarsening ( ipart, ndime, nlev, coarse, npdir, nsckt, id_parts, num_parts )
  !  implicit none
  !  ! Parameters
  !  integer(ip), intent(in)  :: ipart, ndime, nlev
  !  integer(ip), intent(in)  :: coarse(3,nlev)
  !  integer(ip), intent(in)  :: npdir(3), nsckt(3)
  !  integer(ip), intent(inout) :: id_parts(nlev)
  !  integer(ip), intent(inout) :: num_parts(nlev)

  !  ! Locals
  !  integer(ip) :: ijkpart(3)
  !  integer(ip) :: ijkcoarse(3,nlev)
  !  integer(ip) :: npdirmlev(3,nlev)
  !  integer(ip) :: i
  !  
  !  assert(ndime==1 .or. ndime==2 .or. ndime==3)
  ! 
  !  do i=1,ndime
  !     assert(nsckt(i)==1)
  !  end do

  !  call global_to_ijk(ipart,nsckt,npdir,ndime,ijkpart) 
  !  
  !  ! write(*,*) coarse
  !  ! write(*,*) ijkpart
  !  
  !  ijkcoarse(:,1) = ijkpart
  !  npdirmlev(:,1) = npdir
  !  call mlevel_ijk(ijkcoarse, ndime, nlev, coarse)
  !  call mlevel_ijk(npdirmlev, ndime, nlev, coarse)

  !  do i=1, nlev
  !     id_parts(i)=ijk_to_global(ijkcoarse(:,i),ndime,npdirmlev(:,i))
  !     num_parts(i)=ijk_to_global(npdirmlev(:,i),ndime,npdirmlev(:,i))
  !  end do

  !  ! Only higher level MPI tasks perform the recursive call
  !  ! as long as the number of levels is larger than two
  !  if ( ipart > num_parts(1) .and. nlev > 2 ) then
  !     call mlevel_coarsening ( ipart-num_parts(1), ndime, nlev-1, coarse(1,2), npdirmlev(1,2), nsckt, id_parts(2), num_parts(2) )
  !  else if ( nlev == 2 ) then
  !     ! Only in the coarse-grid task we have to fix id_parts(2)
  !     ! (in the fine-grid tasks id_parts(2) was already properly 
  !     ! by the code above this if-elseif block
  !     if ( ipart > num_parts(1) ) then
  !        id_parts(2) = 1
  !     endif
  !  end if
  !  
  !end subroutine mlevel_coarsening

  !subroutine mlevel_ijk(ijkc,d,l,coarse)
  !  implicit none
  !  integer(ip), intent(in) :: d,l,coarse(3,l)
  !  integer(ip), intent(inout) :: ijkc(3,l)
  !  integer(ip) :: i,lc
  !  do lc = 2,l
  !     do i=1,d
  !        ijkc(i,lc) = int((ijkc(i,lc-1)-1)/coarse(i,lc))+1
  !     end do
  !  end do    
  !end subroutine mlevel_ijk
  
  !integer(ip) function ijk_to_global(i,d,n)
  !  implicit none
  !  integer(ip) :: d,i(3),n(3),c,k
  !  ijk_to_global = 1
  !  c=1
  !  do k = d,1,-1
  !     ijk_to_global = ijk_to_global + (i(k)-1)*c
  !     c = c*n(k)
  !  end do
  !end function  ijk_to_global
  
end program par_test_poisson

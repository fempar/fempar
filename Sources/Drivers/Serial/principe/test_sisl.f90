program test_sisl
  use types_names
  use iso_c_binding
  use sisl_names
  implicit none
  type(c_ptr) :: curve,wcurve
  integer(ip), parameter :: ndime=3
  integer(ip), parameter :: n=5
  integer(ip), parameter :: p=2
  integer(ip), parameter :: num_knots=n+p+1
  real(rp) :: control_points( (ndime+1)*n)
  real(rp) :: knots(num_knots)
  real(rp) :: tol,length, parvalue

  integer(ip) :: p_shape(1)
  real(rp), pointer :: param(:)
  type(c_ptr)       :: p_param

  real(rp) :: point(ndime),derive(2*ndime)
  integer(ip) :: stat,i,j, num_int, num_curves, leftknot

  ! This is an arc taken from GiD (bounding half the base of the cylinder)
  control_points(1)=0.0_rp
  control_points(2)=10.0_rp
  control_points(3)=0.0_rp
  control_points(4)=1.0_rp

  control_points(5)=-10.0_rp*0.70710678118654757_rp
  control_points(6)=10.0_rp*0.70710678118654757_rp
  control_points(7)=0.0_rp
  control_points(8)=0.70710678118654757_rp
  
  control_points(9)=-10.0_rp
  control_points(10)=0.0_rp
  control_points(11)=0.0_rp
  control_points(12)=1.0_rp
  
  control_points(13)=-10.0_rp*0.70710678118654757_rp
  control_points(14)=-10.0_rp*0.70710678118654757_rp
  control_points(15)=0.0_rp
  control_points(16)=0.70710678118654757_rp
  
  control_points(17)= 0.0_rp
  control_points(18)=-10.0_rp
  control_points(19)=0.0_rp
  control_points(20)=1.0_rp

  knots( 1 )=0.0_rp
  knots( 2 )=0.0_rp
  knots( 3 )=0.0_rp
  knots( 4 )=0.5_rp
  knots( 5 )=0.5_rp
  knots( 6 )=1.0_rp
  knots( 7 )=1.0_rp
  knots( 8 )=1.0_rp

  write(*,*) 'Creating curve'
  curve = new_curve(n,p+1,knots,control_points,2,ndime,0)

  tol = 1.0e-8
  call curve_length(curve, tol , length, stat)

  write(*,*) 'Curve length status', stat
  write(*,*) 'Curve length', length

  point(1) = -10.0_rp
  point(2) = 0.0_rp
  point(3) = 0.0_rp
  call point_intersection(curve, point, ndime, tol, num_int, p_param, num_curves, wcurve, stat)
  p_shape(1) = num_int
  call c_f_pointer(p_param,param,p_shape)

  write(*,*) 'Point_intersection status', stat
  write(*,*) 'Point_intersection num_int', num_int
  write(*,*) 'Point_intersection param', param(1:num_int)
  write(*,*) 'Point_intersection num_curves', num_curves

  parvalue = 0.75_rp
  leftknot = 0
  call curve_left_evaluation(curve, 1, parvalue, leftknot, derive, stat) 
  write(*,*) 'Curve_left_evaluation status', stat
  write(*,*) 'Curve_left_evaluation leftknot', leftknot
  write(*,*) 'Curve_left_evaluation derive', derive

end program test_sisl

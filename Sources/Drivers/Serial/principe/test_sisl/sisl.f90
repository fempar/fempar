module sisl_names
  use types
  implicit none
  interface

     type(c_ptr) function new_curve(number, order, knots, coef, kind, dim, copy) & 
          & bind(c,NAME='newCurve')
       use iso_c_binding
       integer(c_int), value, intent(in) :: number, order, kind, dim, copy
       real(c_double)       , intent(in) :: knots(*),coef(*)
     end function new_curve

     subroutine point_intersection(pc1, pt1, idim, aepsge, jpt, gpar1, jcrv, wcurve, jstat) & 
          & bind(c,NAME='s1871')
       use iso_c_binding
       type(c_ptr)   , value, intent(in) :: pc1
       real(c_double)       , intent(in) :: pt1(idim)
       integer(c_int), value, intent(in) :: idim
       real(c_double), value, intent(in) :: aepsge
       integer(c_int)       , intent(out) :: jpt
       type(c_ptr)          , intent(out) :: gpar1
       integer(c_int)       , intent(out) :: jcrv
       type(c_ptr)          , intent(out) :: wcurve
       integer(c_int)       , intent(out) :: jstat
     end subroutine point_intersection
     ! void s1871(pc1, pt1, idim, aepsge, jpt, gpar1, jcrv, wcurve, jstat)
     ! SISLCurve *pc1;
     ! double *pt1;
     ! int idim;
     ! double aepsge;
     ! int *jpt;
     ! double **gpar1;
     ! int *jcrv;
     ! SISLIntcurve ***wcurve;
     ! int *jstat;

     subroutine curve_left_evaluation(curve, der, parvalue, leftknot, derive, stat) & 
          & bind(c,NAME='s1227')
       use iso_c_binding
       type(c_ptr)   , value, intent(in)    :: curve
       integer(c_int), value, intent(in)    :: der
       real(c_double), value, intent(in)    :: parvalue
       integer(c_int)       , intent(inout) :: leftknot
       real(c_double)       , intent(out)   :: derive(*)
       integer(c_int)       , intent(out)   :: stat
     end subroutine curve_left_evaluation
     ! void s1227(curve, der, parvalue, leftknot, derive, stat)
     ! SISLCurve *curve;
     ! int der;
     ! double parvalue;
     ! int *leftknot;
     ! double derive[ ];
     ! int *stat;

     subroutine curve_length(curve, epsge, length, stat) &
          & bind(c,NAME='s1240')
       use iso_c_binding
       type(c_ptr)   , value, intent(in) :: curve
       real(c_double), value, intent(in) :: epsge
       real(c_double), intent(out) :: length
       integer(c_int), intent(out) :: stat
     end subroutine curve_length
     ! void s1240(curve, epsge, length, stat)
     ! SISLCurve *curve;
     ! double epsge;
     ! double *length;
     ! int *stat;
  end interface

end module sisl_names

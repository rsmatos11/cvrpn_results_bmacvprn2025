!Fonte: https://people.math.sc.edu/Burkardt/f_src/geometry/geometry.f90

subroutine segments_int_2d ( p1, p2, q1, q2, flag, r )
  
  !*****************************************************************************80
  !
  !! SEGMENTS_INT_2D computes the intersection of two line segments in 2D.
  !
  !  Discussion:
  !
  !    A line segment is the finite portion of a line that lies between
  !    two points P1 and P2.
  !
  !    In 2D, two line segments might not intersect, even though the
  !    lines, of which they are portions, intersect.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    05 January 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the first
  !    segment.
  !
  !    Input, real ( kind = 8 ) Q1(2), Q2(2), the endpoints of the second
  !    segment.
  !
  !    Output, integer ( kind = 4 ) FLAG, records the results.
  !    0, the line segments do not intersect.
  !    1, the line segments intersect.
  !
  !    Output, real ( kind = 8 ) R(2), an intersection point, if there is one.
  !
  implicit none
  
  integer ( kind = 4 ), parameter :: dim_num = 2
  
  integer ( kind = 4 ) flag
  integer ( kind = 4 ) ival
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) q1(dim_num)
  real ( kind = 8 ) q2(dim_num)
  real ( kind = 8 ) r(dim_num)
  real ( kind = 8 ), parameter :: tol = 0.001D+00
  real ( kind = 8 ) u(dim_num)
  !
  !  Find the intersection of the two lines.
  !
  r(1:dim_num) = (/ 0.0D+00, 0.0D+00 /)
  
  call lines_exp_int_2d ( p1, p2, q1, q2, ival, r )
  
  if ( ival == 0 ) then
     flag = 0
     return
  end if
  
  if ( ival == 2 ) then
     if ( all( q1(1:dim_num) == p1(1:dim_num) ) .or. all( q1(1:dim_num) == p2(1:dim_num) ) ) then
        r(1:dim_num) = q1(1:dim_num)
        flag = 1
        return
     end if
     
     if ( all( q2(1:dim_num) == p1(1:dim_num) ) .or. all( q2(1:dim_num) == p2(1:dim_num) ) ) then
        r(1:dim_num) = q2(1:dim_num)
        flag = 1
        return
     end if
     !
     !  Is the q1 in the line segment given by p1 and p2?
     !
     call segment_contains_point_2d ( p1, p2, q1, u )
     if ( 0.0d0 .le. u(1) .and. u(1) .le. 1.0d+00 .and. u(2) .le. tol ) then
        r(1:dim_num) = q1(1:dim_num)
        flag = 1
        return
     end if
     !
     !  Is the q2 in the line segment given by p1 and p2?
     !
     call segment_contains_point_2d ( p1, p2, q2, u )
     if ( 0.0d0 .le. u(1) .and. u(1) .le. 1.0d+00 .and. u(2) .le. tol ) then
        r(1:dim_num) = q2(1:dim_num)
        flag = 1
        return
     end if
     !
     !  Is the p1 in the line segment given by q1 and q2?
     !
     call segment_contains_point_2d ( q1, q2, p1, u )
     if ( 0.0d0 .le. u(1) .and. u(1) .le. 1.0d+00 .and. u(2) .le. tol ) then
        r(1:dim_num) = p1(1:dim_num)
        flag = 1
        return
     end if
     !
     !  Is the p2 in the line segment given by q1 and q2?
     !
     call segment_contains_point_2d ( q1, q2, p2, u )
     if ( 0.0d0 .le. u(1) .and. u(1) .le. 1.0d+00 .and. u(2) .le. tol ) then
        r(1:dim_num) = p2(1:dim_num)
        flag = 1
        return
     end if
     
     flag = 0
     return
  end if  
  !
  !  Is the intersection point part of the first line segment?
  !
  call segment_contains_point_2d ( p1, p2, r, u )
  if ( u(1) < 0.0D+00 .or. 1.0D+00 < u(1) .or. tol < u(2) ) then
     flag = 0
     return
  end if
  !
  !  Is the intersection point part of the second line segment?
  !
  call segment_contains_point_2d ( q1, q2, r, u )
  if ( u(1) < 0.0D+00 .or. 1.0D+00 < u(1) .or. tol < u(2) ) then
     flag = 0
     return
  end if
  
  flag = 1
  
  return
end subroutine segments_int_2d

!================================================================================
!================================================================================

subroutine lines_exp_int_2d ( p1, p2, q1, q2, ival, p )
  
  !*****************************************************************************80
  !
  !! LINES_EXP_INT_2D determines where two explicit lines intersect in 2D.
  !
  !  Discussion:
  !
  !    The explicit form of a line in 2D is:
  !
  !      the line through the points P1 and P2.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    02 January 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the first line.
  !
  !    Input, real ( kind = 8 ) Q1(2), Q2(2), two points on the second line.
  !
  !    Output, integer ( kind = 4 ) IVAL, reports on the intersection:
  !    0, no intersection, the lines may be parallel or degenerate.
  !    1, one intersection point, returned in P.
  !    2, infinitely many intersections, the lines are identical.
  !
  !    Output, real ( kind = 8 ) P(2), if IVAl = 1, P is
  !    the intersection point.  Otherwise, P = 0.
  !
  implicit none
  
  integer ( kind = 4 ), parameter :: dim_num = 2
  
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) ival
  logical ( kind = 4 ) point_1
  logical ( kind = 4 ) point_2
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) q1(2)
  real ( kind = 8 ) q2(2)
  
  ival = 0
  p(1:dim_num) = 0.0D+00
  !
  !  Check whether either line is a point. Verifica se as retas são degeneradas.
  !
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then
     point_1 = .true.
  else
     point_1 = .false.
  end if
  
  if ( all ( q1(1:dim_num) == q2(1:dim_num) ) ) then
     point_2 = .true.
  else
     point_2 = .false.
  end if  
  !
  !  Convert the lines to ABC format.
  !
  if ( .not. point_1 ) then
     call line_exp2imp_2d ( p1, p2, a1, b1, c1 )
  end if
  
  if ( .not. point_2 ) then
     call line_exp2imp_2d ( q1, q2, a2, b2, c2 )
  end if
  !
  !  Search for intersection of the lines.
  !
  if ( point_1 .and. point_2 ) then
     if ( all ( p1(1:dim_num) == q1(1:dim_num) ) ) then
        ival = 1
        p(1:dim_num) = p1(1:dim_num)
     end if
  else if ( point_1 ) then
     if ( a2 * p1(1) + b2 * p1(2) == c2 ) then
        ival = 1
        p(1:dim_num) = p1(1:dim_num)
     end if
  else if ( point_2 ) then
     if ( a1 * q1(1) + b1 * q1(2) == c1 ) then
        ival = 1
        p(1:dim_num) = q1(1:dim_num)
     end if
  else
     call lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )
  end if
  
  return
end subroutine lines_exp_int_2d

!================================================================================
!================================================================================

subroutine segment_contains_point_2d ( p1, p2, p, u )
  
  !*****************************************************************************80
  !
  !! SEGMENT_CONTAINS_POINT_2D reports if a line segment contains a point in 2D.
  !
  !  Discussion:
  !
  !    A line segment is the finite portion of a line that lies between
  !    two points P1 and P2.
  !
  !    In exact arithmetic, point P is on the line segment between
  !    P1 and P2 if and only if 0 <= U <= 1 and V = 0.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    17 August 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
  !
  !    Input, real ( kind = 8 ) P(2), a point to be tested.
  !
  !    Output, real ( kind = 8 ) U(2), the components of P, with the first
  !    component measured along the axis with origin at P1 and unit at P2, 
  !    and second component the magnitude of the off-axis portion of the
  !    vector P-P1, measured in units of (P2-P1).
  !
  implicit none
  
  integer ( kind = 4 ), parameter :: dim_num = 2
  
  real ( kind = 8 ) normsq
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) u(dim_num)

  normsq = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )
  
  if ( normsq == 0.0D+00 ) then
     
     if ( all ( p(1:dim_num) == p1(1:dim_num) ) ) then
        u(1) = 0.5D+00
        u(2) = 0.0D+00
     else
        u(1) = 0.5D+00
        u(2) = huge ( u(2) )
     end if
     
  else
     
     u(1) = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
          * ( p2(1:dim_num) - p1(1:dim_num) ) ) / normsq
     
     u(2) = sqrt ( ( ( u(1) - 1.0D+00 ) * p1(1) - u(1) * p2(1) + p(1) )**2 &
          + ( ( u(1) - 1.0D+00 ) * p1(2) - u(1) * p2(2) + p(2) )**2 ) &
          / sqrt ( normsq )
     
  end if
  
  return
end subroutine segment_contains_point_2d

!================================================================================
!================================================================================

subroutine line_exp2imp_2d ( p1, p2, a, b, c )
  
  !*****************************************************************************80
  !
  !! LINE_EXP2IMP_2D converts an explicit line to implicit form in 2D.
  !
  !  Discussion:
  !
  !    The explicit form of a line in 2D is:
  !
  !      the line through the points P1 and P2.
  !
  !    The implicit form of a line in 2D is:
  !
  !      A * X + B * Y + C = 0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    06 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) P1(2), P2(2), two points on the line.
  !
  !    Output, real ( kind = 8 ) A, B, C, the implicit form of the line.
  !
  implicit none
  
  integer ( kind = 4 ), parameter :: dim_num = 2
  
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) c
  logical ( kind = 4 ) line_exp_is_degenerate_nd
  real ( kind = 8 ) norm
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  !
  !  Take care of degenerate cases.
  !
  if ( line_exp_is_degenerate_nd ( dim_num, p1, p2 ) ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'LINE_EXP2IMP_2D - Warning!'
     write ( *, '(a)' ) '  The line is degenerate.'
  end if
  
  a = p2(2) - p1(2)
  b = p1(1) - p2(1)
  c = p2(1) * p1(2) - p1(1) * p2(2)
  
  norm = a * a + b * b + c * c
  
  if ( 0.0D+00 < norm ) then
     a = a / norm
     b = b / norm
     c = c / norm
  end if
  
  if ( a < 0.0D+00 ) then
     a = -a
     b = -b
     c = -c
  end if
  
  return
end subroutine line_exp2imp_2d

!================================================================================
!================================================================================

subroutine lines_imp_int_2d ( a1, b1, c1, a2, b2, c2, ival, p )
  
  !*****************************************************************************80
  !
  !! LINES_IMP_INT_2D determines where two implicit lines intersect in 2D.
  !
  !  Discussion:
  !
  !    The implicit form of a line in 2D is:
  !
  !      A * X + B * Y + C = 0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    25 February 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) A1, B1, C1, define the first line.
  !    At least one of A1 and B1 must be nonzero.
  !
  !    Input, real ( kind = 8 ) A2, B2, C2, define the second line.
  !    At least one of A2 and B2 must be nonzero.
  !
  !    Output, integer ( kind = 4 ) IVAL, reports on the intersection.
  !
  !    -1, both A1 and B1 were zero.
  !    -2, both A2 and B2 were zero.
  !     0, no intersection, the lines are parallel.
  !     1, one intersection point, returned in P.
  !     2, infinitely many intersections, the lines are identical.
  !
  !    Output, real ( kind = 8 ) P(2), if IVAL = 1, then P is
  !    the intersection point.  Otherwise, P = 0.
  !
  implicit none
  
  integer ( kind = 4 ), parameter :: dim_num = 2
  
  real ( kind = 8 ) a(dim_num,dim_num+1)
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) c1
  real ( kind = 8 ) c2
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ival
  logical ( kind = 4 ) line_imp_is_degenerate_2d
  real ( kind = 8 ) p(dim_num)

  p(1:dim_num) = 0.0D+00
  !
  !  Refuse to handle degenerate lines.
  !
  if ( line_imp_is_degenerate_2d ( a1, b1 ) ) then
     ival = -1
     return
  end if
  
  if ( line_imp_is_degenerate_2d ( a2, b2 ) ) then
     ival = -2
     return
  end if
  !
  !  Set up and solve a linear system.
  !
  a(1,1) = a1
  a(1,2) = b1
  a(1,3) = -c1

  a(2,1) = a2
  a(2,2) = b2
  a(2,3) = -c2
  
  call r8mat_solve ( 2, 1, a, info )
  !
  !  If the inverse exists, then the lines intersect at the solution point.
  !
  if ( info == 0 ) then
     
     ival = 1
     p(1:dim_num) = a(1:dim_num,3)
     !
     !  If the inverse does not exist, then the lines are parallel
     !  or coincident.  Check for parallelism by seeing if the
     !  C entries are in the same ratio as the A or B entries.
     !
  else
     
     ival = 0
     
     if ( a1 == 0.0D+00 ) then
        if ( b2 * c1 - c2 * b1 < 1.0D-14 ) then
           ival = 2
        end if
     else
        if ( a2 * c1 - c2 * a1 < 1.0D-14 ) then
           ival = 2
        end if
     end if
  end if
  
  return
  
end subroutine lines_imp_int_2d

!================================================================================
!================================================================================

subroutine r8mat_solve ( n, rhs_num, a, info )
  
  !*****************************************************************************80
  !
  !! R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    29 August 2003
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input, integer ( kind = 4 ) RHS_NUM, the number of right hand sides.  
  !    RHS_NUM must be at least 0.
  !
  !    Input/output, real ( kind = 8 ) A(N,N+rhs_num), contains in rows and
  !    columns 1 to N the coefficient matrix, and in columns N+1 through
  !    N+rhs_num, the right hand sides.  On output, the coefficient matrix
  !    area has been destroyed, while the right hand sides have
  !    been overwritten with the corresponding solutions.
  !
  !    Output, integer ( kind = 4 ) INFO, singularity flag.
  !    0, the matrix was not singular, the solutions were computed;
  !    J, factorization failed on step J, and the solutions could not
  !    be computed.
  !
  implicit none
  
  integer ( kind = 4 ) n
  integer ( kind = 4 ) rhs_num
  
  real ( kind = 8 ) a(n,n+rhs_num)
  real ( kind = 8 ) apivot
  real ( kind = 8 ) factor
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) ipivot
  integer ( kind = 4 ) j
  
  info = 0

  do j = 1, n
     !
     !  Choose a pivot row.
     !
     ipivot = j
     apivot = a(j,j)
     
     do i = j+1, n
        if ( abs ( apivot ) < abs ( a(i,j) ) ) then
           apivot = a(i,j)
           ipivot = i
        end if
     end do
          
     !if ( apivot == 0.0D+00 ) then
     if ( abs(apivot) < 1.0D-14 ) then
        info = j
        return
     end if
     !
     !  Interchange.
     !
     do i = 1, n + rhs_num
        call r8_swap ( a(ipivot,i), a(j,i) )
     end do
     !
     !  A(J,J) becomes 1.
     !
     a(j,j) = 1.0D+00
     a(j,j+1:n+rhs_num) = a(j,j+1:n+rhs_num) / apivot
     !
     !  A(I,J) becomes 0.
     !
     do i = 1, n
        
        if ( i /= j ) then
           
           factor = a(i,j)
           a(i,j) = 0.0D+00
           a(i,j+1:n+rhs_num) = a(i,j+1:n+rhs_num) - factor * a(j,j+1:n+rhs_num)
           
        end if
        
     end do

  end do

  return
end subroutine r8mat_solve

!================================================================================
!================================================================================

function line_exp_is_degenerate_nd ( dim_num, p1, p2 )
  
  !*****************************************************************************80
  !
  !! LINE_EXP_IS_DEGENERATE_ND finds if an explicit line is degenerate in ND.
  !
  !  Discussion:
  !
  !    The explicit form of a line in ND is:
  !
  !      the line through the points P1 and P2.
  !
  !    An explicit line is degenerate if the two defining points are equal.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    06 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, integer ( kind = 4 ) DIM_NUM, the spatial dimension.
  !
  !    Input, real ( kind = 8 ) P1(DIM_NUM), P2(DIM_NUM), two points on the line.
  !
  !    Output, logical ( kind = 4 ) LINE_EXP_IS_DEGENERATE_ND, is TRUE if the line
  !    is degenerate.
  !
  implicit none
  
  integer ( kind = 4 ) dim_num

  logical ( kind = 4 ) line_exp_is_degenerate_nd
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)

  line_exp_is_degenerate_nd = ( all ( p1(1:dim_num) == p2(1:dim_num) ) )
  
  return
end function line_exp_is_degenerate_nd

!================================================================================
!================================================================================

function line_imp_is_degenerate_2d ( a, b )
  
  !*****************************************************************************80
  !
  !! LINE_IMP_IS_DEGENERATE_2D finds if an implicit point is degenerate in 2D.
  !
  !  Discussion:
  !
  !    The implicit form of a line in 2D is:
  !
  !      A * X + B * Y + C = 0
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    06 May 2005
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input, real ( kind = 8 ) A, B, C, the implicit line parameters.
  !
  !    Output, logical ( kind = 4 ) LINE_IMP_IS_DEGENERATE_2D, is true if the
  !    line is degenerate.
  !
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2
  
  real ( kind = 8 ) a
  real ( kind = 8 ) b
  !real ( kind = 8 ) c
  logical ( kind = 4 ) line_imp_is_degenerate_2d
  
  line_imp_is_degenerate_2d = ( a * a + b * b == 0.0D+00 )

  return
end function line_imp_is_degenerate_2d

!================================================================================
!================================================================================

subroutine r8_swap ( x, y )
  
  !*****************************************************************************80
  !
  !! R8_SWAP switches two R8's.
  !
  !  Licensing:
  !
  !    This code is distributed under the GNU LGPL license. 
  !
  !  Modified:
  !
  !    01 May 2000
  !
  !  Author:
  !
  !    John Burkardt
  !
  !  Parameters:
  !
  !    Input/output, real ( kind = 8 ) X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !
  implicit none
  
  real ( kind = 8 ) x
  real ( kind = 8 ) y
  real ( kind = 8 ) z
  
  z = x
  x = y
  y = z

  return
end subroutine r8_swap

!================================================================================
!================================================================================

subroutine polygon_point_near_2d ( n, v, p, pn, dist )

!*****************************************************************************80
!
!! POLYGON_POINT_NEAR_2D computes the nearest point on a polygon in 2D.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    28 February 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V(2,N), the polygon vertices.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest polygon point
!    is to be determined.
!
!    Output, real ( kind = 8 ) PN(2), the nearest point to P.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the
!    polygon.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) dist
  real ( kind = 8 ) dist2
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jp1
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) pn2(dim_num)
  real ( kind = 8 ) tval
  real ( kind = 8 ) v(dim_num,n)
!
!  Find the distance to each of the line segments that make up the edges
!  of the polygon.
!
  dist = huge ( dist )
  pn(1:dim_num) = 0.0D+00

  do j = 1, n

    jp1 = i4_wrap ( j+1, 1, n )

    call segment_point_near_2d ( v(1:dim_num,j), v(1:dim_num,jp1), p, &
      pn2, dist2, tval )

    if ( dist2 < dist ) then
      dist = dist2
      pn(1:dim_num) = pn2(1:dim_num)
    end if

  end do

  return
end subroutine polygon_point_near_2d

!================================================================================
!================================================================================

subroutine segment_point_near_2d ( p1, p2, p, pn, dist, t )

!*****************************************************************************80
!
!! SEGMENT_POINT_NEAR_2D: nearest point on line segment to point in 2D.
!
!  Discussion:
!
!    A line segment is the finite portion of a line that lies between
!    two points P1 and P2.
!
!    The nearest point will satisfy the condition
!
!      PN = (1-T) * P1 + T * P2.
!
!    T will always be between 0 and 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    03 May 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the endpoints of the line segment.
!
!    Input, real ( kind = 8 ) P(2), the point whose nearest neighbor
!    on the line segment is to be determined.
!
!    Output, real ( kind = 8 ) PN(2), the point on the line segment which is
!    nearest the point P.
!
!    Output, real ( kind = 8 ) DIST, the distance from the point to the 
!    nearest point on the line segment.
!
!    Output, real ( kind = 8 ) T, the relative position of the point PN
!    to the points P1 and P2.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  real ( kind = 8 ) bot
  real ( kind = 8 ) dist
  real ( kind = 8 ) p(dim_num)
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pn(dim_num)
  real ( kind = 8 ) t
!
!  If the line segment is actually a point, then the answer is easy.
!
  if ( all ( p1(1:dim_num) == p2(1:dim_num) ) ) then

    t = 0.0D+00

  else

    bot = sum ( ( p2(1:dim_num) - p1(1:dim_num) )**2 )

    t = sum ( ( p(1:dim_num)  - p1(1:dim_num) ) &
            * ( p2(1:dim_num) - p1(1:dim_num) ) ) / bot

    t = max ( t, 0.0D+00 )
    t = min ( t, 1.0D+00 )

  end if

  pn(1:dim_num) = p1(1:dim_num) + t * ( p2(1:dim_num) - p1(1:dim_num) )

  dist = sqrt ( sum ( ( p(1:dim_num) - pn(1:dim_num) )**2 ) )

  return
end subroutine segment_point_near_2d

!================================================================================
!================================================================================

function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    19 August 2003
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer
!    value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) jhi
  integer ( kind = 4 ) jlo
  integer ( kind = 4 ) wide

  jlo = min ( ilo, ihi )
  jhi = max ( ilo, ihi )

  wide = jhi - jlo + 1

  if ( wide == 1 ) then
    i4_wrap = jlo
  else
    i4_wrap = jlo + i4_modp ( ival - jlo, wide )
  end if

  return
end function i4_wrap

!================================================================================
!================================================================================

function i4_modp ( i, j )

!*****************************************************************************80
!
!! I4_MODP returns the nonnegative remainder of integer division.
!
!  Discussion:
!
!    If
!      NREM = I4_MODP ( I, J )
!      NMULT = ( I - NREM ) / J
!    then
!      I = J * NMULT + NREM
!    where NREM is always nonnegative.
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
!
!  Example:
!
!        I     J     MOD  I4_MODP    Factorization
!
!      107    50       7       7    107 =  2 *  50 + 7
!      107   -50       7       7    107 = -2 * -50 + 7
!     -107    50      -7      43   -107 = -3 *  50 + 43
!     -107   -50      -7      43   -107 =  3 * -50 + 43
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    02 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) I, the number to be divided.
!
!    Input, integer ( kind = 4 ) J, the number that divides I.
!
!    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
!    divided by J.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) j

  if ( j == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'I4_MODP - Fatal error!'
    write ( *, '(a,i8)' ) '  I4_MODP ( I, J ) called with J = ', j
    stop 1
  end if

  i4_modp = mod ( i, j )

  if ( i4_modp < 0 ) then
    i4_modp = i4_modp + abs ( j )
  end if

  return
end function i4_modp

!================================================================================
!================================================================================

subroutine polygon_contains_point_2d ( n, v, p, inside )

!*****************************************************************************80
!
!! POLYGON_CONTAINS_POINT_2D finds if a point is inside a polygon.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 November 2016
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of nodes or vertices in 
!    the polygon.  N must be at least 3.
!
!    Input, real ( kind = 8 ) V(2,N), the vertices of the polygon.
!
!    Input, real ( kind = 8 ) P(2), the coordinates of the point to be tested.
!
!    Output, logical ( kind = 4 ) INSIDE, is TRUE if the point is inside 
!    the polygon.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  logical ( kind = 4 ) inside
  !integer ( kind = 4 ) ip1
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) px1
  real ( kind = 8 ) px2
  real ( kind = 8 ) py1
  real ( kind = 8 ) py2
  real ( kind = 8 ) v(2,n)
  real ( kind = 8 ) xints

  inside = .false.

  px1 = v(1,1)
  py1 = v(2,1)
  xints = p(1) - 1.0D+00

  do i = 1, n

    px2 = v(1,mod(i,n)+1)
    py2 = v(2,mod(i,n)+1)

    if ( min ( py1, py2 ) < p(2) ) then
      if ( p(2) <= max ( py1, py2 ) ) then
        if ( p(1) <= max ( px1, px2 ) ) then
          if ( py1 /= py2 ) then
            xints = ( p(2) - py1 ) * ( px2 - px1 ) / ( py2 - py1 ) + px1
          end if
          if ( px1 == px2 .or. p(1) <= xints ) then
            inside = .not. inside
          end if
        end if
      end if
    end if

    px1 = px2
    py1 = py2

  end do

  return
end subroutine polygon_contains_point_2d

!================================================================================
!================================================================================

subroutine points_hull_2d ( node_num, node_xy, hull_num, hull )

!*****************************************************************************80
!
!! POINTS_HULL_2D computes the convex hull of 2D points.
!
!  Discussion:
!
!    The work involved is N*log(H), where N is the number of points, and H is
!    the number of points that are on the hull.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
!
!    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
!
!    Output, integer ( kind = 4 ) HULL_NUM, the number of nodes that lie on 
!    the convex hull.
!
!    Output, integer ( kind = 4 ) HULL(NODE_NUM).  Entries 1 through HULL_NUM 
!    contain the indices of the nodes that form the convex hull, in order.
!
  implicit none

  integer ( kind = 4 ) node_num

  real ( kind = 8 ) angle
  real ( kind = 8 ) angle_max
  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) di
  real ( kind = 8 ) dr
  integer ( kind = 4 ) first
  integer ( kind = 4 ) hull(node_num)
  integer ( kind = 4 ) hull_num
  integer ( kind = 4 ) i
  real ( kind = 8 ) node_xy(2,node_num)
  real ( kind = 8 ) p_xy(2)
  integer ( kind = 4 ) q
  real ( kind = 8 ) q_xy(2)
  integer ( kind = 4 ) r
  real ( kind = 8 ) r_xy(2)

  if ( node_num < 1 ) then
    hull_num = 0
    return
  end if
!
!  If NODE_NUM = 1, the hull is the point.
!
  if ( node_num == 1 ) then
    hull_num = 1
    hull(1) = 1
    return
  end if
!
!  If NODE_NUM = 2, then the convex hull is either the two distinct points,
!  or possibly a single (repeated) point.
!
  if ( node_num == 2 ) then

    if ( node_xy(1,1) /= node_xy(1,2) .or. node_xy(2,1) /= node_xy(2,2) ) then
      hull_num = 2
      hull(1) = 1
      hull(2) = 2
    else
      hull_num = 1
      hull(1) = 1
    end if

    return

  end if
!
!  Find the leftmost point and call it "Q".
!  In case of ties, take the bottom-most.
!
  q = 1
  do i = 2, node_num
    if ( node_xy(1,i) < node_xy(1,q) .or. &
       ( node_xy(1,i) == node_xy(1,q) .and. node_xy(2,i) < node_xy(2,q) ) ) then
      q = i
    end if
  end do

  q_xy(1:2) = node_xy(1:2,q)
!
!  Remember the starting point, so we know when to stop!
!
  first = q
  hull_num = 1
  hull(1) = q
!
!  For the first point, make a dummy previous point, 1 unit south,
!  and call it "P".
!
  p_xy(1) = q_xy(1)
  p_xy(2) = q_xy(2) - 1.0D+00
!
!  Now, having old point P, and current point Q, find the new point R
!  so the angle PQR is maximal.
!
!  Watch out for the possibility that the two nodes are identical.
!
  do

    r = 0
    angle_max = 0.0D+00

    do i = 1, node_num

      if ( i /= q .and. &
           ( node_xy(1,i) /= q_xy(1) .or. node_xy(2,i) /= q_xy(2) ) ) then

        angle = angle_rad_2d ( p_xy, q_xy, node_xy(1:2,i) )

        if ( r == 0 .or. angle_max < angle ) then

          r = i
          r_xy(1:2) = node_xy(1:2,r)
          angle_max = angle
!
!  In case of ties, choose the nearer point.
!
        else if ( r /= 0 .and. angle == angle_max ) then

          di = ( node_xy(1,i) - q_xy(1) )**2 + ( node_xy(2,i) - q_xy(2) )**2
          dr = ( r_xy(1)      - q_xy(1) )**2 + ( r_xy(2)      - q_xy(2) )**2

          if ( di < dr ) then
            r = i
            r_xy(1:2) = node_xy(1:2,r)
            angle_max = angle
          end if

        end if

      end if

    end do
!
!  We are done when we have returned to the first point on the convex hull.
!
    if ( r == first ) then
      exit
    end if

    hull_num = hull_num + 1

    if ( node_num < hull_num ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POINTS_HULL_2D - Fatal error!'
      write ( *, '(a)' ) '  The algorithm has failed.'
      stop 1
    end if
!
!  Add point R to convex hull.
!
    hull(hull_num) = r
!
!  Set P := Q, Q := R, and prepare to search for next point R.
!
    q = r

    p_xy(1:2) = q_xy(1:2)
    q_xy(1:2) = r_xy(1:2)

  end do

  return
end subroutine points_hull_2d

!================================================================================
!================================================================================

function angle_rad_2d ( p1, p2, p3 )

!*****************************************************************************80
!
!! ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
!
!  Discussion:
!
!    Except for the zero angle case, it should be true that
!
!      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
!
!        P1
!        /
!       /    
!      /     
!     /  
!    P2--------->P3
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    15 January 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), P3(2), define the rays
!    P1 - P2 and P3 - P2 which define the angle.
!
!    Output, real ( kind = 8 ) ANGLE_RAD_2D, the angle swept out by the rays,
!    in radians.  0 <= ANGLE_RAD_2D < 2 * PI.  If either ray has zero
!    length, then ANGLE_RAD_2D is set to 0.
!
  implicit none

  real ( kind = 8 ) angle_rad_2d
  real ( kind = 8 ) p(2)
  real ( kind = 8 ) p1(2)
  real ( kind = 8 ) p2(2)
  real ( kind = 8 ) p3(2)
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00

  p(1) = ( p3(1) - p2(1) ) * ( p1(1) - p2(1) ) &
       + ( p3(2) - p2(2) ) * ( p1(2) - p2(2) )

  p(2) = ( p3(1) - p2(1) ) * ( p1(2) - p2(2) ) &
       - ( p3(2) - p2(2) ) * ( p1(1) - p2(1) )

  if ( all ( p(1:2) == 0.0D+00)  ) then
    angle_rad_2d = 0.0D+00
    return
  end if

  angle_rad_2d = atan2 ( p(2), p(1) )

  if ( angle_rad_2d < 0.0D+00 ) then
    angle_rad_2d = angle_rad_2d + 2.0D+00 * r8_pi
  end if

  return
end function angle_rad_2d

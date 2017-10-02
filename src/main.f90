program PME
   use precision
   use constants
   use routines
   use mesh_routines
   use quadrature
   use basis_functions
   !*********************************************************************************
   !**                                                                             **
   !**     Program to solve the 2-D porous media equation using moving finite      **
   !**     elements.  The porous media equation is given as:-                      **
   !**                                                                             **
   !**                             u_t + (u^m u_x)_x =0                            **
   !**     Variables :-                                                            **
   !**                                                                             **
   !**             u (solution) - Double precision array                           **
   !**             x (mesh position in x-direction)- Double precision array        **
   !**             x (mesh position in y-direction)- Double precision array        **
   !**             x_dot (mesh velocity in x-direction) - Double precision array   **
   !**             y_dot (mesh velocity in y-direction) - Double precision array   **
   !**                                                                             **
   !**             nodes - number of nodes in mesh. Integer                        **
   !**             m - power in equation. Double precision                         **
   !**             output_t - Time of solution output. Double precision            **
   !**                                                                             **
   !**             delta_t - Size of time-step. Double precision                   **
   !**             total_t - total time of numerical simulation. Double precision  **
   !**                                                                             **
   !**                                                                             **
   !*********************************************************************************

   implicit none
!------------------------------------------------------------------------------
   integer :: nodes, no_of_tris, N, M, max_tris, mesh, order, nbdy, reports
   integer, dimension(:, :), allocatable :: tri, con
   integer, dimension(:), allocatable :: bdy
   real(kind=DP), dimension(:), allocatable :: u, x, y, x_dot, y_dot, theta, jac, w
   real(kind=DP), dimension(:, :), allocatable :: mass, ints, coeffs
   real(kind=DP) :: total_t, delta_t, output_t, Q, t_init, mpower
   real(kind=DP) :: report_step, report_time
   real(kind=DP) :: rzero, lambda, tzero, e, gamman, gammad
   real(kind=DP) :: error, xpos, ypos
   ! Computational timing variables
   integer :: System_Time_Start, System_Time_Stop, System_Time_Rate
   real :: CPU_Time_Start, CPU_Time_Stop
   real(kind=DP), external :: exact_solution
   integer :: i, j, count
   character(LEN=32) :: meshfile
   logical :: writesol

   character(LEN=10) :: numbers
   integer :: test_number, hundreds, tens, units
   character(LEN=32) :: filename
!------------------------------------------------------------------------------
   numbers = "0123456789"; filename = "solution"

   ! Print simulator banner.
   write (6, *)
   write (6, *) '---------------------------------------------------------------------'
   write (6, *) '                      Porous Media Equation'
   write (6, *) '---------------------------------------------------------------------'
   write (6, *)

   ! Start timing procedure
   call system_clock(System_Time_Start, System_Time_Rate)
   call cpu_time(CPU_Time_Start)

   open (unit=10, file='variables.data', status='old', form='formatted')

   read (10, *) N
   read (10, *) M
   read (10, *) Q
   read (10, *) t_init
   read (10, *) mpower
   read (10, *) output_t
   read (10, *) mesh
   read (10, *) reports
   read (10, *) order
   read (10, *) meshfile

   allocate (w(1:order), coeffs(1:order, 1:2)); w = 0.0d0; coeffs = 0.0d0; 
   call triquad(w, coeffs, order)

   if (mesh == 1) then
      nodes = N*M + 1; no_of_tris = M*(2*N - 1)

      allocate (tri(1:no_of_tris, 1:3)); tri = 0
      allocate (u(1:nodes)); u = 0
      allocate (x(1:nodes)); x = 0
      allocate (y(1:nodes)); y = 0
      allocate (x_dot(1:nodes)); x_dot = 0
      allocate (y_dot(1:nodes)); y_dot = 0
      allocate (theta(1:nodes)); theta = 0
      allocate (jac(1:no_of_tris)); jac = 0
      allocate (ints(1:3, 1:3)); ints = 0

      call initial_conditions(u, x, y, tri, mpower, Q, t_init, nodes, no_of_tris, N, M, output_t, reports)
   elseif (mesh == 2) then

      open (unit=20, file=meshfile, status='old', form='formatted')
      read (20, *) nodes, no_of_tris, nbdy
      allocate (tri(1:no_of_tris, 1:3)); tri = 0
      allocate (jac(1:no_of_tris)); jac = 0
      allocate (u(1:nodes)); u = 0
      allocate (x(1:nodes)); x = 0
      allocate (y(1:nodes)); y = 0
      allocate (x_dot(1:nodes)); x_dot = 0
      allocate (y_dot(1:nodes)); y_dot = 0
      allocate (theta(1:nodes)); theta = 0
      allocate (ints(1:3, 1:3)); ints = 0
      allocate (bdy(1:nbdy)); bdy = 0
      do i = 1, nodes
         read (20, *) x(i), y(i)
      end do
      do i = 1, no_of_tris
         read (20, *) tri(i, 1), tri(i, 2), tri(i, 3)
      end do
      do i = 1, nbdy
         read (20, *) bdy(i)
      end do

      ! sort the boundary node indexes
      call bubble(bdy, nbdy)

      ! set-up the initial data given by the self-similar solution
      gamman = gamma((1.0d0/mpower) + 2.0d0)
      gammad = gamma((1.0d0/mpower) + 1.0d0)

      rzero = sqrt(Q*gamman/(SQRT(pi)*gammad))
      tzero = ((rzero**2)*mpower)/(4.0d0*(mpower + 1.0d0))
      lambda = (t_init/tzero)**(1/(2.0d0*(mpower + 1.0d0)))

      x = rzero*lambda*x
      y = rzero*lambda*y

      do i = 1, nodes
         e = 1.0d0 - ((SQRT((x(i)**2 + y(i)**2))/(rzero*lambda))**2)
         if (e < 0.0d0) e = 0
         u(i) = (1.0d0/lambda**2)*(e**(1.0d0/mpower))
      end do

      ! write the solution variables to file
      open (unit=60, file='variables.m')
      write (60, *) mpower
      write (60, *) rzero
      write (60, *) tzero
      write (60, *) t_init
      write (60, *) output_t
      write (60, *) reports
      close (60)
   endif

   call max_no_of_tris(tri, nodes, no_of_tris, max_tris)
   max_tris = max_tris + 1
   allocate (con(1:nodes, 1:max_tris + 1)); con = 0
   allocate (mass(1:nodes, 1:max_tris + 1)); mass = 0

   call connectivity(con, tri, nodes, no_of_tris, max_tris)

   call mass_setup(u, x, y, nodes, no_of_tris, max_tris, tri, con, mass, theta, ints, jac)

   total_t = t_init; 
   delta_t = 0.001d0
   report_step = (output_t - t_init)/reports
   report_time = t_init + report_step
   writesol = .false.

   count = 1; 
   test_number = count
   hundreds = test_number/100
   test_number = test_number - 100*hundreds
   tens = test_number/10
   test_number = test_number - 10*tens
   units = test_number

   open (unit=100, file=trim(filename)//numbers(hundreds + 1:hundreds + 1)//numbers(tens + 1:tens + 1)// &
      &numbers(units + 1:units + 1)//".m")
   do i = 1, nodes
      write (100, *) x(i), y(i), u(i)
   end do
   close (100)
   count = count + 1

   do while (total_t < output_t)

      call adaptive_timestep(u, mpower, delta_t, total_t, nodes, no_of_tris, jac)

      if (total_t + delta_t > report_time) then
         delta_t = report_time - total_t
         report_time = report_time + report_step
         writesol = .true.
      endif

      if (total_t + delta_t > output_t) then
         delta_t = output_t - total_t
      endif

      total_t = total_t + delta_t

      ! Write out the current time and time-step size
      write (6, '(2f20.8)') total_t, delta_t

      call mesh_movement(u, x, y, x_dot, y_dot, tri, con, mass, nodes, no_of_tris, max_tris, delta_t, mpower, jac, w, coeffs, order)

      call u_calc(u, x, y, nodes, no_of_tris, max_tris, tri, con, mass, theta, ints, jac, bdy, nbdy)

      if (writesol) then

         test_number = count
         hundreds = test_number/100
         test_number = test_number - 100*hundreds
         tens = test_number/10
         test_number = test_number - 10*tens
         units = test_number

         open (unit=100, file=trim(filename)//numbers(hundreds + 1:hundreds + 1)//numbers(tens + 1:tens + 1)// &
            &numbers(units + 1:units + 1)//".m")
         do i = 1, nodes
            write (100, *) x(i), y(i), u(i)
         end do
         close (100)
         count = count + 1
         writesol = .false.
      endif
   end do

   open (unit=50, file='triangles.m')
   do i = 1, no_of_tris
      write (50, *) tri(i, :)
   end do
   close (50)

   ! Stop the timing
   call system_clock(System_Time_Stop)
   call cpu_time(CPU_Time_Stop)

   write (*, '(/)')
   write (*, '(1x,a)') 'Timing report:'
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(1x,a)') '    Elapsed Time   CPU Time', &
      '        (s)           (s)'
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(1x,2e15.4)') dble(System_Time_Stop - System_Time_Start)/dble(System_Time_Rate)&
                        &, CPU_Time_Stop - CPU_Time_Start
   write (6, *) '---------------------------------------------------------------------'
   write (*, '(/)')

end program PME

function exact_solution(x, y, lambda, rzero, m) result(uexact)

   use precision

   implicit none
!------------------------------------------------------------------------------
   real(kind=DP), intent(IN) :: x, y, lambda, rzero, m
   real(kind=DP) :: uexact, r
!------------------------------------------------------------------------------
   r = sqrt(x**2 + y**2)

   if (r > rzero*lambda) then
      uexact = 0.0d0
   else
      uexact = (1.0d0/(lambda**2))*(1.0d0 - ((r/(rzero*lambda))**2))**(1.0d0/m)
   endif

   return

end function exact_solution

subroutine initial_conditions(u, x, y, tri, mpower, Q, t_init, nodes, no_of_tris, N, M, output_t, reports)

   use precision
   use constants
   use routines

   implicit none
!------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, no_of_tris, N, M, reports
   integer, intent(INOUT), dimension(1:no_of_tris, 1:3) :: tri
   real(kind=DP), intent(INOUT), dimension(1:nodes) :: u, x, y
   real(kind=DP), intent(IN) :: Q, t_init, mpower, output_t
!------------------------------------------------------------------------------
   real(kind=DP), dimension(1:N) :: r
   real(kind=DP), dimension(1:M + 1) :: angle
   real(kind=DP) :: gamman, gammad, rzero, tzero, lambda, delta_x, e
   integer :: i, j
!------------------------------------------------------------------------------
   i = 0; j = 0; angle = 0; r = 0; e = 0

   gamman = gamma((1.0d0/mpower) + 2.0d0)
   gammad = gamma((1.0d0/mpower) + 1.0d0)

   rzero = sqrt(Q*gamman/(SQRT(pi)*gammad))
   tzero = ((rzero**2)*mpower)/(4.0d0*(mpower + 1.0d0))
   lambda = (t_init/tzero)**(1/(2.0d0*(mpower + 1.0d0)))

   delta_x = rzero*lambda/N

   do i = 1, N
      r(i) = (delta_x*i)
   end do

   do i = 1, M + 1
      angle(i) = (2.0d0*pi*(i - 1))/M
   end do

   x(1) = 0; y(1) = 0; u(1) = 1.0d0/(lambda**2)

   do i = 1, M
      do j = 1, N
         x(j + 1 + (i - 1)*N) = r(j)*dcos(angle(i))
         y(j + 1 + (i - 1)*N) = r(j)*dsin(angle(i))
         e = 1.0d0 - ((SQRT((x(j + 1 + (i - 1)*N)**2 + y(j + 1 + (i - 1)*N)**2))/(rzero*lambda))**2)
         if (e < 0.0d0) e = 0
         u(j + 1 + (i - 1)*N) = (1.0d0/lambda**2)*(e**(1.0d0/mpower))
      end do
   end do

   do i = 1, M - 1
      tri(i, 1) = 1
      tri(i, 2) = 2 + (i - 1)*N
      tri(i, 3) = 2 + i*N
   end do

   tri(i, 1) = 1
   tri(i, 2) = 2 + (M - 1)*N
   tri(i, 3) = 2

   do i = 1, N - 1
      do j = 1, M, 2
         if (j < M) then
            tri(j + (j - 1) + (2*i - 1)*M, 1) = 1 + (j - 1)*N + i
            tri(j + (j - 1) + (2*i - 1)*M, 2) = 2 + (j - 1)*N + i
            tri(j + (j - 1) + (2*i - 1)*M, 3) = 2 + j*N + i
            tri(j + (j - 1) + 1 + (2*i - 1)*M, 1) = 1 + (j - 1)*N + i
            tri(j + (j - 1) + 1 + (2*i - 1)*M, 2) = 2 + j*N + i
            tri(j + (j - 1) + 1 + (2*i - 1)*M, 3) = 1 + j*N + i
            tri(j + (j - 1) + 2 + (2*i - 1)*M, 1) = 1 + j*N + i
            tri(j + (j - 1) + 2 + (2*i - 1)*M, 2) = 2 + j*N + i
            tri(j + (j - 1) + 2 + (2*i - 1)*M, 3) = 1 + (j + 1)*N + i
            tri(j + (j - 1) + 3 + (2*i - 1)*M, 1) = 1 + (j + 1)*N + i
            tri(j + (j - 1) + 3 + (2*i - 1)*M, 2) = 2 + j*N + i
            tri(j + (j - 1) + 3 + (2*i - 1)*M, 3) = 2 + (j + 1)*N + i
         endif
      end do
   end do

   do i = 1, no_of_tris
      if (tri(i, 1) > (1 + M*N)) then
         tri(i, 1) = tri(i, 1) - M*N
      endif
      if (tri(i, 2) > (1 + M*N)) then
         tri(i, 2) = tri(i, 2) - M*N
      endif
      if (tri(i, 3) > (1 + M*N)) then
         tri(i, 3) = tri(i, 3) - M*N
      endif
   end do

   open (unit=60, file='variables.m')
   write (60, *) mpower
   write (60, *) rzero
   write (60, *) tzero
   write (60, *) t_init
   write (60, *) output_t
   write (60, *) reports
   close (60)

   return

end subroutine initial_conditions

subroutine mass_setup(u, x, y, nodes, no_of_tris, max_tris, tri, con, mass, theta, ints, jac)

   use precision
   use mesh_routines

   implicit none
!------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, no_of_tris, max_tris
   integer, intent(IN), dimension(1:no_of_tris, 1:3) :: tri
   integer, intent(IN), dimension(1:nodes, 1:max_tris + 1) :: con
   real(kind=DP), intent(INOUT), dimension(1:no_of_tris) :: jac
   real(kind=DP), intent(IN), dimension(1:nodes) :: u, x, y
   real(kind=DP), intent(INOUT), dimension(1:nodes, 1:max_tris + 1) :: mass
   real(kind=DP), intent(INOUT), dimension(1:nodes) :: theta
   real(kind=DP), intent(INOUT), dimension(1:3, 1:3) :: ints
!------------------------------------------------------------------------------
   integer :: i, j, k, L, loc
!------------------------------------------------------------------------------
   i = 0; j = 0; k = 0; L = 0; mass = 0.0d0; theta = 0.0d0; ints = 0.0d0; jac = 0.0d0

   ints(1, 1) = 1.0d0/12.0d0; ints(1, 2) = 1.0d0/24.0d0; ints(1, 3) = 1.0d0/24.0d0
   ints(2, 1) = 1.0d0/24.0d0; ints(2, 2) = 1.0d0/12.0d0; ints(2, 3) = 1.0d0/24.0d0
   ints(3, 1) = 1.0d0/24.0d0; ints(3, 2) = 1.0d0/24.0d0; ints(3, 3) = 1.0d0/12.0d0

   do i = 1, no_of_tris
      jac(i) = jacobian(x(tri(i, 1)), x(tri(i, 2)), x(tri(i, 3)), y(tri(i, 1)), y(tri(i, 2)), y(tri(i, 3)))

      if (jac(i) < 0.0d0) print *, 'triangle not ordered anti-clockwise.'

      do j = 1, 3
         do k = 1, 3
            loc = 0
            do L = 1, max_tris + 1
               if (tri(i, k) == con(tri(i, j), L)) then
                  loc = L
                  exit
               endif
            end do
            mass(tri(i, j), loc) = mass(tri(i, j), loc) + jac(i)*ints(j, k)
         end do
      end do
   end do

   do i = 1, nodes
      do j = 1, max_tris + 1
         if (con(i, j) /= 0) then
            theta(i) = theta(i) + mass(i, j)*u(con(i, j))
         endif
      end do
   end do

   return

end subroutine mass_setup

subroutine mesh_movement(u, x, y, x_dot, y_dot, tri, con, mass, nodes, no_of_tris, max_tris, delta_t, mpower, jac, w, coeffs, order)

   use precision
   use matrix_solvers
   use mesh_routines
   use basis_functions

   implicit none
!------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, no_of_tris, max_tris, order
   integer, intent(IN), dimension(1:no_of_tris, 1:3) :: tri
   integer, intent(IN), dimension(1:nodes, 1:max_tris + 1) :: con
   real(kind=DP), intent(IN), dimension(1:nodes, 1:max_tris + 1) :: mass
   real(kind=DP), intent(IN), dimension(1:no_of_tris) :: jac
   real(kind=DP), intent(INOUT), dimension(1:nodes) :: x, y, x_dot, y_dot
   real(kind=DP), intent(IN), dimension(1:nodes) :: u
   real(kind=DP), intent(IN), dimension(1:order, 1:2) :: coeffs
   real(kind=DP), intent(IN), dimension(1:order) :: w
   real(kind=DP), intent(IN) :: delta_t, mpower
!------------------------------------------------------------------------------
   real(kind=DP), dimension(1:nodes, 1:max_tris + 1) :: B
   real(kind=DP), dimension(1:nodes) :: psi, B_1, B_2, rhs
   real(kind=DP), dimension(1:no_of_tris, 1:2, 1:3) :: grad
   real(kind=DP) :: integral_1, integral_2
   integer :: i, j, k, L, loc
!------------------------------------------------------------------------------
   B = 0; psi = 0; i = 0; j = 0; k = 0; L = 0; B_1 = 0; B_2 = 0; rhs = 0; grad = 0

   do i = 1, no_of_tris

      call gradient(grad(i, :, :), x(tri(i, 1)), x(tri(i, 2)), x(tri(i, 3)), &
                                 & y(tri(i, 1)), y(tri(i, 2)), y(tri(i, 3)), jac(i))

      integral_1 = jac(i)*((u(tri(i, 1)) + u(tri(i, 2)) + u(tri(i, 3)))/6.0d0)

      do j = 1, 3
         do k = 1, 3
            loc = 0
            do L = 1, max_tris + 1
               if (tri(i, k) == con(tri(i, j), L)) then
                  loc = L
                  exit
               endif
            end do
            B(tri(i, j), loc) = B(tri(i, j), loc) + integral_1*DOT_PRODUCT(grad(i, :, j), grad(i, :, k))
         end do
         integral_2 = 0.0d0
         do k = 1, order
            integral_2 = integral_2 + 0.5d0*jac(i)*w(k)*( &
               &(u(tri(i, 1))*phi_A(coeffs(k, 1), coeffs(k, 2)) + &
               & u(tri(i, 2))*phi_B(coeffs(k, 1), coeffs(k, 2)) + &
               & u(tri(i, 3))*phi_C(coeffs(k, 1), coeffs(k, 2)))**mpower)
         end do

         rhs(tri(i, j)) = rhs(tri(i, j)) - &
            &(u(tri(i, 1))*dot_product(grad(i, :, j), grad(i, :, 1)) + &
            & u(tri(i, 2))*dot_product(grad(i, :, j), grad(i, :, 2)) + &
            & u(tri(i, 3))*dot_product(grad(i, :, j), grad(i, :, 3)))*integral_2
      end do
   end do

   call sparse_cg(B, Con, psi, rhs, nodes, max_tris, 1)

   !will need to include velocity field V.

   do i = 1, no_of_tris
      do j = 1, 3
         B_1(tri(i, j)) = B_1(tri(i, j)) + jac(i)*(psi(tri(i, 1))*grad(i, 1, 1) +  &
            &psi(tri(i, 2))*grad(i, 1, 2) + psi(tri(i, 3))*grad(i, 1, 3))/6.0d0
         B_2(tri(i, j)) = B_2(tri(i, j)) + jac(i)*(psi(tri(i, 1))*grad(i, 2, 1) +  &
            &psi(tri(i, 2))*grad(i, 2, 2) + psi(tri(i, 3))*grad(i, 2, 3))/6.0d0
      end do
   end do

   call sparse_cg(mass, Con, x_dot, B_1, nodes, max_tris, 1)
   call sparse_cg(mass, Con, y_dot, B_2, nodes, max_tris, 1)

   !time_step

   do i = 1, nodes
      x(i) = x(i) + (delta_t*x_dot(i))
      y(i) = y(i) + (delta_t*y_dot(i))
   end do

   return

end subroutine mesh_movement

subroutine u_calc(u, x, y, nodes, no_of_tris, max_tris, tri, con, mass, theta, ints, jac, bdy, nbdy)

   use precision
   use mesh_routines
   use matrix_solvers

   implicit none
!------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, no_of_tris, max_tris, nbdy
   integer, intent(IN), dimension(1:no_of_tris, 1:3) :: tri
   integer, intent(IN), dimension(1:nodes, 1:max_tris + 1) :: con
   integer, intent(IN), dimension(1:nbdy) :: bdy
   real(kind=DP), intent(INOUT), dimension(1:no_of_tris) :: jac
   real(kind=DP), intent(INOUT), dimension(1:nodes) :: u
   real(kind=DP), intent(IN), dimension(1:nodes) :: x, y, theta
   real(kind=DP), intent(INOUT), dimension(1:nodes, 1:max_tris + 1) :: mass
   real(kind=DP), intent(IN), dimension(1:3, 1:3) :: ints
!------------------------------------------------------------------------------
   integer :: i, j, k, L, loc
!------------------------------------------------------------------------------
   i = 0; j = 0; k = 0; L = 0; mass = 0.0d0; jac = 0.0d0

   do i = 1, no_of_tris
      jac(i) = jacobian(x(tri(i, 1)), x(tri(i, 2)), x(tri(i, 3)), y(tri(i, 1)), y(tri(i, 2)), y(tri(i, 3)))
      do j = 1, 3
         do k = 1, 3
            loc = 0
            do L = 1, max_tris + 1
               if (tri(i, k) == con(tri(i, j), L)) then
                  loc = L
                  exit
               endif
            end do
            mass(tri(i, j), loc) = mass(tri(i, j), loc) + jac(i)*ints(j, k)
         end do
      end do
   end do

   call sparse_cg(mass, Con, u, theta, nodes, max_tris, 2)

   do i = 1, nbdy
      u(bdy(i)) = 0.0d0
   end do

   return

end subroutine u_calc

subroutine adaptive_timestep(u, m, delta_t, t, nodes, no_of_tris, jac)

   use precision
   use mesh_routines

   implicit none
!------------------------------------------------------------------------------
   integer, intent(IN) :: nodes, no_of_tris
   real(kind=DP), intent(IN), dimension(1:nodes) :: u
   real(kind=DP), intent(IN), dimension(1:no_of_tris) :: jac
   real(kind=DP), intent(IN) :: t, m
   real(kind=DP), intent(INOUT) :: delta_t
!------------------------------------------------------------------------------
   real(kind=DP) :: min_delta_x
   integer :: i
!------------------------------------------------------------------------------
   min_delta_x = 1000.0d0

   do i = 1, no_of_tris
      min_delta_x = dmin1((0.5d0*jac(i))**2, min_delta_x)
   end do

   delta_t = (50.0d0/(maxval(u)**m))*(t**(1.0d0/(2.0d0 + 2.0d0*m)))*min_delta_x
!   delta_t = (1.0d0/(maxval(u)**m))*min_delta_x
   return

end subroutine adaptive_timestep

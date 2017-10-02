module precision
!******************************************************************************
!**                                                                          **
!**     Module defining an integer kind DP or DDP.  DPP is accurate to 18    **
!**     significant digits over the exponent range of +-4932 and DP  is      **
!**     accurate to 15 significant digits over the exponent range of +-307.  **
!**                                                                          **
!******************************************************************************
   implicit none
!-------------------------------------------------------------------------------
   integer, parameter :: DP = selected_real_kind(15, 307)
   integer, parameter :: DDP = selected_real_kind(18, 4932)
!-------------------------------------------------------------------------------

end module precision

module constants
!******************************************************************************
!**                                                                          **
!**     Module defining a number of mathematical constants.  The constants   **
!**     are held to the precision DP which is defined in the module          **
!**     precision.                                                           **
!**                                                                          **
!******************************************************************************
   use precision

   implicit none
!-------------------------------------------------------------------------------
   real(DP), parameter :: pi = 3.141592653589793238462643_DP
!-------------------------------------------------------------------------------

end module constants

module routines

   implicit none

contains

   double precision function gamma(xx)

      implicit none
!---------------------------------------------------------------------------------
      double precision, intent(IN) :: xx
!---------------------------------------------------------------------------------
      double precision, dimension(1:6) :: cof
      double precision :: ser, stp, tmp, x, y
      integer :: j
!---------------------------------------------------------------------------------

      cof(1) = 76.18009172947146d0
      cof(2) = -86.50532032941677d0
      cof(3) = 24.01409824083091d0
      cof(4) = -1.231739572450155d0
      cof(5) = 0.1208650973866179d-2
      cof(6) = -0.5395239384953d-5
      stp = 2.5066282746310005d0

      x = xx; y = x; tmp = x + 5.50d0

      tmp = (x + 0.50d0)*dlog(tmp) - tmp

      ser = 1.000000000190015d0

      do j = 1, 6
         y = y + 1.0d0
         ser = ser + cof(j)/y
      end do

      gamma = tmp + dlog(stp*ser/x)
      gamma = dexp(gamma)

      return

   end function gamma

   subroutine bubble(vec, n)

      implicit none
      !==============================================================================
      integer, intent(IN) :: n
      integer, intent(INOUT), dimension(1:n) :: vec
      !==============================================================================
      integer :: i, j, swap, low, lowi
      !==============================================================================

      do i = 1, n - 1
         low = vec(i)
         lowi = i

         do j = i + 1, n
            if (vec(j) < low) then
               low = vec(j)
               lowi = j
            endif
         end do

         swap = vec(i)
         vec(i) = vec(lowi)
         vec(lowi) = swap

      end do

      return

   end subroutine bubble

end module routines

module matrix_solvers

   implicit none

contains

   subroutine sparse_cg(A, C, X, B, nodes, max_tris, pre)

      implicit none
      !==============================================================================
      integer, intent(IN) :: nodes, max_tris, pre
      integer, intent(IN), dimension(1:nodes, 1:max_tris + 1) :: C
      double precision, intent(IN), dimension(1:nodes, 1:max_tris + 1) :: A
      double precision, intent(INOUT), dimension(1:nodes) :: X
      double precision, intent(IN), dimension(1:nodes) :: B
      !==============================================================================
      double precision, dimension(1:nodes, 1:max_tris + 1) :: PA
      double precision, dimension(1:nodes) :: residual, direction, product_1, product_2, PB
      double precision :: alpha, beta, TOL, ans, denom
      integer :: i, j, count
      !------------------------------------------------------------------------------
      residual = 0; direction = 0; alpha = 0; beta = 0; TOL = 0.00000000001d0; 
      count = 0; product_1 = 0; product_2 = 0

      !Use value of X at previous time-step as the initial guess for the CG method.
      !precondition

      if (pre == 1) then
         PA = A; PB = B
      elseif (pre == 2) then
         PA = 0; PB = 0
         do i = 1, nodes
            !find diagonal elememt in A
            do j = 1, max_tris + 1
               if (C(i, j) == i) then
                  denom = A(i, j)
                  exit
               endif
            end do

            PB(i) = B(i)/denom
            PA(i, :) = A(i, :)/denom
         end do
      endif

      do i = 1, nodes
         ans = 0
         do j = 1, max_tris + 1
            if (c(i, j) /= 0) then
               ans = ans + PA(i, j)*X(C(i, j))
            endif
         end do

         residual(i) = PB(i) - ans
      end do

      do while (dabs(maxval(residual)) > TOL)
         if ((dabs(MAXVAL(direction)) < TOL) .AND. (dabs(MINVAL(direction)) < TOL)) then
            beta = 0.0d0
         else
            product_1 = 0; product_2 = 0
            do i = 1, nodes
               do j = 1, max_tris + 1
                  if (c(i, j) /= 0) then
                     product_1(i) = product_1(i) + PA(i, j)*direction(C(i, j))
                     product_2(i) = product_2(i) + PA(i, j)*residual(C(i, j))
                  endif
               end do
            end do
            beta = -((DOT_PRODUCT(direction, product_2))/(DOT_PRODUCT(direction, product_1)))
         endif

         direction = residual + beta*direction

         if ((dabs(MAXVAL(direction)) < TOL) .AND. (dabs(MINVAL(direction)) < TOL)) then
            alpha = 0.0d0
         else
            product_1 = 0

            do i = 1, nodes
               do j = 1, max_tris + 1
                  if (c(i, j) /= 0) then
                     product_1(i) = product_1(i) + PA(i, j)*direction(C(i, j))
                  endif
               end do
            end do

            alpha = ((DOT_PRODUCT(direction, residual))/(DOT_PRODUCT(direction, product_1)))
         endif

         X = X + (alpha*direction)

         do i = 1, nodes
            ans = 0
            do j = 1, max_tris + 1
               if (c(i, j) /= 0) then
                  ans = ans + PA(i, j)*X(C(i, j))
               endif
            end do
            residual(i) = PB(i) - ans
         end do

         count = count + 1

         if (count >= 10000) then
            print *, 'residual bounds', maxval(residual), minval(residual)
            exit
         endif
      end do

      return

   end subroutine sparse_cg

end module matrix_solvers

module mesh_routines

   implicit none

contains

   subroutine max_no_of_tris(tri, nodes, no_of_tris, max_tris)
      !******************************************************************************
      !**                                                                          **
      !**     Subroutine calculates the maximum number of triangles surrounding    **
      !**     any given node for a mesh.                                           **
      !**                                                                          **
      !**     GLOBAL VARIABLES :-                                                  **
      !**                     tri - integer array - contains data about triangles  **
      !**                     nodes - integer - number of nodes in mesh            **
      !**                     no_of_tris - integer - no. of triangles in mesh      **
      !**                     max_tris - integer - maximum no. of triangles        **
      !**                                          surrounding nodes               **
      !**     LOCAL VARIABLES :-                                                   **
      !**                     array - integer array - contains number of triangles **
      !**                                         surrounding each node.           **
      !**                     i, j - integer - loop counters                       **
      !**                                                                          **
      !******************************************************************************
      implicit none
      !------------------------------------------------------------------------------
      integer, intent(IN) :: nodes, no_of_tris
      integer, intent(IN), dimension(1:no_of_tris, 1:3) :: tri
      integer, intent(INOUT) :: max_tris
      !------------------------------------------------------------------------------
      integer, dimension(1:nodes) :: array
      integer :: i
      !------------------------------------------------------------------------------
      array = 0; i = 0

      do i = 1, no_of_tris
         array(tri(i, 1)) = array(tri(i, 1)) + 1
         array(tri(i, 2)) = array(tri(i, 2)) + 1
         array(tri(i, 3)) = array(tri(i, 3)) + 1
      end do

      max_tris = MAXVAL(array)

      return

   end subroutine max_no_of_tris

   subroutine connectivity(con, tri, nodes, no_of_tris, max_tris)
      !******************************************************************************
      !**                                                                          **
      !**     Subroutine calculates the connectivity of the mesh using the         **
      !**     subroutine search_tri.  The connectivity is stored in array con.     **
      !**     GLOBAL VARIABLES :-                                                  **
      !**                     con - integer array - connectivity of mesh           **
      !**                     nodes - integer - no. of nodes in mesh               **
      !**                     no_of_tris - integer - no. of triangles in mesh      **
      !**                     tri - integer array - contains data about triangles  **
      !**                     max_tris - integer - maximum no. of triangles        **
      !**                                          surrounding nodes               **
      !**     LOCAL VARIABLES :-                                                   **
      !**                     i - integer - loop counters                          **
      !**                                                                          **
      !******************************************************************************
      implicit none
      !------------------------------------------------------------------------------
      integer, intent(IN) :: nodes, no_of_tris, max_tris
      integer, dimension(1:no_of_tris, 1:3), intent(IN) :: tri
      integer, dimension(1:nodes, 1:max_tris + 1), intent(INOUT) :: con
      !------------------------------------------------------------------------------
      integer, dimension(1:nodes) :: counts
      logical :: location
      integer :: i, j, k, l
      !------------------------------------------------------------------------------
      counts = 1

      do i = 1, no_of_tris
         do L = 1, 3
            do j = 1, 3
               location = .false.
               do k = 1, max_tris + 1
                  if (tri(i, j) == con(tri(i, L), k)) then
                     location = .true.
                  endif
               end do

               if (location) then
               else
                  con(tri(i, L), counts(tri(i, L))) = tri(i, j)
                  counts(tri(i, L)) = counts(tri(i, L)) + 1
               endif
            end do
         end do
      end do

      return

   end subroutine connectivity

   function jacobian(x_1, x_2, x_3, y_1, y_2, y_3) result(jac)

      use precision

      implicit none
      !------------------------------------------------------------------------------
      real(kind=DP), intent(IN) :: x_1, x_2, x_3, y_1, y_2, y_3
      real(kind=DP) :: jac
      !------------------------------------------------------------------------------

      jac = ((x_2 - x_1)*(y_3 - y_1) - (x_3 - x_1)*(y_2 - y_1))

      return

   end function jacobian

end module mesh_routines

module basis_functions

   implicit none

contains

   subroutine gradient(x, x_1, x_2, x_3, y_1, y_2, y_3, jac)

      implicit none
      !------------------------------------------------------------------------------
      double precision, intent(INOUT), dimension(1:2, 1:3) :: x
      double precision, intent(IN) :: x_1, x_2, x_3, y_1, y_2, y_3, jac
      !------------------------------------------------------------------------------

      if (dabs(jac) .lt. 1d-10) then
         x(1, 1) = 0.0d0
         x(2, 1) = 0.0d0
         x(1, 2) = 0.0d0
         x(2, 2) = 0.0d0
         x(1, 3) = 0.0d0
         x(2, 3) = 0.0d0
      else
         x(1, 1) = -(y_3 - y_2)/jac
         x(2, 1) = (x_3 - x_2)/jac
         x(1, 2) = -(y_1 - y_3)/jac
         x(2, 2) = (x_1 - x_3)/jac
         x(1, 3) = -(y_2 - y_1)/jac
         x(2, 3) = (x_2 - x_1)/jac
      endif

      return

   end subroutine gradient

   function phi_A(eta, nu) result(phi)

      use precision

      implicit none
      !------------------------------------------------------------------------------
      real(kind=DP), intent(IN) :: eta, nu
      real(kind=DP) :: phi
      !------------------------------------------------------------------------------

      phi = 1.0-eta - nu

      return

   end function phi_A

   function phi_B(eta, nu) result(phi)

      use precision

      implicit none
      !------------------------------------------------------------------------------
      real(kind=DP), intent(IN) :: eta, nu
      real(kind=DP) :: phi
      !------------------------------------------------------------------------------

      phi = eta

      return

   end function phi_B

   function phi_C(eta, nu) result(phi)

      use precision

      implicit none
      !------------------------------------------------------------------------------
      real(kind=DP), intent(IN) :: eta, nu
      real(kind=DP) :: phi
      !------------------------------------------------------------------------------

      phi = nu

      return

   end function phi_c

end module basis_functions

module quadrature

   implicit none

contains

   subroutine triquad(w, coeffs, order)

      use precision

      implicit none
      !------------------------------------------------------------------------------
      integer, intent(IN) :: order
      real(kind=DP), intent(INOUT), dimension(1:order, 1:2) :: coeffs
      real(kind=DP), intent(INOUT), dimension(1:order) :: w
      !------------------------------------------------------------------------------

      if (order == 6) then
         ! integrates polynominals of degree 4 exactly.
         w(1) = 0.109951743655322d0; w(2) = 0.109951743655322d0
         w(3) = 0.109951743655322d0; w(4) = 0.223381589678011d0
         w(5) = 0.223381589678011d0; w(6) = 0.223381589678011d0
         coeffs(1, 1) = 0.816847572980459d0; coeffs(1, 2) = 0.091576213509771d0
         coeffs(2, 1) = 0.091576213509771d0; coeffs(2, 2) = 0.816847572980459d0
         coeffs(3, 1) = 0.091576213509771d0; coeffs(3, 2) = 0.091576213509771d0
         coeffs(4, 1) = 0.108103018168070d0; coeffs(4, 2) = 0.445948490915965d0
         coeffs(5, 1) = 0.445948490915965d0; coeffs(5, 2) = 0.108103018168070d0
         coeffs(6, 1) = 0.445948490915965d0; coeffs(6, 2) = 0.445948490915965d0
      elseif (order == 13) then
         ! integrates polynominals of degree 7 exactly.
         w(1) = -0.149570044467670d0; w(2) = 0.175615257433204d0
         w(3) = 0.175615257433204d0; w(4) = 0.175615257433204d0
         w(5) = 0.053347235608839d0; w(6) = 0.053347235608839d0
         w(7) = 0.053347235608839d0; w(8) = 0.077113760890257d0
         w(9) = 0.077113760890257d0; w(10) = 0.077113760890257d0
         w(11) = 0.077113760890257d0; w(12) = 0.077113760890257d0
         w(13) = 0.077113760890257d0; 
         coeffs(1, 1) = 0.333333333333333d0; coeffs(1, 2) = 0.333333333333333d0
         coeffs(2, 1) = 0.479308067841923d0; coeffs(2, 2) = 0.260345966079038d0
         coeffs(3, 1) = 0.260345966079038d0; coeffs(3, 2) = 0.479308067841923d0
         coeffs(4, 1) = 0.260345966079038d0; coeffs(4, 2) = 0.260345966079038d0
         coeffs(5, 1) = 0.869739794195568d0; coeffs(5, 2) = 0.065130102902216d0
         coeffs(6, 1) = 0.065130102902216d0; coeffs(6, 2) = 0.869739794195568d0
         coeffs(7, 1) = 0.065130102902216d0; coeffs(7, 2) = 0.065130102902216d0
         coeffs(8, 1) = 0.638444188569809d0; coeffs(8, 2) = 0.312865496004875d0
         coeffs(9, 1) = 0.638444188569809d0; coeffs(9, 2) = 0.048690315425316d0
         coeffs(10, 1) = 0.312865496004875d0; coeffs(10, 2) = 0.638444188569809d0
         coeffs(11, 1) = 0.312865496004875d0; coeffs(11, 2) = 0.048690315425316d0
         coeffs(12, 1) = 0.048690315425316d0; coeffs(12, 2) = 0.312865496004875d0
         coeffs(13, 1) = 0.048690315425316d0; coeffs(13, 2) = 0.638444188569809d0
      endif

      return

   end subroutine triquad

end module quadrature


module routines

   implicit none

contains

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


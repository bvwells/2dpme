module linear_solvers

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

end module linear_solvers

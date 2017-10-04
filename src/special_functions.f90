module special_functions

contains

   double precision function gamma(xx)
      !*********************************************************************************
      !**                                                                             **
      !**  This function calculates the gamma funtion for the value xx.               **
      !**                                                                             **
      !*********************************************************************************
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

end module special_functions

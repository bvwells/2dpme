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

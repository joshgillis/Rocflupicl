!-----------------------------------------------------------------------
!
! Created Feb. 1, 2024
!
! Subroutine to set user-defined values at time t=0
!
!-----------------------------------------------------------------------
!
      subroutine ppiclf_user_InitZero
!
      implicit none
!
      include "PPICLF"
!
! Internal:
!
      integer*4 i,j,k

!
! Code:
!
      ppiclf_TimeBH = 0.0d0

      ppiclf_drudtMixt = 0.0d0
      ppiclf_drudtPlag = 0.0d0


      return
      end


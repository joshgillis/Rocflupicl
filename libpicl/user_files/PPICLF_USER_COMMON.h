#include "PPICLF_USER.h"
#include "PPICLF_STD.h"


!
! General useage
!
      integer*4 :: stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag
      real*8 :: rmu_ref, tref, suth, ksp, erest
      common /RFLU_ppiclF/ stationary, qs_flag, am_flag, pg_flag,
     >   collisional_flag, heattransfer_flag, feedback_flag,
     >   qs_fluct_flag, ppiclf_debug, rmu_flag, rmu_ref, tref, suth,
     >   rmu_fixed_param, rmu_suth_param, qs_fluct_filter_flag,
     >   qs_fluct_filter_adapt_flag, ksp, erest,
     >   ViscousUnsteady_flag, ppiclf_nUnsteadyData,ppiclf_nTimeBH,
     >   sbNearest_flag, burnrate_flag, flow_model, pseudoTurb_flag

      real*8 rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,
     >   rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,
     >   vx,vy,vz,
     >   rcp_part,rpr,
     >   phi, mp, re, rem
      common /RFLU_user/ rpi,rmu,rkappa,rmass,vmag,rhof,dp,rep,rphip,
     >   rphif,asndf,rmachp,rhop,rhoMixt,reyL,rnu,fac,
     >   vx,vy,vz,
     >   rcp_part,rpr,
     >   phi,mp,re,rem

!
! For misc values
!
      real*8 OneThird
      common /ppiclf_misc01/ OneThird
 

!
! For ppiclf_user_Fluctuations.f
!
      integer*4 icpmean
      real*8 upmean, vpmean, wpmean, phipmean
      real*8 u2pmean, v2pmean, w2pmean
      common /user_fluct01/ icpmean
      common /user_fluct02/ upmean, vpmean, wpmean, phipmean
      common /user_fluct03/ u2pmean, v2pmean, w2pmean

      real*8 UnifRnd(6), Rsg(3,3)
      common /user_fluct02/ UnifRnd, Rsg

      real*8 k_tilde, k_Mach, b_tilde,
     >	     b_Mach, b_par, b_perp,
     >       Rmean_par, Rmean_perp
      common /user_fluct03/ k_tilde, k_Mach,
     >                      b_Mach, b_par, b_perp,
     >                      Rmean_par, Rmean_perp

      real*8 C1P, C2P, C3P, C4P, C5P,
     >       D1P, D2P, D3P, D4P, D5P, D6P, D7P, D8P,
     >       E1P, E2P, E3P, E4P
      real*8 F1P, F2P, F3P, F4P, F5P, F6P, F7p, F8P, F9P, F10P, F11P,
     >       G1P, G2P, G3P, G4P, G5P, G6P, G7P, G8P,
     >       A1P, A2P, A3P

      ! Constants taken from Osnes PseudoTurbulent paper, Table 1
      parameter C1P = -1.2152, C2P = -7.6314, C3P=0.2889, C4P=0.6143,
     >          C5P = 0.3082, 
     >          D1P = -0.0462, D2P = -0.1068, D3P=0.6793, D4P=1.1461,
     >          D5P = -2.6886, D6P = -2.1376, D7P=0.4873, D8P=0.2395,
     >          E1P = -0.2906, E2P=1.1899, E3P=0.5218, E4P=0.0699

      ! Constants taken from Osnes PseudoTurbulent paper, Table 2
      parameter F1P= -0.0022, F2P= -0.0219, F3P=0.0932, F4P= -0.0135,
     >          F5P=0.0361, F6P=0.0403, F7P= -0.0761, F8P = 0.0599, 
     >          F9P=0.0164, F10P=0.0453, F11P= -0.0265,
     >          G1P= -0.2867, G2P=0.2176, G3P=0.2826, G4P=-0.0644,
     >          G5P=0.0466, G6P=0.0973, G7P= -0.0081, G8P=-0.0235
      
      ! Testing if I need this instead for calculation
      real*8 cd_average
      common /user_fluct04/ cd_average

!
! For ppiclf_user_debug.f
!
      real*8 phimax,
     >         fqsx_max,fqsy_max,fqsz_max,
     >         famx_max,famy_max,famz_max,
     >         fdpdx_max,fdpdy_max,fdpdz_max,
     >         fcx_max,fcy_max,fcz_max,
     >         umean_max,vmean_max,wmean_max,
     >         fqs_mag,fam_mag,fdp_mag,fc_mag,
     >         fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,
     >         fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >         fvux_max,fvuy_max,fvuz_max,
     >         qq_max,tau_max,lift_max
      common /user_debug/ phimax,
     >         fqsx_max,fqsy_max,fqsz_max,
     >         famx_max,famy_max,famz_max,
     >         fdpdx_max,fdpdy_max,fdpdz_max,
     >         fcx_max,fcy_max,fcz_max,
     >         umean_max,vmean_max,wmean_max,
     >         fqs_mag,fam_mag,fdp_mag,fc_mag,
     >         fqsx_fluct_max,fqsy_fluct_max,fqsz_fluct_max,
     >         fqsx_total_max,fqsy_total_max,fqsz_total_max,
     >         fvux_max,fvuy_max,fvuz_max,
     >         qq_max,tau_max,lift_max

!
! For ppiclf_user_AddedMass.f
!
      integer*4 nneighbors
      real*8 Fam(3), FamUnary(3), FamBinary(3),
     >       Wdot_neighbor_mean(3), R_pair(6,6)

      common /user_AddedMass01/ nneighbors
      common /user_AddedMass02/ Fam, FamUnary, FamBinary,
     >                          Wdot_neighbor_mean, R_pair

!
! For ppiclf_solve_InitAngularPeriodic
!
      integer*4 x_per_flag, y_per_flag, z_per_flag, ang_per_flag,
     >          ang_case 
      real*8 ang_per_angle, ang_per_xangle,
     >       ang_per_rin, ang_per_rout,
     >       xrot(3) , vrot(3)
      real*8 x_per_min, x_per_max,
     >  y_per_min, y_per_max, z_per_min, z_per_max

      common /solve_InitAngularPeriodic01/ x_per_flag, y_per_flag, 
     >                                z_per_flag, ang_per_flag,
     >                                ang_case
      common /solve_InitAngularPeriodic02/ ang_per_angle, 
     >                                ang_per_xangle,
     >                                ang_per_rin, 
     >                                ang_per_rout,
     >                                xrot, vrot
      common /solve_InitAngularPeriodic03/ x_per_min, x_per_max,
     >  y_per_min, y_per_max, z_per_min, z_per_max
     

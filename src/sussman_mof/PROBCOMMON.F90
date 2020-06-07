#undef BL_LANG_CC
#define BL_LANG_FORT

#include "REAL.H"
#include "CONSTANTS.H"

#if (BL_SPACEDIM==3)
#define SDIM 3
#elif (BL_SPACEDIM==2)
#define SDIM 2
#else  
print *,"dimension bust"
stop
#endif

module probcommon_module

implicit none

! fort_initial_temperature added: April 10, 2018
! nucleation_init_time added: May 5, 2018
#include "probdataf95.H"

      REAL_T, PARAMETER :: GAMMA_SIMPLE_PARMS=1.4

         ! R=CP-CV
         ! CP=1.007D+7
         ! GAMMA=CP/CV=1.39861
      REAL_T, PARAMETER :: R_AIR_PARMS=0.287D+7
      REAL_T, PARAMETER :: CV_AIR_PARMS=0.72D+7
      REAL_T, PARAMETER :: PCAV_TAIT=220.2726D0  ! dyne/cm^2
      REAL_T, PARAMETER :: PCAV_TAIT_VACUUM=220.2726D0  ! dyne/cm^2
      REAL_T, PARAMETER :: A_TAIT=1.0D+6  ! dyne/cm^2
      REAL_T, PARAMETER :: B_TAIT=3.31D+9  ! dyne/cm^2
      REAL_T, PARAMETER :: RHOBAR_TAIT=1.0D0  ! g/cm^3
      REAL_T, PARAMETER :: GAMMA_TAIT=7.15D0

      INTEGER_T, PARAMETER :: visual_RT_transform=1

      INTEGER_T, PARAMETER :: bubbleInPackedColumn=1001
      INTEGER_T, PARAMETER :: DO_SANITY_CHECK=0
      INTEGER_T, PARAMETER :: COARSE_FINE_VELAVG=1
      REAL_T, PARAMETER :: MASK_FINEST_TOL=1.0D-3
      REAL_T, PARAMETER :: FACETOL=1.0D-3
      REAL_T, PARAMETER :: LSTOL=1.0D-2
      REAL_T, PARAMETER :: VOFTOL_SLOPES=1.0D-2
      REAL_T, PARAMETER :: VOFTOL=1.0D-8
      REAL_T, PARAMETER :: VOFTOL_AREAFRAC=1.0D-1
      INTEGER_T, PARAMETER :: POLYGON_LIST_MAX=400

       ! variables for comparing with Fred Stern's experiment
      INTEGER_T, PARAMETER :: SHALLOW_M=100
      INTEGER_T, PARAMETER :: SHALLOW_N=1000
      REAL_T, PARAMETER :: SHALLOW_TIME=12.0

      REAL_T shallow_water_data(0:SHALLOW_M,0:SHALLOW_N,2)
      REAL_T inflow_time(3000)
      REAL_T inflow_elevation(3000)
      REAL_T inflow_velocity(3000)
      REAL_T outflow_time(3000)
      REAL_T outflow_elevation(3000)
      REAL_T outflow_velocity(3000)
      INTEGER_T inflow_count,outflow_count
      INTEGER_T last_inflow_index,last_outflow_index

      INTEGER_T, PARAMETER :: recalesce_num_state=6
      INTEGER_T recalesce_material(100)
      REAL_T recalesce_state_old(recalesce_num_state*100)

       ! variables from "rfiledata"
      REAL_T zstatic(0:300),rstatic(0:300) 

       ! variables from "pressure_bcs"
      real*8  dt_pressure_bcs
      real*8  time_pressure_bcs(0:100) ,  pressbc_pressure_bcs(0:100,1:3)
      INTEGER_T selectpress
       ! variables from "vel_bcs"
      real*8  timehist_velbc(0:100), &
        zpos_velbc(1:50),velbc_velbc(0:100,1:50), &
        period_velbc,rigidwall_velbc
      INTEGER_T itime_velbc,ipos_velbc

        ! level,index,dir
      REAL_T, allocatable, dimension(:,:,:) :: grid_cache
      INTEGER_T cache_index_low,cache_index_high,cache_max_level
      INTEGER_T :: grid_cache_allocated=0

contains

end module probcommon_module


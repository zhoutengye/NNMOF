
!------------------------------------------------------RUDMAN INITIALIZATION---------------------------------------

module module_rudman_init
 use module_grid
 use module_flow, only:s
 implicit none
 real(8), dimension(:,:,:), allocatable, target :: cvof_f ! fine grid vof stored
 integer, dimension(:,:,:), allocatable :: vof_flag_f     ! fine grid vof flags stored
 real(8), dimension(:), allocatable :: x_f, xh_f, dx_f, dxh_f
 real(8), dimension(:), allocatable :: y_f, yh_f, dy_f, dyh_f
 real(8), dimension(:), allocatable :: z_f, zh_f, dz_f, dzh_f
 integer, dimension(3,2) :: idx_f , idx_c ! save fine and coarse grid indexes
 logical :: Rudman
 integer, parameter :: coarse = 1 , fine = 2
 logical :: subgrid  = .false.
 real(8), dimension(2) :: cell_size = 0.d0 ! cell_size(1) --> coarse , cell_size(2) --> fine
 save
! write functions to check global flags and place error checks

contains 

  function interpole_quad(val1,val2,val3,delta_x)
  implicit none
  real (8) :: interpole_quad, val1, val2, val3, delta_x
  real (8) :: a,b,c
  
  a = 0.d0 ; b = 0.d0 ; c = 0.d0 
  
  a = 0.5d0*((val1 + val3) - 2*val2)
  b = 0.5d0*(val3 - val1)
  c = val2
  
  interpole_quad = a*(delta_x**2) + b*(delta_x) + c 
  
  end function interpole_quad 



 subroutine setup_coord() ! sets fine or coarse coordinates according to current bounds     
    implicit none
    integer :: i,j,k
     
     do i=1,Nxt; s=dfloat(i-Ng)/dfloat(Nx); xh(i)=xLength*(xform*s*(0.5d0-s)*(1d0-s)+s); enddo
     do j=1,Nyt; s=dfloat(j-Ng)/dfloat(Ny); yh(j)=yLength*(yform*s*(0.5d0-s)*(1d0-s)+s); enddo
     do k=1,Nzt; s=dfloat(k-Ng)/dfloat(Nz); zh(k)=zLength*(zform*s*(0.5d0-s)*(1d0-s)+s); enddo

     do i=2,Nxt; x(i)=0.5d0*(xh(i)+xh(i-1)); enddo; x(1)=2d0*xh(1)-x(2)
     do j=2,Nyt; y(j)=0.5d0*(yh(j)+yh(j-1)); enddo; y(1)=2d0*yh(1)-y(2)
     do k=2,Nzt; z(k)=0.5d0*(zh(k)+zh(k-1)); enddo; z(1)=2d0*zh(1)-z(2)

     do i=1,Nxt-1; dxh(i)=x(i+1)-x(i); enddo; dxh(Nxt)=dxh(Nxt-1)
     do j=1,Nyt-1; dyh(j)=y(j+1)-y(j); enddo; dyh(Nyt)=dyh(Nyt-1)
     do k=1,Nzt-1; dzh(k)=z(k+1)-z(k); enddo; dzh(Nzt)=dzh(Nzt-1)
     
     do i=2,Nxt; dx(i)=xh(i)-xh(i-1); enddo; dx(1)=dx(2);
     do j=2,Nyt; dy(j)=yh(j)-yh(j-1); enddo; dy(1)=dy(2);
     do k=2,Nzt; dz(k)=zh(k)-zh(k-1); enddo; dz(1)=dz(2);

 end subroutine setup_coord


 subroutine init_rudman_coord(val)    
    implicit none
    integer :: i,j,k
    real(8) :: s
    real(8), intent(out) :: val

     val = 0.d0
    
     call update_index_fine() ! changes Nx,Nxt etc and indices to fine grid
     
     allocate(x_f(Nxt), xh_f(Nxt), dx_f(Nxt), dxh_f(Nxt), &
              y_f(Nyt), yh_f(Nyt), dy_f(Nyt), dyh_f(Nyt), &   ! allocation of co-ordinate arrays for fine grid
              z_f(Nzt), zh_f(Nzt), dz_f(Nzt), dzh_f(Nzt)  )   
     
     do i=1,Nxt; s=dfloat(i-Ng)/dfloat(Nx); xh_f(i)=xLength*(xform*s*(0.5d0-s)*(1d0-s)+s); enddo
     do j=1,Nyt; s=dfloat(j-Ng)/dfloat(Ny); yh_f(j)=yLength*(yform*s*(0.5d0-s)*(1d0-s)+s); enddo
     do k=1,Nzt; s=dfloat(k-Ng)/dfloat(Nz); zh_f(k)=zLength*(zform*s*(0.5d0-s)*(1d0-s)+s); enddo

     do i=2,Nxt; x_f(i)=0.5d0*(xh_f(i)+xh_f(i-1)); enddo; x_f(1)=2d0*xh_f(1)-x_f(2)
     do j=2,Nyt; y_f(j)=0.5d0*(yh_f(j)+yh_f(j-1)); enddo; y_f(1)=2d0*yh_f(1)-y_f(2)
     do k=2,Nzt; z_f(k)=0.5d0*(zh_f(k)+zh_f(k-1)); enddo; z_f(1)=2d0*zh_f(1)-z_f(2)

     do i=1,Nxt-1; dxh_f(i)=x_f(i+1)-x_f(i); enddo; dxh_f(Nxt)=dxh_f(Nxt-1)
     do j=1,Nyt-1; dyh_f(j)=y_f(j+1)-y_f(j); enddo; dyh_f(Nyt)=dyh_f(Nyt-1)
     do k=1,Nzt-1; dzh_f(k)=z_f(k+1)-z_f(k); enddo; dzh_f(Nzt)=dzh_f(Nzt-1)
     
     do i=2,Nxt; dx_f(i)=xh_f(i)-xh_f(i-1); enddo; dx_f(1)=dx_f(2);
     do j=2,Nyt; dy_f(j)=yh_f(j)-yh_f(j-1); enddo; dy_f(1)=dy_f(2);
     do k=2,Nzt; dz_f(k)=zh_f(k)-zh_f(k-1); enddo; dz_f(1)=dz_f(2);

     val = dxh_f(is)

     call revert_index_coarse() ! sets back Nx etc and indices to coarse grid

 end subroutine init_rudman_coord

 
 subroutine init_rudman_vof()
    implicit none
   
      allocate(cvof_f(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)))       

      allocate(vof_flag_f(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)))       

      cvof_f = 0.d0 ; vof_flag_f = 0   

 end subroutine init_rudman_vof

end module module_rudman_init

!----------------------------------------------------------END INITIALIZATION---------------------------------------


!----------------------------------------------RUDMAN GHOST PROCEDURES-------------------------------------------
module module_rudman_ghost
 use module_grid
 use module_rudman_init , only: idx_f, idx_c
 implicit none

contains

 subroutine do_all_ghost_fine(var)
   implicit none
   real(8), dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: var
   integer, parameter :: ngh=2
   include 'mpif.h'
   integer :: req(48),sta(MPI_STATUS_SIZE,48)
   integer :: ierr
 
      call ghost_x_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ghost_y_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ghost_z_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
 
 end subroutine do_all_ghost_fine
 
 
 subroutine do_all_ighost_fine(var)
   implicit none
   integer, dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: var 
   integer, parameter :: ngh=2
   include 'mpif.h'
   integer :: req(48),sta(MPI_STATUS_SIZE,48)
   integer :: ierr

      call ighost_x_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_y_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
      call ighost_z_fine(var,ngh,req(1:4)); call MPI_WAITALL(4,req(1:4),sta(:,1:4),ierr)
 
 end subroutine do_all_ighost_fine

end module module_rudman_ghost
!-----------------------------------------------------------END GHOST PROCEDURES----------------------------------


!--------------------------------------------------------SUBGRID VOF ADVECTION--------------------------------------
module module_rudman_advection
 use module_grid  
 use module_flow
 use module_bc
 use module_rudman_init , only: idx_f,idx_c,cvof_f,vof_flag_f 
 implicit none
 real(8), dimension(:,:,:,:), allocatable :: vel_fine
 real(8), dimension(:,:,:,:), allocatable :: vel_coarse
 real(8), dimension(:,:,:), allocatable :: vofmask
 logical :: t0 = .true.
 save

contains

 subroutine init_adv(vof)
   implicit none
   include 'mpif.h'
   real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: vof
   integer :: i,j,k,iini,iend,jini,jend,kini,kend

   if(t0) then
    t0 = .false.  
    allocate(vel_fine(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2),6))
    allocate(vel_coarse(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),6))
    allocate(vofmask(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)))
   endif 
  
   ! Million dollar droplet initialization

  ! iini=idx_c(1,1) + Ng ; iend=idx_c(1,2) - Ng
  ! jini=idx_c(2,1) + Ng ; jend=idx_c(2,2) - Ng
  ! kini=idx_c(3,1) + Ng ; kend=idx_c(3,2) - Ng 

  ! do i=iini,iend ; do j=jini,jend ; do k=kini,kend
  ! 
  !    if(vof(i,j,k)==0.d0 ) then
  !      u(i,j,k) = 0.d0      
  !      v(i,j,k) = 0.d0
  !      w(i,j,k) = 0.d0
  !    else
  !      u(i,j,k) = 1.d0      
  !      v(i,j,k) = 0.d0
  !      w(i,j,k) = 0.d0
  !    endif

  ! enddo ; enddo ; enddo

   vel_coarse(:,:,:,1) = u  ;  vel_coarse(:,:,:,4) = umask 
   vel_coarse(:,:,:,2) = v  ;  vel_coarse(:,:,:,5) = vmask    
   vel_coarse(:,:,:,3) = w  ;  vel_coarse(:,:,:,6) = wmask 

   vel_fine = 0.d0
   vofmask = 0.d0

   

 end subroutine init_adv


 subroutine project_coarse2fine()
    use module_rudman_ghost
    implicit none
    include 'mpif.h'
    integer :: dir,i

    ! u,v,w are coarse grid and inherited from module_flow through rudman_advection
    ! Fine grid containers for u,v,w, --> vel_fine are initialized with values of u,v,w
    ! by projecting them onto the subgrid

   ! print * , "BEFORE PROJECTION--------------"
   ! print * , u(22,17,22), v(22,17,22) , w(22,17,22)
   ! print * , u(22,28,22), v(22,28,22) , w(22,28,22)

    call compute_vel_coarse2fine(1,vel_fine(:,:,:,1),u)
    call compute_vel_coarse2fine(2,vel_fine(:,:,:,2),v)
    call compute_vel_coarse2fine(3,vel_fine(:,:,:,3),w)

   ! All solid masks for the fine grid are set to 1 so we don't deal with embedded solids    

    vel_fine(:,:,:,4) = 1.d0 ; vel_fine(:,:,:,5) = 1.d0 ; vel_fine(:,:,:,6) = 1.d0 

   ! Implementing fine grid ghost conditions on fine grid velocities

    do i=1,6
      call do_all_ghost_fine(vel_fine(:,:,:,i))
    enddo

    ! boundary conditions for fine grid velocities applied

    call vel_bc_rudman(vel_fine,time,dt/2.d0,0,idx_f)
  
 end subroutine project_coarse2fine


 subroutine vofsweeps_fine(tswap,t)
     implicit none
     include 'mpif.h'
     integer, intent(in) :: tswap
     real(8) , intent(in) :: t
     real(8) , dimension(idx_f(1,1):idx_f(1,2),idx_f(2,1):idx_f(2,2),idx_f(3,1):idx_f(3,2)) :: uf,vf,wf
     save  

     vofmask = 0.d0
     call compute_wy_mask(cvof_f,vofmask,idx_f)

     uf = vel_fine(:,:,:,1) ; vf = vel_fine(:,:,:,2) ; wf = vel_fine(:,:,:,3) 
     !uf = 1.d0 ; vf = 0.d0 ; wf = 0.d0 

     if (MOD(tswap,3).eq.0) then  ! do z x y 
        call advect_vof_fine_centered(wf,cvof_f,vof_flag_f,3,vofmask,dt)
        call advect_vof_fine_centered(uf,cvof_f,vof_flag_f,1,vofmask,dt)
        call advect_vof_fine_centered(vf,cvof_f,vof_flag_f,2,vofmask,dt)

     elseif (MOD(tswap,3).eq.1) then ! do y z x
        call advect_vof_fine_centered(vf,cvof_f,vof_flag_f,2,vofmask,dt)
        call advect_vof_fine_centered(wf,cvof_f,vof_flag_f,3,vofmask,dt)
        call advect_vof_fine_centered(uf,cvof_f,vof_flag_f,1,vofmask,dt)

     else ! do x y z
        call advect_vof_fine_centered(uf,cvof_f,vof_flag_f,1,vofmask,dt)
        call advect_vof_fine_centered(vf,cvof_f,vof_flag_f,2,vofmask,dt)
        call advect_vof_fine_centered(wf,cvof_f,vof_flag_f,3,vofmask,dt)

    endif


    
 end subroutine vofsweeps_fine


 subroutine vofandmomsweepsstaggered_fine(tswap,t)
    use module_rudman_init
    implicit none
    include 'mpif.h'
    integer :: dir
    integer, intent(in) :: tswap
    real(8), intent(in) :: t
    integer, dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)) :: v_flag ! staggered vof flag
    real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)):: mom ! staggered
    real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2)):: vof ! staggered
    real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),2):: flux ! stg coarse vof fluxes
    real(8), dimension(idx_c(1,1):idx_c(1,2),idx_c(2,1):idx_c(2,2),idx_c(3,1):idx_c(3,2),2):: wy_mask  ! wy coefficients

   ! flux(:,:,:,1) -> vof1 flux
   ! flux(:,:,:,2) -> vof3 flux
   ! wy_mask(:,:,:,1) --> VOF mask wy  ,  wy_mask(:,:,:,2) --> MOM mask wy

   do dir=1,3
 
      flux = 0.d0 ; mom = 0.d0 ; vof = 0.d0 ; v_flag = 0

      wy_mask = 0.d0 
 
      call compute_mom_vof_coarse_stg(mom,vof,v_flag,dir)
        

      wy_mask(:,:,:,2) = vel_coarse(:,:,:,dir)  ! assignment of vel at previous time-step as WY MOM coefficient

      call compute_wy_mask(vof,wy_mask(:,:,:,1),idx_c) ! compute WY VOF coefficient 


      if (MOD(tswap,3).eq.0) then ! z x y

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,3),cvof_f,flux,3,dt/2.d0)  
         call advect_mom_vof_coarse_stg(w,wmask,mom,vof,v_flag,flux,wy_mask,3,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,1),cvof_f,flux,1,dt/2.d0)  
         call advect_mom_vof_coarse_stg(u,umask,mom,vof,v_flag,flux,wy_mask,1,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,2),cvof_f,flux,2,dt/2.d0)  
         call advect_mom_vof_coarse_stg(v,vmask,mom,vof,v_flag,flux,wy_mask,2,dir,t,dt)        


      elseif (MOD(tswap,3).eq.1) then ! y z x

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,2),cvof_f,flux,2,dt/2.d0)  
         call advect_mom_vof_coarse_stg(v,vmask,mom,vof,v_flag,flux,wy_mask,2,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,3),cvof_f,flux,3,dt/2.d0)  
         call advect_mom_vof_coarse_stg(w,wmask,mom,vof,v_flag,flux,wy_mask,3,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,1),cvof_f,flux,1,dt/2.d0)  
         call advect_mom_vof_coarse_stg(u,umask,mom,vof,v_flag,flux,wy_mask,1,dir,t,dt)        

      else ! x y z

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,1),cvof_f,flux,1,dt/2.d0)  
         call advect_mom_vof_coarse_stg(u,umask,mom,vof,v_flag,flux,wy_mask,1,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,2),cvof_f,flux,2,dt/2.d0)  
         call advect_mom_vof_coarse_stg(v,vmask,mom,vof,v_flag,flux,wy_mask,2,dir,t,dt)        

         call compute_vof_flux_coarse_stg(vel_fine(:,:,:,3),cvof_f,flux,3,dt/2.d0)  
         call advect_mom_vof_coarse_stg(w,wmask,mom,vof,v_flag,flux,wy_mask,3,dir,t,dt)        

      endif

     if (dir.eq.1) then
        call rudman_get_velocity_coarse(mom,vof,rho1,rho2,u,du,dt,1)
     elseif (dir.eq.2) then
        call rudman_get_velocity_coarse(mom,vof,rho1,rho2,v,dv,dt,2)
     elseif (dir.eq.3) then
        call rudman_get_velocity_coarse(mom,vof,rho1,rho2,w,dw,dt,3)
     endif

   enddo


 end subroutine vofandmomsweepsstaggered_fine


end module module_rudman_advection



!----------------------------------------------------END SUBGRID VOF ADVECTION---------------------------------------------

!================================================================================================
!=================================================================================================
! Paris-0.1
!
! Thermal energy conservation and phase change extensions
! written by Leon Malan
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation; either version 2 of the
! License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
! General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
! 02111-1307, USA.  
!=================================================================================================
!=================================================================================================
#ifdef PHASE_CHANGE
subroutine init_phase_change
  use module_boil
  implicit none
  call ReadHeatParameters
  call initialize_heat
end subroutine init_phase_change

!Read parameters
subroutine ReadHeatParameters
  use module_boil
  use module_grid
  use module_2phase
  use module_IO
  use module_flow
  implicit none
  integer :: ierr,BC,line
  real(8) :: R0, Ja, T_inf, t0 ! Sup_liquid test 
  namelist /heatparameters/ Cp1, Cp2, h_fg, kc1, kc2, T_sat, T0_1, T0_2,&
       BDRY_T, BC_T, ec1, ec2, q_set, type, T_ref, setup, adv_method,&
       perturb_T, epsT, epsBL, EtAdvScheme, HYPRE_Ps, HYPRE_Timp

  Cp1 = 1.0d0; Cp2 = 10.0; h_fg=1.0d30; 
  kc1=1.0d0; kc2=1.0d3; T_sat=1.0d0
  T0_1=1.0d1; T0_2=1.0d1
  BDRY_T = 1
  BC_T = 0.0d0
  q_set = 1.0d0
  type = 0
  T_ref = T0_1
  setup = 'standard'
  EtAdvScheme = 'WENO5'
  super_vapour=.false.
  diff_tsat_1D=.false.
  ngradT_test=.false.
  TdiffMods=.false.
  adv_method=1 ! 0:straight discretisation w WENO. 1:Consistent VOF approach
  perturb_T = .false.
  epsT = 1.0d0
  epsBL= 1.0d-2
  HYPRE_Ps=.false.
  HYPRE_Timp=.false.
  
  open(unit=10, file='inputheat', status='old', action='read', iostat=ierr)
  if (ierr .ne. 0) call err_no_out_dir("ReadHeatParameters: error opening 'inputheat' file --- perhaps it does not exist ?")
  read(10,NML=heatparameters)
  close(10)

  if (.not.(HYPRE) .and. (HYPRE_Ps .or. HYPRE_Timp)) then
     write(*,*)'Main HYPRE must be selected to use HYPRE for phase change pressure and implicit temperature'
     HYPRE_Ps=.false.
     HYPRE_Timp=.false.
  endif

  do BC=1,6
     if (.not.(BDRY_T(BC)==0 .or. BDRY_T(BC)==1 .or. BDRY_T(BC)==2) ) then
        write(*,'("Error, Temp BC must be either 0 (Dirichlet) or 1 (Von Neumann) or 2 (periodic)")')
        call pariserror("Invalid temperature BC")
     endif
  enddo

  if (.not.(setup=='standard' &
       .or. setup=='sup_vapour' &
       .or. setup=='diff_tsat' &
       .or. setup=='normal_gradT' &
       .or. setup=='ImpTMods' &
       .or. setup=='sup_liq')) then
     write(*,*)'Energy setup type in inputheat invalid, reverting to standard'
     setup='standard'
  endif
  
  if (.not.(EtAdvScheme=='WENO5' &
       .or. EtAdvScheme=='WENO' &
       .or. EtAdvScheme=='ENO' &
       .or. EtAdvScheme=='Superbee')) then
     write(*,*)'Energy advection scheme type invalid, reverting to WENO5'
     EtAdvScheme='WENO5'
  endif
       
  if (setup=='sup_vapour') super_vapour=.true.
  if (setup=='ImpTMods') TdiffMods=.true.
  if (setup=='sup_liq') then
     super_liquid=.true.
     if (.not. restart) then
        open(unit=11, file='inputsupliq', status='old', action='read', iostat=ierr)
        if (ierr .ne. 0) call err_no_out_dir("Error opening 'inputsupliq' file --- perhaps it does not exist ?")
        read(11,*)R0, t_0, Ja, T_inf
        close(11)
        if (ABS(T0_2-T_inf)>1e-10) then
           if (rank==0) write(*,'("Inconsistent temp for sup liquid, changing T0_2 to T_inf: ",e14.5)')T_inf
           T0_2 = T_inf
        endif
        if (ABS(rad(1)-R0)>1e-10) then
           if (rank==0) write(*,'("Inconsistent initial radius for sup liquid, changing rad(1) to R0: ",e14.5)')R0
           rad(1) = R0 
        endif
        if (rank==0) then
           write(*,'("Running bubble in superheated liquid with Ja = ",e14.5," and T_inf: ",e14.5," with R0: ",e14.5)')&
                Ja,T0_2, rad(1)
        endif
        allocate(Tprof(1:101),R_tprof(1:101))
        open(unit=12, file='t_profile', status='old', action='read', iostat=ierr)
        if (ierr .ne. 0) call err_no_out_dir("Error opening 't_profile' file --- perhaps it does not exist ?")
        do line=1,101
           read(12,*)R_tprof(line), Tprof(line)
        enddo
        close(12)
     endif
  endif
  if (setup=='normal_gradT') ngradT_test=.true.
  if (setup=='diff_tsat') then
     diff_tsat_1D =.true.
     type=0
  endif

  !add adv_method logical check
  if (.not.(adv_method==0 .or. adv_method==1)) &
       call pariserror("Invalid adv_method, use '0' for standard discretisation and '1' for VOF consistent scheme.")

end subroutine ReadHeatParameters

subroutine initialize_heat
  use module_grid
  use module_flow    
  use module_boil
  implicit none
  
  ! Allocate variables
  allocate ( Te(imin:imax,jmin:jmax,kmin:kmax),&
       s_v(imin:imax,jmin:jmax,kmin:kmax),kc(imin:imax,jmin:jmax,kmin:kmax),&
       rho_cp(imin:imax,jmin:jmax,kmin:kmax), energy(imin:imax,jmin:jmax,kmin:kmax),&
       topo_mask(imin:imax,jmin:jmax,kmin:kmax), mdot(imin:imax,jmin:jmax,kmin:kmax),&
       dc(imin:imax,jmin:jmax,kmin:kmax), Tumask(imin:imax,jmin:jmax,kmin:kmax), &
       Tvmask(imin:imax,jmin:jmax,kmin:kmax), Twmask(imin:imax,jmin:jmax,kmin:kmax),&
       u_l(imin:imax,jmin:jmax,kmin:kmax), v_l(imin:imax,jmin:jmax,kmin:kmax), &
       w_l(imin:imax,jmin:jmax,kmin:kmax)  )

  !test stuff: u_l
  if (type>0 .or. diff_tsat_1D ) then
     allocate(  u_s(imin:imax,jmin:jmax,kmin:kmax), v_s(imin:imax,jmin:jmax,kmin:kmax), &
          w_s(imin:imax,jmin:jmax,kmin:kmax), p_s(imin:imax,jmin:jmax,kmin:kmax) )
     u_s=0.0d0; v_s=0.0d0; w_s=0.0d0
     p_s=0.0d0
  endif

  if (itime_scheme==2) then
     allocate( Te_old(imin:imax,jmin:jmax,kmin:kmax), kc_old(imin:imax,jmin:jmax,kmin:kmax),&
          rho_cp_old(imin:imax,jmin:jmax,kmin:kmax), energy_old(imin:imax,jmin:jmax,kmin:kmax) )
  endif

  u_l=0.0d0; v_l=0.0d0; w_l=0.0d0
  Tumask = 1.0d0; Tvmask = 1.0d0; Twmask = 1.0d0
  inv_drho = 1.d0/rho1 - 1.d0/rho2
  mdot=0.d0; s_v=0.0d0; dc=0.0d0
end subroutine initialize_heat

subroutine initcondition_heat
  use module_grid
  use module_VOF
  use module_flow
  use module_boil
  use module_2phase
  implicit none   
  integer :: i,j,k,dir,dor,line
  real(8) :: x_scale, T_th, c1, c2, gap, T_sup, r, dr, Tint
  logical :: found_dir, not_interpolated

  if (super_vapour) then !need to make general for all orientations
     found_dir=.false.
     if (.not.test_plane) call pariserror("For superheated vapour test, the VOF test type must be set to 'plane'.")
     !get dir from plane constants
     do dir=1,3
        if (ABS(ABS(n_p(dir))-1.0d0)<EPSC ) then
           found_dir=.true.
           dor=dir
        endif
     enddo
     if (ABS(ABS(n_p(1))+ABS(n_p(2))+ABS(n_p(3))-1.0d0)>1.d-14 .or. .not.(found_dir)) &
          call pariserror("For superheated vapour test, set plane normal to a single coordinate direction")
     Select case(dor)
     case(1)
        do i=is,ie
           gap=MAX(n_p(dor)*xLength,0.d0)-plane
           T_sup=BC_T(dor+3*INT( 0.5d0*(n_p(dor)+1)) )-T_sat
           if (ABS(gap)<EPSC) call pariserror("Geometry error in gap for superheated vapour test case")
           x_scale= MAX((n_p(dor)*x(i)-plane)/gap,0.d0)
           Te(i,js:je,ks:ke) = T_sat+x_scale*T_sup
           if (plane-n_p(dor)*xh(i)>0.d0) then
              u(i,js:je,ks:ke)=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
           else
              u(i,js:je,ks:ke)=0.d0
           endif
           u_l=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
        enddo
        if (perturb_T) then
           do j=js,je; do k=ks,ke
              Te(:,j,k) = Te(:,j,k) + &
                   epsT*0.5d0*(2.0d0-COS(2.0d0*pi*y(j)/yLength)-COS(2.0d0*pi*z(k)/zLength) )
           enddo; enddo
        endif
     case(2)
        do j=js,je
           gap=MAX(n_p(dor)*yLength,0.d0)-plane
           T_sup=BC_T(dor+3*INT( 0.5d0*(n_p(dor)+1)) )-T_sat
           if (ABS(gap)<EPSC) call pariserror("Geometry error in gap for superheated vapour test case")
           x_scale= MAX((n_p(dor)*y(j)-plane)/gap,0.d0)
           Te(is:ie,j,ks:ke) = T_sat+x_scale*T_sup
           if (plane-n_p(dor)*yh(j)>0.d0) then
              v(is:ie,j,ks:ke)=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
           else
              v(is:ie,j,ks:ke)=0.d0
           endif
           v_l=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
        enddo
        if (perturb_T) then
           do i=is,ie; do k=ks,ke
              Te(i,:,k) = Te(i,:,k) + &
                   epsT*0.5d0*(2.0d0-COS(2.0d0*pi*x(i)/xLength) - COS(2.0d0*pi*z(k)/zLength) )
           enddo; enddo
        endif
     case(3)
        do k=ks,ke
           gap=MAX(n_p(dor)*zLength,0.d0)-plane
           T_sup=BC_T(dor+3*INT( 0.5d0*(n_p(dor)+1)) )-T_sat
           if (ABS(gap)<EPSC) call pariserror("Geometry error in gap for superheated vapour test case")
           x_scale= MAX((n_p(dor)*z(k)-plane)/gap,0.d0)
           Te(is:ie,js:je,k) = T_sat+x_scale*T_sup
           if (plane-n_p(dor)*zh(k)>0.d0) then
              w(is:ie,js:je,k)=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
           else
              w(is:ie,js:je,k)=0.d0
           endif
           w_l=-1.0d0*n_p(dor)*inv_drho*kc1*T_sup/(gap*h_fg)
        enddo
        if (perturb_T) then
           do i=is,ie; do j=js,je
              Te(i,j,:) = Te(i,j,:) + &
                   epsT*0.5d0*(2.0d0-COS(2.0d0*pi*x(i)/xLength)-COS(2.0d0*pi*y(j)/yLength) )
           enddo; enddo
        endif
     end select
  else if (ngradT_test) then
     if (test_droplet) then
        type=2
        kc1=1.0d0
        kc2=1.0d0
        h_fg=1.0d0
        T_sat=1.0d1
        do i=is,ie; do j=js,je;  do k=ks,ke
           r = sqrt( (x(i)-xc(1))**2.0d0 + (y(j)-yc(1))**2.0d0 +  (z(k)-zc(1))**2.0d0 )
           Te(i,j,k) = r/rad(1) * T_sat
        enddo; enddo; enddo
     else if (test_plane) then
        call pariserror('Still to implement test for plane')
     else
        call pariserror('VOF test type for interface normal temp gradient test must be plane or droplet')
     endif
  else if (super_liquid) then
     T_sup=T0_2-T_sat  
     if (rank==0) then
        write(*,'("Running superheated liquid test")') 
        write(*,'("T_sup: ",e14.5," T_sat: ",e14.5)')T_sup,T_sat
     endif
     !dr=rad(1)*epsBL
     do i=is,ie; do j=js,je;  do k=ks,ke
        r = sqrt((x(i)-xc(1))**2d0 + (y(j)-yc(1))**2d0 + (z(k)-zc(1))**2d0)
        if (r <= rad(1)) then
           Te(i,j,k) = T_sat
        else
           not_interpolated=.true.
           line=1
           do while (not_interpolated)
              if (r >= R_tprof(line) .and. r <= R_tprof(line+1) ) then
                 Tint = Tprof(line) + (r-R_tprof(line)) / (R_tprof(line+1)-R_tprof(line)) * &
                      (Tprof(line+1)-Tprof(line))
                 not_interpolated = .false.
              else
                 line=line+1
              endif
              if (line==101) then
                 write(*,*)r,line,R_tprof(line-1),R_tprof(line)
                 call pariserror('Error in Sup liquid T0 profile generation')
              endif
           enddo
           Te(i,j,k) = Tint
        endif
     enddo; enddo; enddo
  else
     if (TwoPhase) then
        if (DoVof) then
           do i=is,ie; do j=js,je;  do k=ks,ke
              Te(i,j,k) = cvof(i,j,k)*T0_2 + (1.0d0-cvof(i,j,k))*T0_1
           enddo; enddo; enddo
        endif
     else
        Te = T0_1
     endif
  endif

  if (diff_tsat_1D) then
     if (rank==0) then
        OPEN(unit=86,file='ref_solution.txt')
        c1=T_sat-BC_T(1)*plane
        c2=T_sat-BC_T(4)*plane
        do i = is,nx+ng
           if (x(i)<plane) then
              T_th = BC_T(1)*x(i)+c1
           else
              T_th = BC_T(4)*x(i)+c2
           endif
           write(86,'(5e17.8)')x(i),T_th
        enddo
        close(86)
     endif
  endif

end subroutine initcondition_heat

!Heat diffusion
subroutine TempDiffusion(dt) !! FIXME: Need to add masks, Poisson BC's
  use module_grid
  use module_boil
  use module_VOF
  implicit none
  real(8), intent(in) :: dt
  integer :: d,i,j,k

  do i=is,ie; do j=js,je;  do k=ks,ke
     energy(i,j,k) = energy(i,j,k) + dt*(&
          (0.5d0*( kc(i+1,j,k)+kc(i,j,k) )*(Te(i+1,j,k)-Te(i,j,k))/dxh(i) -&
          0.5d0*(kc(i,j,k)+kc(i-1,j,k))*(Te(i,j,k)-Te(i-1,j,k))/dxh(i-1) )/dx(i) + &
          ( 0.5d0*(kc(i,j+1,k)+kc(i,j,k))*(Te(i,j+1,k)-Te(i,j,k))/dyh(j) -&
          0.5d0*(kc(i,j,k)+kc(i,j-1,k))*(Te(i,j,k)-Te(i,j-1,k))/dyh(j-1) )/dy(j) + &
          ( 0.5d0*(kc(i,j,k+1)+kc(i,j,k))*(Te(i,j,k+1)-Te(i,j,k))/dzh(k) -&
          0.5d0*(kc(i,j,k)+kc(i,j,k-1))*(Te(i,j,k)-Te(i,j,k-1))/dzh(k-1) )/dz(k) ) 
  enddo; enddo; enddo
end subroutine TempDiffusion

subroutine TempDiffusionImplicit 
  use module_grid
  use module_boil
  use module_VOF
  use module_flow
  use module_poisson
  use module_mgsolver
  use module_BC
  use module_surface_tension
  implicit none
  include 'mpif.h'
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: Tpass
  real(8), dimension(is:ie,js:je,ks:ke,8) :: coeff
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,0:1) :: xmod,ymod,zmod
  real(8), dimension(-1:0,1:3) :: dmod
  real(8) :: dh, kcc(0:1),lim_cut, delta, l2var!, delta_rhocp
  integer :: i,j,k,phase,d,b
  integer :: iters, ierr, stag_calls
  real(8) :: res_tabTimp(3), DivTolTimp, residual !, MaxEtLoc, MaxEt
  integer :: reqd(24),stat(MPI_STATUS_SIZE,24)
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: cstagx,cstagy,cstagz
  lim_cut=1.0d-2
  DivTolTimp=1.0d4
  dh = dx(is) !Structured, equi-sized grids only for now
  xmod=dh; ymod=dh; zmod=dh
  kcc(0)=kc1; kcc(1)=kc2
  !delta_rhocp = rho2*cp2-rho1*cp1

  cstagx=0.0d0; cstagy=0.0d0; cstagz=0.0d0
  do d=1,3
     select case(d)
     case(1)
        call get_staggered_fractions(cstagx,d)
     case(2)
        call get_staggered_fractions(cstagy,d)
     case(3)
        call get_staggered_fractions(cstagz,d)
     end select
  enddo
  
  !call get_all_heights(itimestep) ! Try to recalculate instead of passing in mem
  
  if (TdiffMods) then
     stag_calls=0
     OPEN(unit=70,file='mods_xy.txt',position='append')
     OPEN(unit=71,file='mods_xz.txt',position='append')
     OPEN(unit=74,file='mods_xy_h3.txt',position='append')
  endif
  
  do k=ks,ke; do j=js,je; do i=is,ie ! Get interface distances using Heights 
     ! or staggered cell reconstructions

     !Check neighbours
     if(vof_phase(i,j,k)==1) then
        if (vof_phase(i+1,j,k)==0) then
           if ( (height(i,j,k,1) >= lim_cut) .and. (height(i,j,k,1) <= 1.0d0-lim_cut) ) then
              xmod(i,j,k,1)=height(i,j,k,1)*dh
              xmod(i,j,k,0)=dh-xmod(i,j,k,1)
           else if ( (-1.0d0*height(i,j,k,2))>lim_cut .and. -1.0d0*height(i,j,k,2)<1.0d0-lim_cut) then
              xmod(i,j,k,1)=-1.0d0*height(i,j,k,2)*dh
              xmod(i,j,k,0)=dh-xmod(i,j,k,1)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              !write(*,'("Doing staggered cut in x",3I8)')i,j,k
              call staggered_cut_gen(cstagx,i,j,k,xmod(i,j,k,1),1,lim_cut)
              xmod(i,j,k,0)=MIN(MAX(dh-xmod(i,j,k,1),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif

           if (TdiffMods) then
              if (k==(ks+ke)/2) then
                 write(70,410)x(i),y(j),xmod(i,j,k,1),0.d0
                 write(70,410)x(i+1),y(j),-xmod(i,j,k,0),0.d0
              endif
              if (j==(js+je)/2) then
                 write(71,410)x(i),z(k),xmod(i,j,k,1),0.d0
                 write(71,410)x(i+1),z(k),-xmod(i,j,k,0),0.d0 
              endif
           endif
        endif
        if (vof_phase(i,j+1,k)==0) then
           if (height(i,j,k,3)>lim_cut .and. height(i,j,k,3)<1.0d0-lim_cut) then
              ymod(i,j,k,1)=height(i,j,k,3)*dh
              ymod(i,j,k,0)=dh-ymod(i,j,k,1)
           elseif (-height(i,j,k,4)>lim_cut .and. -height(i,j,k,4)<1.0d0-lim_cut) then
              ymod(i,j,k,1)=-height(i,j,k,4)*dh
              ymod(i,j,k,0)=dh-ymod(i,j,k,1)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              call staggered_cut_gen(cstagy,i,j,k,ymod(i,j,k,1),2,lim_cut)
              ymod(i,j,k,0)=MIN(MAX(dh-ymod(i,j,k,1),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif
           if (TdiffMods) then
              if (k==(ks+ke)/2) then
                 write(70,410)x(i),y(j),0.d0,ymod(i,j,k,1)
                 write(70,410)x(i),y(j+1),0.d0,-ymod(i,j,k,0)
              endif
           endif
        endif
        if (vof_phase(i,j,k+1)==0) then
           if (height(i,j,k,5)>lim_cut .and. height(i,j,k,5)<1.0d0-lim_cut) then
              zmod(i,j,k,1)=height(i,j,k,5)*dh
              zmod(i,j,k,0)=dh-zmod(i,j,k,1)
           elseif (-height(i,j,k,6)>lim_cut .and. -height(i,j,k,6)<1.0d0-lim_cut) then
              zmod(i,j,k,1)=-height(i,j,k,6)*dh
              zmod(i,j,k,0)=dh-zmod(i,j,k,1)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              call staggered_cut_gen(cstagz,i,j,k,zmod(i,j,k,1),3,lim_cut)
              zmod(i,j,k,0)=MIN(MAX(dh-zmod(i,j,k,1),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif
           if (TdiffMods) then
              if (j==(js+je)/2) then
                 write(71,410)x(i),z(k),0.d0,zmod(i,j,k,1)
                 write(71,410)x(i),z(k+1),0.d0,-zmod(i,j,k,0) 
              endif
           endif
        endif
     endif ! Liquid cell

     if(vof_phase(i,j,k)==0) then
        if (vof_phase(i+1,j,k)==1) then
           if (height(i,j,k,2)>lim_cut .and. height(i,j,k,2)<1.0d0-lim_cut) then
              xmod(i,j,k,0)=height(i,j,k,2)*dh
              xmod(i,j,k,1)=dh-xmod(i,j,k,0)
           else if (-height(i,j,k,1)>lim_cut .and. -height(i,j,k,1)<1.0d0-lim_cut) then
              xmod(i,j,k,0)=-height(i,j,k,1)*dh
              xmod(i,j,k,1)=dh-xmod(i,j,k,0)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              call staggered_cut_gen(cstagx,i,j,k,xmod(i,j,k,0),1,lim_cut)
              xmod(i,j,k,1)=MIN(MAX(dh-xmod(i,j,k,0),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif
           if (TdiffMods) then
              if (k==(ks+ke)/2) then
                 write(70,410)x(i),y(j),xmod(i,j,k,0),0.d0
                 write(70,410)x(i+1),y(j),-xmod(i,j,k,1),0.d0
              endif
              if (j==(js+je)/2) then
                 write(71,410)x(i),z(k),xmod(i,j,k,0),0.d0
                 write(71,410)x(i+1),z(k),-xmod(i,j,k,1),0.d0
              endif
           endif
        endif
        if (vof_phase(i,j+1,k)==1) then
           if (height(i,j,k,4)>lim_cut .and. height(i,j,k,4)<1.0d0-lim_cut) then
              ymod(i,j,k,0)=height(i,j,k,4)*dh
              ymod(i,j,k,1)=dh-ymod(i,j,k,0)
           else if (-height(i,j,k,3)>lim_cut .and. -height(i,j,k,3)<1.0d0-lim_cut) then
              ymod(i,j,k,0)=-height(i,j,k,3)*dh
              ymod(i,j,k,1)=dh-ymod(i,j,k,0)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              call staggered_cut_gen(cstagy,i,j,k,ymod(i,j,k,0),2,lim_cut)
              ymod(i,j,k,1)=MIN(MAX(dh-ymod(i,j,k,0),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif
           if (TdiffMods) then
              if (k==(ks+ke)/2) then
                 write(70,410)x(i),y(j),0.d0,ymod(i,j,k,0)
                 write(70,410)x(i),y(j+1),0.d0,-ymod(i,j,k,1)
              endif
           endif
        endif
        if (vof_phase(i,j,k+1)==1) then
           if (height(i,j,k,6)>lim_cut .and. height(i,j,k,6)<1.0d0-lim_cut) then
              zmod(i,j,k,0)=height(i,j,k,6)*dh
              zmod(i,j,k,1)=dh-zmod(i,j,k,0)
           else if (-height(i,j,k,5)>lim_cut .and. -height(i,j,k,5)<1.0d0-lim_cut) then
              zmod(i,j,k,0)=-height(i,j,k,5)*dh
              zmod(i,j,k,1)=dh-zmod(i,j,k,0)
           else
              if (TdiffMods) stag_calls=stag_calls+1
              call staggered_cut_gen(cstagz,i,j,k,zmod(i,j,k,0),3,lim_cut)
              zmod(i,j,k,1)=MIN(MAX(dh-zmod(i,j,k,0),lim_cut*dh),(1.0d0-lim_cut)*dh)
           endif
           if (TdiffMods) then
              if (j==(js+je)/2) then
                 write(71,410)x(i),z(k),0.d0,zmod(i,j,k,0)
                 write(71,410)x(i),z(k+1),0.d0,-zmod(i,j,k,1) 
              endif
           endif
        endif
     endif ! Vapour cell

  enddo; enddo; enddo
  if (TdiffMods) then
     write(*,'("Number of staggered cut calls for mods: ",I8)')stag_calls
     close(70); close(71)
  endif
410 format(4e14.5)
  call ghost_x(xmod(:,:,:,0),1,reqd( 1: 4));   call ghost_x(xmod(:,:,:,1),1,reqd( 5: 8))
  call ghost_y(ymod(:,:,:,0),1,reqd( 9: 12));  call ghost_y(ymod(:,:,:,1),1,reqd( 13: 16))
  call ghost_z(zmod(:,:,:,0),1,reqd( 17: 20)); call ghost_z(zmod(:,:,:,1),1,reqd( 21: 24))
  call MPI_WAITALL(24,reqd,stat,ierr)

  do phase=0,1

     ! set initial phase temps for guess values 
     do k=ks-1,ke+1; do j=js-1,je+1; do i=is-1,ie+1
        if (vof_phase(i,j,k)==phase) then
           Tpass(i,j,k)=energy(i,j,k)/rho_cp(i,j,k) ! Setting initial guess value
        else
           Tpass(i,j,k)=T_sat
        endif
     enddo;enddo;enddo

     !set mod BC's for Dirichlet condition
     if (coords(1)==0) then
        if (BDRY_T(1)==0) xmod(is-1,jmin:jmax,kmin:kmax,phase)=0.5d0*dh
     endif

     !BC on x+
     if (coords(1)==nPx-1) then
        if (BDRY_T(4)==0) xmod(ie,jmin:jmax,kmin:kmax,phase)=0.5d0*dh
     endif

     !BC on y-
     if (coords(2)==0) then
        if (BDRY_T(2)==0) ymod(imin:imax,js-1,kmin:kmax,phase)=0.5d0*dh
     endif

     !BC on y+
     if (coords(2)==nPy-1) then     
        if (BDRY_T(5)==0) ymod(imin:imax,je,kmin:kmax,phase)=0.5d0*dh
     endif

     !BC on z-
     if (coords(3)==0) then
        if (BDRY_T(3)==0) zmod(imin:imax,jmin:jmax,ks-1,phase)=0.5d0*dh
     endif

     !BC on z+
     if (coords(3)==nPz-1) then
        if (BDRY_T(6)==0) zmod(imin:imax,jmin:jmax,ke,phase)=0.5d0*dh
     endif

     do k=ks,ke; do j=js,je; do i=is,ie
        if (.not.(vof_phase(i,j,k)==0 .or. vof_phase(i,j,k)==1)) &
             call pariserror('Phase error implicit temp diffusion') 
        if (vof_phase(i,j,k)==phase) then

           coeff(i,j,k,1) = 2.d0*kcc(phase)*Tumask(i-1,j,k)*dt &
                /((xmod(i-1,j,k,phase)+xmod(i,j,k,phase))*xmod(i-1,j,k,phase))
           coeff(i,j,k,2) = 2.d0*kcc(phase)*Tumask(i,j,k)*dt &
                /((xmod(i,j,k,phase)+xmod(i-1,j,k,phase))*xmod(i,j,k,phase))
           coeff(i,j,k,3) = 2.d0*kcc(phase)*Tvmask(i,j-1,k)*dt &
                /((ymod(i,j-1,k,phase)+ymod(i,j,k,phase))*ymod(i,j-1,k,phase))
           coeff(i,j,k,4) = 2.d0*kcc(phase)*Tvmask(i,j,k)*dt &
                /((ymod(i,j,k,phase)+ymod(i,j-1,k,phase))*ymod(i,j,k,phase))
           coeff(i,j,k,5) = 2.d0*kcc(phase)*Twmask(i,j,k-1)*dt &
                /((zmod(i,j,k-1,phase)+zmod(i,j,k,phase))*zmod(i,j,k-1,phase))
           coeff(i,j,k,6) = 2.d0*kcc(phase)*Twmask(i,j,k)*dt &
                /((zmod(i,j,k,phase)+zmod(i,j,k-1,phase))*zmod(i,j,k,phase))
           coeff(i,j,k,7) = rho_cp(i,j,k) + sum(coeff(i,j,k,1:6))
           coeff(i,j,k,8) = energy(i,j,k)
        else
           coeff(i,j,k,1:6) = 0.0d0
           coeff(i,j,k,7) = 1.0d0
           coeff(i,j,k,8) = T_sat
        endif
     enddo;enddo;enddo

     !Set Poisson BC's
     !BC on x-
     if (coords(1)==0) then
        do j=js,je; do k=ks,ke
           if (vof_phase(is,j,k)==phase) then
              if (BDRY_T(1)==0) then !Dirichlet
                 coeff(is,j,k,7)=coeff(is,j,k,7)+2.d0*kcc(phase)*dt&
                      /((xmod(is-1,j,k,phase)+xmod(is,j,k,phase))*xmod(is-1,j,k,phase))
                 coeff(is,j,k,8)=coeff(is,j,k,8)+2.d0*kcc(phase)*dt*BC_T(1)&
                      /((xmod(is-1,j,k,phase)+xmod(is,j,k,phase))*xmod(is-1,j,k,phase))
              endif
              if (BDRY_T(1)==1) then !Neumann
                 coeff(is,j,k,8)=coeff(is,j,k,8)-2.d0*kcc(phase)*dt*BC_T(1)&
                      /(xmod(is-1,j,k,phase)+xmod(is,j,k,phase))
              endif
           endif
        enddo; enddo
     endif

     !BC on x+
     if (coords(1)==nPx-1) then
        do j=js,je; do k=ks,ke
           if (vof_phase(ie,j,k)==phase) then
              if (BDRY_T(4)==0) then 
                 coeff(ie,j,k,7)=coeff(ie,j,k,7)+2.d0*kcc(phase)*dt&
                      /((xmod(ie,j,k,phase)+xmod(ie-1,j,k,phase))*xmod(ie,j,k,phase)) 
                 coeff(ie,j,k,8)=coeff(ie,j,k,8)+2.d0*kcc(phase)*dt*BC_T(4)&
                      /((xmod(ie,j,k,phase)+xmod(ie-1,j,k,phase))*xmod(ie,j,k,phase))
              endif
              if (BDRY_T(4)==1) then !Neumann
                 coeff(ie,j,k,8)=coeff(ie,j,k,8)+2.d0*kcc(phase)*dt*BC_T(4)&
                      /(xmod(ie,j,k,phase)+xmod(ie-1,j,k,phase))
              endif
           endif
        enddo; enddo
     endif

     !BC on y-
     if (coords(2)==0) then
        do i=is,ie; do k=ks,ke
           if (vof_phase(i,js,k)==phase) then
              if (BDRY_T(2)==0) then !Dirichlet
                 coeff(i,js,k,7)=coeff(i,js,k,7)+2.d0*kcc(phase)*dt&
                      /((ymod(i,js-1,k,phase)+ymod(i,js,k,phase))*ymod(i,js-1,k,phase))
                 coeff(i,js,k,8)=coeff(i,js,k,8)+2.d0*kcc(phase)*dt*BC_T(2)&
                      /((ymod(i,js-1,k,phase)+ymod(i,js,k,phase))*ymod(i,js-1,k,phase))
              endif
              if (BDRY_T(2)==1) then !Neumann
                 coeff(i,js,k,8)=coeff(i,js,k,8)-2.d0*kcc(phase)*dt*BC_T(2)&
                      /(ymod(i,js-1,k,phase)+ymod(i,js,k,phase))
              endif
           endif
        enddo; enddo
     endif

     !BC on y+
     if (coords(2)==nPy-1) then
        do i=is,ie; do k=ks,ke
           if (vof_phase(i,je,k)==phase) then
              if (BDRY_T(5)==0) then 
                 coeff(i,je,k,7)=coeff(i,je,k,7)+2.d0*kcc(phase)*dt&
                      /((ymod(i,je,k,phase)+ymod(i,je-1,k,phase))*ymod(i,je,k,phase))
                 coeff(i,je,k,8)=coeff(i,je,k,8)+2.d0*kcc(phase)*dt*BC_T(5)&
                      /((ymod(i,je,k,phase)+ymod(i,je-1,k,phase))*ymod(i,je,k,phase))
              endif
              if (BDRY_T(5)==1) then !Neumann 
                 coeff(i,je,k,8)=coeff(i,je,k,8)+2.d0*kcc(phase)*dt*BC_T(5)&
                      /(ymod(i,je,k,phase)+ymod(i,je-1,k,phase))
              endif
           endif
        enddo; enddo
     endif

     !BC on z-
     if (coords(3)==0) then
        do i=is,ie; do j=js,je
           if (vof_phase(i,j,ks)==phase) then
              if (BDRY_T(3)==0) then !Dirichlet
                 coeff(i,j,ks,7)=coeff(i,j,ks,7)+2.d0*kcc(phase)*dt&
                      /((zmod(i,j,ks-1,phase)+zmod(i,j,ks,phase))*zmod(i,j,ks-1,phase))
                 coeff(i,j,ks,8)=coeff(i,j,ks,8)+2.d0*kcc(phase)*dt*BC_T(3)&
                      /((zmod(i,j,ks-1,phase)+zmod(i,j,ks,phase))*zmod(i,j,ks-1,phase))
              endif
              if (BDRY_T(3)==1) then !Neumann
                 coeff(i,j,ks,8)=coeff(i,j,ks,8)-2.d0*kcc(phase)*dt*BC_T(3)&
                      /(zmod(i,j,ks-1,phase)+zmod(i,j,ks,phase))
              endif
           endif
        enddo; enddo
     endif

     !BC on z+
     if (coords(3)==nPz-1) then
        do i=is,ie; do j=js,je
           if (vof_phase(i,j,ke)==phase) then
              if (BDRY_T(6)==0) then 
                 coeff(i,j,ke,7)=coeff(i,j,ke,7)+2.d0*kcc(phase)*dt&
                      /((zmod(i,j,ke,phase)+zmod(i,j,ke-1,phase))*zmod(i,j,ke,phase))
                 coeff(i,j,ke,8)=coeff(i,j,ke,8)+2.d0*kcc(phase)*dt*BC_T(6)&
                      /((zmod(i,j,ke,phase)+zmod(i,j,ke-1,phase))*zmod(i,j,ke,phase))
              endif
              if (BDRY_T(6)==1) then !Neumann  
                 coeff(i,j,ke,8)=coeff(i,j,ke,8)+2.d0*kcc(phase)*dt*BC_T(6)&
                      /(zmod(i,j,ke,phase)+zmod(i,j,ke-1,phase))
              endif
           endif
        enddo; enddo
     endif

     !Solve for phase temps
     if (HYPRE_Timp) then
        call poi_solve(coeff,Tpass,maxError/MaxDt*ErrorScaleHYPRE,maxit,iters,HYPRESolverType)
        call do_all_ghost(Tpass)
     else     
        call NewSolver(coeff,Tpass,maxError/MaxDt,beta,maxit,iters,ierr,ResNormOrderPressure)
     endif
     if (iters>=maxit) then   
        if (.not.(HYPRE_Timp)) then
           ! Get the init values again
           if(rank==0) write(*,'(" Re-solving T with MG ")')
           do k=ks,ke; do j=js,je; do i=is,ie
              if (vof_phase(i,j,k)==phase) then
                 Tpass(i,j,k)=energy(i,j,k)/rho_cp(i,j,k) ! Setting initial guess value
              else
                 Tpass(i,j,k)=T_sat
              endif
           enddo;enddo;enddo   
           call poi_solve(coeff,Tpass,maxError/MaxDt*ErrorScaleHYPRE,maxit,iters,HYPRESolverType)
           call do_all_ghost(Tpass)
        endif
        call calcResiduals(coeff,Tpass,res_tabTimp)
        if (ResNormOrderPressure>2) then
           residual=res_tabTimp(3)
        else
           residual=res_tabTimp(ResNormOrderPressure)
        endif
        if (residual*MaxDt/maxError > DivTolTimp) &
             call pariserror("Solver diverge in implicit T solve!")
     endif

     if(mod(itimestep,termout)==0) then
        call calcResiduals(coeff,Tpass,res_tabTimp)
        if(rank==0) then
           write(*,'("  Implicit T solve phase :   ",I8)')phase
           write(*,'("  Implicit T residuals*Maxdt L1:",e8.1,"         L2:",e8.1,"       Linf:",e8.1)') &
                res_tabTimp*MaxDt
           write(*,'("  Implicit T iterations     :",I8,  " tolerance :",e8.1," norm order:",I4)')   &
                iters,maxError,ResNormOrderPressure        
           write(*,*)
        endif
     endif

     do k=ks,ke; do j=js,je; do i=is,ie
        if (vof_phase(i,j,k)==phase) then
           Te(i,j,k)=Tpass(i,j,k)
        endif
     enddo;enddo;enddo
  enddo

end subroutine TempDiffusionImplicit

!Energy advection
subroutine EnergyAdvection
  use module_grid
  use module_boil
  use module_tmpvar
  use module_BC
  implicit none
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: face_min,face_plus
  integer :: d,i,j,k,i0,j0,k0
  logical :: weno_out
  work=0.d0

  do d=1,3
     call init_i0j0k0 (d,i0,j0,k0)

     do k=ks,ke
        do j=js,je
           do i=is,ie         

              if ((j==js+ny/2) .and. (k==3) .and. (i>=61) .and. (i<=64) .and. (d==1)) then
                 weno_out=.true.
              else
                 weno_out=.false.
              endif

              face_min(i,j,k)=interpoleWENO5(energy(i+i0*2,j+j0*2,k+k0*2),energy(i+i0,j+j0,k+k0),&
                   energy(i,j,k),energy(i-i0,j-j0,k-k0),energy(i-i0*2,j-j0*2,k-k0*2),weno_out)
              face_plus(i,j,k)=interpoleWENO5(energy(i-i0*2,j-j0*2,k-k0*2),energy(i-i0,j-j0,k-k0),&
                   energy(i,j,k),energy(i+i0,j+j0,k+k0),energy(i+i0*2,j+j0*2,k+k0*2),weno_out)
           enddo
        enddo
     enddo

     !comms? Yes...
     call do_all_ghost(face_min)
     call do_all_ghost(face_plus)

     call EnergyConvection_Onedim(d)

  enddo

  do i=is,ie; do j=js,je;  do k=ks,ke
     energy(i,j,k)=energy(i,j,k) - ( (work(i,j,k,1)-work(i-1,j,k,1)) + &
          (work(i,j,k,2)-work(i,j-1,k,2)) + (work(i,j,k,3)-work(i,j,k-1,3)) ) 
  enddo; enddo; enddo

contains
  subroutine  EnergyConvection_Onedim(d)
    use module_grid
    use module_flow
    implicit none
    real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: us
    integer, intent(in) :: d
    integer :: i,j,k,is0,js0,ks0,ie0,je0,ke0
    real(8) :: T_bdry,dh
    dh=dx(is) !uniform grids only


    !Set boundary values
    is0=is; js0=js; ks0=ks
    ie0=ie; je0=je; ke0=ke

    select case (d)
    case(1)
       ie0=ieu
       us=u

       if(coords(1)==0) then
          do j=jmin,jmax
             do k=kmin,kmax
                if (BDRY_T(1)==0) then
                   T_bdry=BC_T(1)
                else
                   T_bdry=Te(is,j,k)-BC_T(1)*dh/2.d0
                endif
                if (u(is-1,j,k)>=0) then
                   work(is-1,j,k,1)=rho_cp(is-1,j,k)*T_bdry*u(is-1,j,k)*dt/dh
                else
                   work(is-1,j,k,1)=face_min(is,j,k)*u(is-1,j,k)*dt/dh !fixme or improve here
                endif
             enddo
          enddo
       endif

       if (coords(1)==nPx-1) then
          do j=jmin,jmax
             do k=kmin,kmax
                if (BDRY_T(4)==0) then 
                   T_bdry=BC_T(4)
                else
                   T_bdry=Te(ie,j,k)+BC_T(4)*dh/2.d0
                endif
                if (u(ie,j,k)>=0) then
                   work(ie,j,k,1)=face_plus(ie,j,k)*u(ie,j,k)*dt/dh !fixme or improve here
                else
                   work(ie,j,k,1)=rho_cp(ie+1,j,k)*T_bdry*u(ie,j,k)*dt/dh
                endif
             enddo
          enddo
       endif

    case(2)
       je0=jev
       us=v
       
       if(coords(2)==0) then
          do i=imin,imax
             do k=kmin,kmax
                if (BDRY_T(2)==0) then
                   T_bdry=BC_T(2)
                else
                   T_bdry=Te(i,js,k)-BC_T(2)*dh/2.d0
                endif
                if (v(i,js-1,k)>=0) then
                   work(i,js-1,k,2)=rho_cp(i,js-1,k)*T_bdry*v(i,js-2,k)*dt/dh
                else
                   work(i,js-1,k,2)=face_min(i,js,k)*v(i,js-1,k)*dt/dh 
                endif
             enddo
          enddo
       endif

       if (coords(2)==nPy-1) then
          do i=imin,imax
             do k=kmin,kmax
                if (BDRY_T(5)==0) then 
                   T_bdry=BC_T(5)
                else
                   T_bdry=Te(i,je,k)+BC_T(5)*dh/2.d0
                endif
                if (v(i,je,k)>=0) then
                   work(i,je,k,2)=face_plus(i,je,k)*v(i,je,k)*dt/dh
                else
                   work(i,je,k,2)=rho_cp(i,je+1,k)*T_bdry*v(i,je,k)*dt/dh
                endif
             enddo
          enddo
       endif

    case(3)
       ke0=kew
       us=w

       if(coords(3)==0) then
          do i=imin,imax
             do j=jmin,jmax
                if (BDRY_T(3)==0) then
                   T_bdry=BC_T(3)
                else
                   T_bdry=Te(i,j,ks)-BC_T(3)*dh/2.d0
                endif
                if (w(i,j,ks-1)>=0) then
                   work(i,j,ks-1,3)=rho_cp(i,j,ks-1)*T_bdry*w(i,j,ks-1)*dt/dh
                else
                   work(i,j,ks-1,3)=face_min(i,j,ks)*w(i,j,ks-1)*dt/dh
                endif
             enddo
          enddo
       endif

       if (coords(3)==nPz-1) then
          do i=imin,imax
             do j=jmin,jmax
                if (BDRY_T(6)==0) then 
                   T_bdry=BC_T(6)
                else
                   T_bdry=Te(i,j,ke)+BC_T(6)*dh/2.d0
                endif
                if (w(i,j,ke)>=0) then
                   work(i,j,ke,3)=face_plus(i,j,ke)*w(i,j,ke)*dt/dh
                else
                   work(i,j,ke,3)=rho_cp(i,j,ke+1)*T_bdry*w(i,j,ke)*dt/dh
                endif
             enddo
          enddo
       endif

    end select

    do k=ks0,ke0; do j=js0,je0; do i=is0,ie0
       if (us(i,j,k)>=0.d0) then
          work(i,j,k,d)=face_plus(i,j,k)*us(i,j,k)*dt/dh
       else
          work(i,j,k,d)=face_min(i+i0,j+j0,k+k0)*us(i,j,k)*dt/dh
       endif
    enddo;enddo;enddo

  end subroutine EnergyConvection_Onedim
end subroutine EnergyAdvection

subroutine get_heat_source(rho1,rho2,tstep)
  use module_grid
  use module_VOF
  use module_boil
  use module_BC
  use module_2phase
  implicit none
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax,0:1) :: ngradT
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: dint, n1, n2, n3
  real(8), intent(in) :: rho1, rho2
  integer, intent(in) :: tstep
  integer :: i,j,k,phase,i0,j0,k0,i1,j1,k1,phase_sign
  real(8) :: m(3), stencil3x3(-1:1,-1:1,-1:1), nr(3), crdshift(1:3), npos(3)
  real(8) :: alpha,al3d,magm,alpha_d,alpha_d_shift, dh, magpos, col_fact, maxcol
  real(8) :: hccf, drel, drel_t, nbr_count, phase_heat(0:1), mdotc
 
  ! Variables for tests and debugging
  real(8) :: e_gradnT, gradnT, maxe_m, L1_m, L2_m, mcount, relerr_m
  real(8) :: e_d ,d_theory, maxe_d, L1_d, L2_d, count, relerr_d

  mdot=0.0d0; s_v=0.0d0
  dh=dx(is) !Assume unstretched grid
  n1=0.0d0; n2=0.0d0; n3=0.0d0
  if (type==2) then
     ! Need to get all normals, since we are taking 2nd neighbours for temp gradients
     do k=ks,ke
        do j=js,je
           do i=is,ie
              do i0=-1,1; do j0=-1,1; do k0=-1,1
                 stencil3x3(i0,j0,k0) = cvof(i+i0,j+j0,k+k0)
              enddo;enddo;enddo
              call mycs(stencil3x3,nr)
              n1(i,j,k) = nr(1)
              n2(i,j,k) = nr(2)
              n3(i,j,k) = nr(3)
           enddo
        enddo
     enddo
     call do_ghost_vector(n1,n2,n3)
     !BCs n??
     dint=2.d20
     ngradT=0.d0

     if (ngradT_test) then
        OPEN(unit=88,file='d_detail.txt',position='append')
        OPEN(unit=89,file='d_stats.txt',position='append')
        WRITE(89,'("Origin of 3D bubble",3e14.5," radius ",e14.5," grid res NX, NY, NZ: ",3I8)')&
             xc(1),yc(1),zc(1),rad(1),Nx,Ny,Nz 
        WRITE(89,'("     Phase      Cell count       L1            L2         L_inf")')
                   !1234567891123412345678911234123456789112341234567891123412345678911234 
        OPEN(unit=90,file='m_detail.txt',position='append')
        OPEN(unit=91,file='m_stats.txt',position='append')
     endif

     call phase_normal_gradient(0)

     ! Set boundary values for ngradT and dint
     call SetBoundaryHeatFlux
     call do_all_ghost(ngradT(:,:,:,0))
     call do_all_ghost(ngradT(:,:,:,1))
     call do_all_ghost(dint)

     if (ngradT_test) then
        mcount=0.d0
        L1_m=0.d0; L2_m=0.d0; maxe_m=0.0d0
     endif

     !Cycle cut cells, get avg gradients from both sides to int
     do i=is,ie; do j=js,je; do k=ks,ke
        if (topo_mask(i,j,k)==0) then
           !loop phase
           do phase=0,1
              phase_sign=2*phase-1
              phase_heat(phase)=0.d0
              if (phase==0) then
                 hccf=kc1/h_fg
              else
                 hccf=kc2/h_fg
              endif

              nbr_count=0.d0
              drel_t=0.d0
              do k0=-2,2; do j0=-2,2; do i0=-2,2
                 if ( vof_flag(i+i0,j+j0,k+k0)==phase .and. &
                      (topo_mask(i+i0,j+j0,k+k0)==phase_sign .or. topo_mask(i+i0,j+j0,k+k0)==2*phase_sign) ) then
                    if (i0==0 .and. j0==0 .and. k0==0) call pariserror('Topological error')
                    m(1)=n1(i,j,k); m(2)=n2(i,j,k); m(3)=n3(i,j,k)  
                    magm=sqrt(m(1)**2.0d0+m(2)**2.0d0+m(3)**2.0d0)
                    magpos=sqrt( (x(i)-x(i+i0))**2.d0 + (y(j)-y(j+j0))**2.d0 + (z(k)-z(k+k0))**2.d0 )
                    npos(1)=(x(i+i0)-x(i))/magpos
                    npos(2)=(y(j+j0)-y(j))/magpos
                    npos(3)=(z(k+k0)-z(k))/magpos
                    drel = ABS(DOT_PRODUCT(npos,m/magm))/(magpos**2.0d0)
                    if (dint(i+i0,j+j0,k+k0)<dh*3.0d0*sqrt(2.0d0)) then
                       drel_t=drel_t + drel
                       phase_heat(phase)=phase_heat(phase)+hccf*drel*ngradT(i+i0,j+j0,k+k0,phase)
                       nbr_count=nbr_count+1.d0
                    endif
                 endif
              enddo; enddo; enddo !neighbour loop
              if (nbr_count>0.9d0) then
                 phase_heat(phase)=phase_heat(phase)/drel_t              
              else
                 write(*,'("Topological error, no bulk neighbours to mixed cell")')
                 write(*,'("CVoF, flag, phase, i,j,k: ",e14.5,5I5)')cvof(i,j,k),vof_flag(i,j,k),phase,i,j,k
                 write(*,'("Topo masks -2:+2 matrix: ")')
                 do j0=-2,2
                    write(*,'("j0: ",I5)')j0
                    do k0=-2,2
                       write(*,'("k0: ",6I5)')k0,topo_mask(i-2:i+2,j+j0,k+k0)
                    enddo
                 enddo
                 write(*,'("dint -2:+2 matrix: ")')
                 do j0=-2,2
                    write(*,'("j0: ",I5)')j0
                    do k0=-2,2
                       write(*,'("k0: ",I5,5e14.5)')k0,dint(i-2:i+2,j+j0,k+k0)
                    enddo
                 enddo
                 call pariserror('Topological error, no bulk neighbours to mixed cell')
              endif
           enddo ! phase
           
           mdot(i,j,k)=sum(phase_heat(0:1))
           if (ngradT_test) then
              mcount=mcount+1   
              if (test_droplet) then
                 gradnT=T_sat/rad(1)
              else if (test_plane) then
                 call pariserror('No plane test yet...drop is tougher test than plane.')
              else
                 call pariserror('WRONG TEST TYPE')
              endif
              e_gradnT=ABS(mdot(i,j,k)) ! Theoretical mdot=0 for this setup!
              relerr_m=ABS(e_gradnT)/gradnT
              L1_m=L1_m+relerr_m
              L2_m=L2_m+relerr_m**2.0d0
              maxe_m = MAX(relerr_m, maxe_m)
              write(90,'(5e14.5)')x(i),y(j),z(k),mdot(i,j,k),relerr_m
           endif
           s_v(i,j,k)=mdot(i,j,k)*inv_drho*cell_area(i,j,k)/dh
        endif ! mixed cells
     enddo;enddo;enddo
     if (ngradT_test) then
        if (mcount>0.5d0) then
           L1_m=L1_m/mcount
           L2_m=sqrt(L2_m/mcount)
           write(91,'(4e14.5)')count,L1_m,L2_m,maxe_m
        else
           write(91,'("Error: no cut cells found ")')
        endif
     endif
  else ! type 1, constant mass transfer rate per unit area
     mdotc=q_set/h_fg
     do k=ks,ke; do j=js,je; do i=is,ie
        if (topo_mask(i,j,k)==0) then
           mdot(i,j,k) = q_set/h_fg
           s_v(i,j,k) = mdot(i,j,k)*inv_drho*cell_area(i,j,k)/dh
        endif
     enddo; enddo; enddo
  endif
  
  if (ngradT_test) then
     s_v=0.0d0
     mdot=0.0d0
     close(88)
     close(89)
     close(90)
     close(91)
  endif
contains
  subroutine SetBoundaryHeatFlux
    implicit none

    if( coords(1)==0 ) then
       if (vofbdry_cond(1) == '90deg' .and. bdry_cond(1)==2) then
          dint(is-1,:,:) = dint(is,:,:)
          dint(is-2,:,:) = dint(is+1,:,:)
          ngradT(is-1,:,:,0:1) = ngradT(is,:,:,0:1)
          ngradT(is-2,:,:,0:1) = ngradT(is+1,:,:,0:1)
       else
          if (BDRY_T(1)<2) call phase_normal_gradient(1)
       endif
    endif

    if( coords(2)==0 ) then
       if (vofbdry_cond(2) == '90deg' .and. bdry_cond(2)==2) then
          dint(:,js-1,:) = dint(:,js  ,:)
          dint(:,js-2,:) = dint(:,js+1,:)
          ngradT(:,js-1,:,0:1) = ngradT(:,js,:,0:1)
          ngradT(:,js-2,:,0:1) = ngradT(:,js+1,:,0:1)
       else
          if (BDRY_T(2)<2) call phase_normal_gradient(2)
       endif
    endif

    if( coords(3)==0 ) then
       if (vofbdry_cond(3) == '90deg' .and. bdry_cond(3)==2) then
          dint(:,:,ks-1) = dint(:,:,ks)
          dint(:,:,ks-2) = dint(:,:,ks+1)
          ngradT(:,:,ks-1,0:1) = ngradT(:,:,ks,0:1)
          ngradT(:,:,ks-2,0:1) = ngradT(:,:,ks+1,0:1)
       else
          if (BDRY_T(3)<2) call phase_normal_gradient(3)
       endif
    endif

    if( coords(1)==nPx-1 ) then
       if (vofbdry_cond(4) == '90deg' .and. bdry_cond(4)==2) then
          dint(ie+1,:,:) = dint(ie,:,:)
          dint(ie+2,:,:) = dint(ie-1,:,:)
          ngradT(ie+1,:,:,0:1) = ngradT(ie,:,:,0:1)
          ngradT(ie+2,:,:,0:1) = ngradT(ie-1,:,:,0:1)
       else
          if (BDRY_T(4)<2) call phase_normal_gradient(4)
       endif
    endif

    if( coords(2)==nPy-1 ) then
       if (vofbdry_cond(5) == '90deg' .and. bdry_cond(5)==2) then
          dint(:,je+1,:) = dint(:,je  ,:)
          dint(:,je+2,:) = dint(:,je-1,:)
          ngradT(:,je+1,:,0:1) = ngradT(:,je  ,:,0:1)
          ngradT(:,je+2,:,0:1) = ngradT(:,je-1,:,0:1)
       else
          if (BDRY_T(5)<2) call phase_normal_gradient(5)
       endif
    endif

    if( coords(3)==nPz-1 ) then
       if (vofbdry_cond(6) == '90deg' .and. bdry_cond(6)==2) then
          dint(:,:,ke+1) = dint(:,:,ke)
          dint(:,:,ke+2) = dint(:,:,ke-1)
          ngradT(:,:,ke+1,0:1) = ngradT(:,:,ke  ,0:1)
          ngradT(:,:,ke+2,0:1) = ngradT(:,:,ke-1,0:1)
       else
          if (BDRY_T(6)<2) call phase_normal_gradient(6)
       endif
    endif

  end subroutine SetBoundaryHeatFlux
  
  subroutine phase_normal_gradient(face)
    use module_VOF
    implicit none
    integer, intent(in) :: face ! 0: Bulk 1: x- 2: y- ... 5: y+ 6: z+
    integer :: isd,ied,jsd,jed,ksd,ked,isn,ien,jsn,jen,ksn,ken
    real(8) :: al3d
    
    !case(0) is the bulk case
    isd=is; jsd=js; ksd=ks; ied=ie; jed=je; ked=ke
    isn=-2; jsn=-2; ksn=-2; ien=2; jen=2; ken=2
    select case(face)
    case(1)
       isd=imin; ied=is-1
       isn=1; ien=2
    case(2)
       jsd=jmin; jed=js-1
       jsn=1; jen=2
    case(3)
       ksd=kmin; ked=ks-1
       ksn=1; ken=2
    case(4)
       isd=ie+1; ied=imax
       isn=-2; ien=-1
    case(5)
       jsd=je+1; jed=jmax
       jsn=-2; jen=-1
    case(6)
       ksd=ke+1; ken=kmax
       ksn=-2; ken=-1
    end select

    do phase=0,1
        phase_sign=2*phase-1
        
        if (ngradT_test) then
           count=0.d0
           L1_d=0.d0; L2_d=0.d0; maxe_d=0.0d0
        endif
        
        do i=isd,ied; do j=jsd,jed; do k=ksd,ked
           if (vof_flag(i,j,k)==phase .and.  &
                (topo_mask(i,j,k)==phase_sign .or. topo_mask(i,j,k)==2*phase_sign) ) then
              maxcol = 0.d0
              do k0=ksn,ken; do j0=jsn,jen; do i0=isn,ien               
                 if (topo_mask(i+i0,j+j0,k+k0)==0) then
                    m(1)=n1(i+i0,j+j0,k+k0); m(2)=n2(i+i0,j+j0,k+k0); m(3)=n3(i+i0,j+j0,k+k0)  
                    magm=sqrt(m(1)**2.0d0+m(2)**2.0d0+m(3)**2.0d0)
                    magpos=sqrt( (x(i)-x(i+i0))**2.d0 + (y(j)-y(j+j0))**2.d0 + (z(k)-z(k+k0))**2.d0 )
                    npos(1)=(x(i)-x(i+i0))/magpos
                    npos(2)=(y(j)-y(j+j0))/magpos
                    npos(3)=(z(k)-z(k+k0))/magpos
                    col_fact=ABS(DOT_PRODUCT(npos,m/magm))
                    if (col_fact>maxcol) then
                       maxcol=col_fact
                       alpha=al3D(m,cvof(i+i0,j+j0,k+k0))
                       alpha_d=alpha/magm !rescale alpha for normal distance
                       !change coordinates
                       crdshift(1)=(x(i)-xh(i+i0)+dh)/dh
                       crdshift(2)=(y(j)-yh(j+j0)+dh)/dh
                       crdshift(3)=(z(k)-zh(k+k0)+dh)/dh
                       alpha_d_shift=alpha_d-(m(1)/magm*crdshift(1)+m(2)/magm*crdshift(2)+m(3)/magm*crdshift(3))
                       !calc normal distance to neighbour node
                       dint(i,j,k)=ABS(alpha_d_shift*dh)
                    endif ! largest collinearity of rel cell pos and plane normal vectors
                 endif ! interface cells
              enddo; enddo; enddo !neighbour loop

              if (dint(i,j,k)<dh*3.0d0*sqrt(2.0d0)) ngradT(i,j,k,phase)=(Te(i,j,k)-T_sat)/dint(i,j,k)

              if (ngradT_test) then
                 count=count+1   
                 if (test_droplet) then
                    d_theory= ABS( sqrt( (x(i)-xc(1))**2.d0 + (y(j)-yc(1))**2.d0 + (z(k)-zc(1))**2.d0) - rad(1))
                 else if (test_plane) then
                    d_theory=ABS(y(j)-plane) !!KYKHIER
                    call pariserror('No plane test yet...drop is tougher test than plane.')
                 else
                    call pariserror('WRONG TEST TYPE')
                 endif
                 e_d=d_theory-dint(i,j,k)
                 relerr_d=ABS(e_d)/rad(1)
                 L1_d=L1_d+relerr_d
                 L2_d=L2_d+relerr_d**2.0d0
                 maxe_d = MAX(relerr_d, maxe_d)
                 write(88,'(I7,7e14.5)')phase,x(i),y(j),z(k),col_fact,dint(i,j,k),d_theory,relerr_d
              endif

           endif ! 1st and 2nd neighbours
        enddo; enddo; enddo

        if (ngradT_test) then
           if (count>0.5d0) then
              L1_d=L1_d/count
              L2_d=sqrt(L2_d/count)
              write(89,'(I14,4e14.5)')phase,count,L1_d,L2_d,maxe_d
           else
              write(89,'("Error: no distance found in phase ",I8)')phase
           endif
        endif
     enddo ! phase loop
    
  end subroutine phase_normal_gradient
end subroutine get_heat_source

subroutine vofandenergysweeps(tswap)
  use module_grid
  use module_VOF
  use module_flow
  use module_tmpvar
  use module_boil
  use module_BC
  implicit none
  integer, intent(in) :: tswap
  integer :: i,j,k

  if (VOF_advect=='Dick_Yue') call c_mask(cvof, work(:,:,:,2))
  if (MOD(tswap,3).eq.0) then  ! do z then x then y 
     call swp_energy(3)
     call swp(w_l,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC

     call swp_energy(1)
     call swp(u_l,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC

     call swp_energy(2)
     call swp(v_l,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC
  elseif (MOD(tswap,3).eq.1) then ! do y z x
     call swp_energy(2)
     call swp(v_l,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC

     call swp_energy(3)
     call swp(w_l,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC
     
     call swp_energy(1)
     call swp(u_l,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC
  else ! do x y z
     call swp_energy(1)
     call swp(u_l,cvof,vof_flag,1,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC

     call swp_energy(2)
     call swp(v_l,cvof,vof_flag,2,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC

     call swp_energy(3)
     call swp(w_l,cvof,vof_flag,3,work(:,:,:,1),work(:,:,:,2),work(:,:,:,3))
     call ThermalEnBC
  endif
contains
  subroutine ThermalEnBC
    implicit none
    if (adv_method==1) then
       call linfunc(rho_cp,rho1*cp1,rho2*cp2,ArithMean)
       call get_Te_from_energy
       call SetTempBC
       call SetBoundaryEnergy
       call do_all_ghost(energy)
    endif
    call get_ghost_masks
  end subroutine ThermalEnBC
end subroutine vofandenergysweeps

subroutine swp_energy(d)
  use module_BC
  use module_boil
  use module_VOF
  use module_flow
  use module_tmpvar
  implicit none
  integer, intent(in) :: d 

  if (VOF_advect=='Dick_Yue') then
     select case (d)
     case(1)
        call swpr_energy(u,u_l,cvof,1,work(:,:,:,2))
     case(2)
        call swpr_energy(v,v_l,cvof,2,work(:,:,:,2))
     case(3)
        call swpr_energy(w,w_l,cvof,3,work(:,:,:,2))
     end select
  elseif (VOF_advect=='CIAM') then
     select case (d)
     case(1)
        call swpz_energy(u,u_l,cvof,1)
     case(2)
        call swpz_energy(v,v_l,cvof,2)
     case(3)
        call swpz_energy(w,w_l,cvof,3)
     end select
  else
     call pariserror("*** unknown vof scheme")
  endif

end subroutine swp_energy

subroutine swpr_energy(us,ul_adv,c,d,cg)
  use module_grid
  use module_boil
  use module_flow
  use module_BC
  implicit none
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us, ul_adv, c, cg
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: ene1,ene3,eni1,eni3,T_min,T_plus
  integer, intent(in) :: d
  real(8) :: rhocp_ll,rhocp_l, rhocp_r, rhocp_rr, rhocp_c
  real(8) :: T_avg, a2, a1, a2_l, a1_l, dh, T_l, T_r, vof, et_comp
  real(8) :: nr(3), stencil3x3(-1:1,-1:1,-1:1), alpha, deltax(3),x0(3),fl3d, al3d
  integer :: i0,j0,k0,i,j,k,i1,j1,k1

  real(8) :: de_adv

  call init_i0j0k0 (d,i0,j0,k0)
  if(ng.lt.2) call pariserror("wrong ng")
  dh = dxh(is) ! Only uniform grids!

  if (EtAdvScheme=='WENO5') then
     do k=ks,ke
        do j=js,je
           do i=is,ie         
              rhocp_ll=(rho2*cp2*c(i-2*i0,j-2*j0,k-2*k0)+rho1*cp1*(1.d0-c(i-2*i0,j-2*j0,k-2*k0)))
              rhocp_l =(rho2*cp2*c(i-  i0,j-  j0,k-  k0)+rho1*cp1*(1.d0-c(i-  i0,j-  j0,k-  k0)))
              rhocp_c =(rho2*cp2*c(i     ,j     ,k     )+rho1*cp1*(1.d0-c(i     ,j     ,k     )))
              rhocp_r =(rho2*cp2*c(i+  i0,j+  j0,k+k0  )+rho1*cp1*(1.d0-c(i+i0  ,j+j0  ,k+k0  )))
              rhocp_rr=(rho2*cp2*c(i+2*i0,j+2*j0,k+2*k0)+rho1*cp1*(1.d0-c(i+2*i0,j+2*j0,k+2*k0)))

              T_avg = energy(i,j,k)/(rhocp_c)

              T_min(i,j,k)=interpoleWENO5(energy(i+i0*2,j+j0*2,k+k0*2)/(rhocp_rr),energy(i+i0,j+j0,k+k0)/(rhocp_r),&
                   T_avg,energy(i-i0,j-j0,k-k0)/(rhocp_l),energy(i-i0*2,j-j0*2,k-k0*2)/(rhocp_ll),.false.)
              T_plus(i,j,k)=interpoleWENO5(energy(i-i0*2,j-j0*2,k-k0*2)/(rhocp_ll),&
                   energy(i-i0,j-j0,k-k0)/(rhocp_l),T_avg,energy(i+i0,j+j0,k+k0)/(rhocp_r),&
                   energy(i+i0*2,j+j0*2,k+k0*2)/(rhocp_rr),.false.)
           enddo
        enddo
     enddo

     call SetBoundFaceTemps(d,T_min,T_plus)
     call do_all_ghost(T_min)
     call do_all_ghost(T_plus)
  endif

  ene1=0.0d0; ene3=0.0d0
  eni1=0.0d0; eni3=0.0d0
  !Fluxes: assume VOF and energy BC's allow for s-1 and e+1 fluxes to be calculated correctly
  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           a2 = us(i,j,k)*dt/dh
           a2_l = ul_adv(i,j,k)*dt/dh
           a1 = us(i-i0,j-j0,k-k0)*dt/dh
           a1_l = ul_adv(i-i0,j-j0,k-k0)*dt/dh

           if (EtAdvScheme=='WENO5') then
              if (us(i,j,k)>0.d0) then
                 T_r=T_plus(i,j,k) 
              else
                 T_r=T_min(i+i0,j+j0,k+k0)
              endif

              if (us(i-i0,j-i0,k-k0)>0.d0) then
                 T_l=T_plus(i-i0,j-j0,k-k0) 
              else
                 T_l=T_min(i,j,k)
              endif
           else
              rhocp_l =(rho2*cp2*c(i-i0,j-j0,k-k0)+rho1*cp1*(1.d0-c(i-i0,j-j0,k-k0)))
              rhocp_c =(rho2*cp2*c(i   ,j   ,k   )+rho1*cp1*(1.d0-c(i   ,j   ,k   )))
              rhocp_r =(rho2*cp2*c(i+i0,j+j0,k+k0)+rho1*cp1*(1.d0-c(i+i0,j+j0,k+k0)))

              T_l = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                   energy(i+i0,j+j0,k+k0)/rhocp_r,EtAdvScheme,-0.5d0*(1.d0+a1))           
              T_r = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                   energy(i+i0,j+j0,k+k0)/rhocp_r,EtAdvScheme,0.5d0*(1.d0-a2))
           endif
           
           if (topo_mask(i,j,k)==0) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                 stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call mycs(stencil3x3,nr)
              alpha = al3d(nr,c(i,j,k))
              x0=0d0
              deltax=1d0

              if (a1_l<0.d0) then
                 deltax(d)=-a1_l
                 if (.not.(EtAdvScheme=='WENO5')) then
                    T_l = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                         energy(i+i0,j+j0,k+k0)/rhocp_r,AdvectionScheme,-0.5d0*(1.d0+a1_l))
                 endif
                 vof = fl3d(nr,alpha,x0,deltax)
                 ene1(i,j,k) = (rho2*cp2*vof + rho1*cp1*(-a1_l - vof))*T_l
                 eni1(i,j,k) = ene1(i,j,k)
              endif

              if(a1<0.d0) then
                 if (topo_mask(i-i0,j-j0,k-k0)/=0) then
                    deltax(d)=-a1
                    vof = fl3d(nr,alpha,x0,deltax)
                    eni1(i,j,k) = (rho2*cp2*vof + rho1*cp1*(-a1 - vof))*T_l
                 endif
              endif
              
              if (a2_l>0.d0) then
                 x0(d)=1.d0-a2_l
                 deltax(d)=a2_l
                 if (.not.(EtAdvScheme=='WENO5')) then
                    T_r = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                         energy(i+i0,j+j0,k+k0)/rhocp_r,AdvectionScheme,0.5d0*(1.d0-a2_l))
                 endif
                 vof = fl3d(nr,alpha,x0,deltax)
                 ene3(i,j,k) = (rho2*cp2*vof + rho1*cp1*(a2_l - vof))*T_r
                 eni3(i,j,k)=ene3(i,j,k)
              endif

              if(a2>0.d0) then
                 if (topo_mask(i+i0,j+j0,k+k0)/=0) then
                    x0(d)=1.d0-a2
                    deltax(d)=a2
                    vof = fl3d(nr,alpha,x0,deltax)
                    eni3(i,j,k) = (rho2*cp2*vof + rho1*cp1*(a2 - vof))*T_r
                 endif
              endif
           else              
              if(a1<0.d0) then
                 vof = -a1*c(i,j,k)
                 ene1(i,j,k) = (rho2*cp2*vof + rho1*cp1*(-a1 - vof))*T_l
                 eni1(i,j,k)=ene1(i,j,k)
              endif

              if (a1_l<0.d0) then
                 if (topo_mask(i-i0,j-j0,k-k0)==0) then
                    vof=-a1_l*c(i,j,k)
                    if (.not.(EtAdvScheme=='WENO5')) then  
                       T_l = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                            energy(i+i0,j+j0,k+k0)/rhocp_r,AdvectionScheme,-0.5d0*(1.d0+a1_l))
                    endif
                    eni1(i,j,k)=(rho2*cp2*vof + rho1*cp1*(-a1_l - vof))*T_l
                 endif
              endif

              if(a2>0.d0) then
                 vof = a2*c(i,j,k)
                 ene3(i,j,k) = (rho2*cp2*vof + rho1*cp1*(a2 - vof))*T_r
                 eni3(i,j,k) = ene3(i,j,k)
              endif

              if (a2_l>0.d0) then
                 if (topo_mask(i+i0,j+j0,k+k0)==0) then
                    vof=a2_l*c(i,j,k)
                    if (.not.(EtAdvScheme=='WENO5')) then
                       T_r = interpole3(energy(i-i0,j-j0,k-k0)/rhocp_l,energy(i,j,k)/rhocp_c,&
                            energy(i+i0,j+j0,k+k0)/rhocp_r,AdvectionScheme,0.5d0*(1.d0-a2_l))
                    endif
                    eni3(i,j,k)=(rho2*cp2*vof + rho1*cp1*(a2_l - vof))*T_r
                 endif
              endif

           endif
        enddo
     enddo
  enddo

  !Combine for new energy
  do k=ks,ke
     do j=js,je
        do i=is,ie
           a2 = us(i,j,k)*dt/dh
           a2_l = ul_adv(i,j,k)*dt/dh
           a1 = us(i-i0,j-j0,k-k0)*dt/dh
           a1_l = ul_adv(i-i0,j-j0,k-k0)*dt/dh

           T_avg = energy(i,j,k)/( rho2*cp2*c(i,j,k) + rho1*cp1*(1.d0-c(i,j,k)) )
           if (topo_mask(i,j,k)==0) then
              et_comp=(rho2*cp2*cg(i,j,k)+(1-cg(i,j,k))*rho1*cp1)*T_avg*(a2_l-a1_l)
           else
              et_comp=(rho2*cp2*cg(i,j,k)+(1-cg(i,j,k))*rho1*cp1)*T_avg*(a2-a1)
           endif

           de_adv = - (ene3(i,j,k) - eni1(i+i0,j+j0,k+k0)) + & 
                (eni3(i-i0,j-j0,k-k0) - ene1(i,j,k)) + et_comp
           energy(i,j,k) = energy(i,j,k) + de_adv

        enddo
     enddo
  enddo

end subroutine swpr_energy

subroutine swpz_energy(us,ul_adv,c,d)
  use module_grid
  use module_VOF
  use module_flow
  use module_BC
  use module_boil
  implicit none
  logical error
  integer i,j,k
  integer i0,j0,k0
  integer i1,j1,k1
  integer, intent(in) :: d
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(in) :: us, ul_adv, c
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax) :: en1,en2,en3,T_min,T_plus
  real(8) :: a1,a2,alpha,T_adv,a1_l,a2_l,vof_c,vol_c
  real(8) :: rhocp_ll, rhocp_l, rhocp_r, rhocp_rr, rhocp_c
  real(8) :: T_l, T_r
  REAL(8) :: deltax(3),x0(3),fl3d
  real(8) :: nr(3), stencil3x3(-1:1,-1:1,-1:1)

  logical :: weno_out

  intrinsic dmax1,dmin1
  !***
  call init_i0j0k0 (d,i0,j0,k0)
  if(ng.lt.2) call pariserror("wrong ng")

  weno_out=.false.

  do k=ks,ke
     do j=js,je
        do i=is,ie         
           rhocp_ll=(rho2*cp2*c(i-2*i0,j-2*j0,k-2*k0)+rho1*cp1*(1.d0-c(i-2*i0,j-2*j0,k-2*k0)))
           rhocp_l =(rho2*cp2*c(i-  i0,j-  j0,k-  k0)+rho1*cp1*(1.d0-c(i-  i0,j-  j0,k-  k0)))
           rhocp_c =(rho2*cp2*c(i     ,j     ,k     )+rho1*cp1*(1.d0-c(i     ,j     ,k     )))
           rhocp_r =(rho2*cp2*c(i+  i0,j+  j0,k+k0  )+rho1*cp1*(1.d0-c(i+i0  ,j+j0  ,k+k0  )))
           rhocp_rr=(rho2*cp2*c(i+2*i0,j+2*j0,k+2*k0)+rho1*cp1*(1.d0-c(i+2*i0,j+2*j0,k+2*k0)))

           weno_out=.false.

           T_adv = energy(i,j,k)/(rhocp_c)

           T_min(i,j,k)=interpoleWENO5(energy(i+i0*2,j+j0*2,k+k0*2)/(rhocp_rr),energy(i+i0,j+j0,k+k0)/(rhocp_r),&
                T_adv,energy(i-i0,j-j0,k-k0)/(rhocp_l),energy(i-i0*2,j-j0*2,k-k0*2)/(rhocp_ll),weno_out)
           T_plus(i,j,k)=interpoleWENO5(energy(i-i0*2,j-j0*2,k-k0*2)/(rhocp_ll),energy(i-i0,j-j0,k-k0)/(rhocp_l),&
                T_adv,energy(i+i0,j+j0,k+k0)/(rhocp_r),energy(i+i0*2,j+j0*2,k+k0*2)/(rhocp_rr),weno_out)

        enddo
     enddo
  enddo

  call SetBoundFaceTemps(d,T_min,T_plus)
  call do_all_ghost(T_min)
  call do_all_ghost(T_plus)

  do k=ks-1,ke+1
     do j=js-1,je+1
        do i=is-1,ie+1
           ! Note on u and u_l usage: The assumption is that u is divergence free everywhere, except for
           ! interface cells. Therefore, a1, a2 are used on all effluxes from mixed cells. a1_l and a2_l 
           ! are used on all compression terms for mixed cells, as well as effluxes from bulk cells into
           ! mixed cells (a1 and a2 are used for all other effluxes). 
           a2 = us(i,j,k)*dt/dxh(is)
           a2_l = ul_adv(i,j,k)*dt/dxh(is)
           a1 = us(i-i0,j-j0,k-k0)*dt/dxh(is)
           a1_l = ul_adv(i-i0,j-j0,k-k0)*dt/dxh(is)

           ! Get new heat capacity after vofsweep
           rhocp_c=(rho2*cp2*c(i   ,j   ,k   )+rho1*cp1*(1.d0-c(i   ,j   ,k   )))
           ! Interpolate for T-1/2 (T_l) and T+1/2 (T_r) using chosen AdvectionScheme
           ! For UpWind, T_l=T_r=T_adv

           if (us(i,j,k)>0.d0) then
              T_r=T_plus(i,j,k) 
           else
              T_r=T_min(i+i0,j+j0,k+k0)
           endif

           if (us(i-i0,j-i0,k-k0)>0.d0) then
              T_l=T_plus(i-i0,j-j0,k-k0) 
           else
              T_l=T_min(i,j,k)
           endif

           T_adv = energy(i,j,k)/(rhocp_c)
           
           en1(i,j,k)=0.d0; en2(i,j,k)=0.d0; en3(i,j,k)=0.d0
           if (topo_mask(i,j,k)==0) then
              do i1=-1,1; do j1=-1,1; do k1=-1,1
                 stencil3x3(i1,j1,k1) = c(i+i1,j+j1,k+k1)
              enddo;enddo;enddo
              call fit_plane_new(c(i,j,k),d,a1,a2,stencil3x3,nr,alpha,error)
              if(error) then
                 write(*,*)"WARNING: new plane error!"
                 cycle
              endif
              x0=0.d0
              deltax=1.d0
              if(a1<0.d0) then
                 x0(d)=a1
                 deltax(d)=-a1
                 vof_c = fl3d(nr,alpha,x0,deltax)
                 en1(i,j,k) = ( rho2*cp2*vof_c+rho1*cp1*(-a1-vof_c) )*T_l
              endif
              if(a2>0.d0) then
                 x0(d)=1.d0
                 deltax(d)=a2
                 vof_c = fl3d(nr,alpha,x0,deltax)
                 en3(i,j,k) = ( rho2*cp2*vof_c+rho1*cp1*(a2-vof_c) )*T_r             
              endif
              vol_c=(1.0d0-max(0.0d0,a1_l)-max(0.0d0,-a2_l));
              call fit_plane_new(c(i,j,k),d,a1_l,a2_l,stencil3x3,nr,alpha,error)
              x0(d) = max(0.0d0,a1_l)
              deltax(d) = vol_c
              vof_c = fl3d(nr,alpha,x0,deltax)
              en2(i,j,k) = ( rho2*cp2*vof_c+rho1*cp1*(vol_c-vof_c) )*T_adv
           else
              vol_c=(1.0d0-max(0.0d0,a1)-max(0.0d0,-a2));
              if (topo_mask(i-i0,j-j0,k-k0)==0) then
                 if(a1_l.lt.0d0) then
                    vof_c = c(i,j,k)*(-a1_l)
                    en1(i,j,k) = (rho2*cp2*vof_c + rho1*cp1*(-a1_l - vof_c))*T_l
                 endif
              else
                 if(a1.lt.0d0) then
                    vof_c = c(i,j,k)*(-a1)
                    en1(i,j,k) = (rho2*cp2*vof_c + rho1*cp1*(-a1 - vof_c))*T_l
                 endif
              endif
              if (topo_mask(i+i0,j+j0,k+k0)==0) then
                 if(a2_l.gt.0d0) then
                    vof_c = c(i,j,k)*a2_l
                    en3(i,j,k) = (rho2*cp2*vof_c + rho1*cp1*(a2_l - vof_c))*T_r
                 endif
              else
                 if(a2.gt.0d0) then
                    vof_c = c(i,j,k)*a2
                    en3(i,j,k) = (rho2*cp2*vof_c + rho1*cp1*(a2 - vof_c))*T_r
                 endif
              endif
              vof_c = c(i,j,k)*vol_c
              en2(i,j,k) = (rho2*cp2*vof_c + rho1*cp1*(vol_c - vof_c))*T_adv 
           endif 

        enddo
     enddo
  enddo

  do k=ks,ke
     do j=js,je
        do i=is,ie
           energy(i,j,k)  = en1(i+i0,j+j0,k+k0)+en2(i,j,k)+en3(i-i0,j-j0,k-k0)
        enddo
     enddo
  enddo
end subroutine swpz_energy

subroutine SetBoundFaceTemps(d,T_min,T_plus)
  use module_grid
  use module_boil
  implicit none
  real(8), DIMENSION(imin:imax,jmin:jmax,kmin:kmax), intent(inout) :: T_min,T_plus
  integer, intent(in) :: d
  integer :: i,j,k

  if (d==1) then
     if (coords(1)==0) then
        do j=js,je; do k=ks,ke
           if (BDRY_T(1)==0) then !Dirichlet
              T_min(is,j,k)=BC_T(1)
              T_plus(is-1,j,k)=BC_T(1)
           endif
           if (BDRY_T(1)==1) then !Neumann
              T_min(is,j,k)=Te(is,j,k)-BC_T(1)*dxh(is)*0.5d0
              T_plus(is-1,j,k)=T_min(is,j,k)
           endif
        enddo; enddo
     endif

     if (coords(1)==nPx-1) then
        do j=js,je; do k=ks,ke
           if (BDRY_T(4)==0) then 
              T_plus(ie,j,k)=BC_T(4)
              T_min(ie+1,j,k)=BC_T(4)
           endif
           if (BDRY_T(4)==1) then 
              T_plus(ie,j,k)=Te(ie,j,k)+BC_T(4)*dxh(ie)*0.5d0
              T_min(ie+1,j,k)=T_plus(ie,j,k)
           endif
        enddo; enddo
     endif
  endif

  if (d==2) then
     if (coords(2)==0) then
        do i=is,ie; do k=ks,ke
           if (BDRY_T(2)==0) then !Dirichlet
              T_min(i,js,k)=BC_T(2)
              T_plus(i,js-1,k)=BC_T(2)
           endif
           if (BDRY_T(2)==1) then
              T_min(i,js,k)=Te(i,js,k)-BC_T(2)*dyh(js-1)*0.5d0
              T_plus(i,js-1,k)=T_min(i,js,k)
           endif
        enddo; enddo
     endif

     !BC on y+
     if (coords(2)==nPy-1) then
        do i=is,ie; do k=ks,ke
           if (BDRY_T(5)==0) then 
              T_plus(i,je,k)=BC_T(5)
              T_min(i,je+1,k)=BC_T(5)
           endif
           if (BDRY_T(5)==1)then
              T_plus(i,je,k)=Te(i,je,k)+BC_T(5)*dyh(je)*0.5d0
              T_min(i,je+1,k)=T_plus(i,je,k)
           endif
        enddo; enddo
     endif
  endif

  if (d==3) then
     if (coords(3)==0) then
        do i=is,ie; do j=js,je
           if (BDRY_T(3)==0) then !Dirichlet
              T_min(i,j,ks)=BC_T(3)
              T_plus(i,j,ks-1)=BC_T(3)
           endif
           if (BDRY_T(3)==1) then
              T_min(i,j,ks)=Te(i,j,ks)-BC_T(3)*dzh(ks-1)*0.5d0
              T_plus(i,j,ks-1)=T_min(i,j,ks)
           endif
        enddo; enddo
     endif

     !BC on z+
     if (coords(3)==nPz-1) then
        do i=is,ie; do j=js,je
           if (BDRY_T(6)==0) then 
              T_plus(i,j,ke)=BC_T(6)
              T_min(i,j,ke+1)=BC_T(6)
           endif
           if (BDRY_T(6)==1) then 
              T_plus(i,j,ke)=Te(i,j,ke)+BC_T(6)*dzh(ke)*0.5d0
              T_min(i,j,ke+1)=T_plus(i,j,ke)
           endif
        enddo; enddo
     endif
  endif

end subroutine SetBoundFaceTemps

subroutine dE_PhaseChange(rho1,rho2)
  use module_grid
  use module_boil
  use module_BC
  implicit none
  real(8), intent(in) :: rho1,rho2
  real(8) :: dE_Evap
  integer :: i,j,k

  do i=is,ie; do j=js,je;  do k=ks,ke
     dE_Evap=T_sat*(rho2*Cp2-rho1*Cp1)*dc(i,j,k)
     energy(i,j,k) = energy(i,j,k)+dE_Evap
  enddo; enddo; enddo

  call do_all_ghost(energy)
end subroutine dE_PhaseChange

subroutine get_Te_from_energy
  use module_grid
  use module_boil
  implicit none
  integer :: i,j,k

  do k=ks,ke
     do j=js,je
        do i=is,ie
           Te(i,j,k)=energy(i,j,k)/rho_cp(i,j,k)
        enddo
     enddo
  enddo

end subroutine get_Te_from_energy

subroutine get_energy
  use module_grid
  use module_boil
  implicit none
  integer :: i,j,k

  do i=is,ie; do j=js,je;  do k=ks,ke
     energy(i,j,k) = rho_cp(i,j,k)*Te(i,j,k)
  enddo; enddo; enddo

end subroutine get_energy

!Fill boundary energies using Te and CVOF in the boundaries
subroutine SetBoundaryEnergy
  use module_grid
  use module_boil
  use module_flow
  use module_VOF
  implicit none
integer :: i,j,k

  if(coords(1)==0) then
     do j=jmin,jmax
        do k=kmin,kmax           
           energy(is-1,j,k)=( rho2*cp2*cvof(is-1,j,k) + rho1*cp1*(1.d0-cvof(is-1,j,k)) ) * Te(is-1,j,k)
           energy(is-2,j,k)=( rho2*cp2*cvof(is-2,j,k) + rho1*cp1*(1.d0-cvof(is-2,j,k)) ) * Te(is-2,j,k)
        enddo
     enddo
  endif


  if(coords(2)==0) then
     do i=imin,imax
        do k=kmin,kmax
           energy(i,js-1,k)=( rho2*cp2*cvof(i,js-1,k) + rho1*cp1*(1.d0-cvof(i,js-1,k)) ) * Te(i,js-1,k)
           energy(i,js-2,k)=( rho2*cp2*cvof(i,js-2,k) + rho1*cp1*(1.d0-cvof(i,js-2,k)) ) * Te(i,js-2,k)
        enddo
     enddo
  endif

  if(coords(3)==0) then
     do i=imin,imax
        do j=jmin,jmax
           energy(i,j,ks-1)=( rho2*cp2*cvof(i,j,ks-1) + rho1*cp1*(1.d0-cvof(i,j,ks-1)) ) * Te(i,j,ks-1)
           energy(i,j,ks-2)=( rho2*cp2*cvof(i,j,ks-2) + rho1*cp1*(1.d0-cvof(i,j,ks-2)) ) * Te(i,j,ks-2) 
        enddo
     enddo
  endif

  if(coords(1)==nPx-1) then
     do j=jmin,jmax
        do k=kmin,kmax
           energy(ie+1,j,k)=( rho2*cp2*cvof(ie+1,j,k) + rho1*cp1*(1.d0-cvof(ie+1,j,k)) ) * Te(ie+1,j,k)
           energy(ie+2,j,k)=( rho2*cp2*cvof(ie+2,j,k) + rho1*cp1*(1.d0-cvof(ie+2,j,k)) ) * Te(ie+2,j,k) 
        enddo
     enddo
  endif

  if(coords(2)==nPy-1) then
     do i=imin,imax
        do k=kmin,kmax
           energy(i,je+1,k)=( rho2*cp2*cvof(i,je+1,k) + rho1*cp1*(1.d0-cvof(i,je+1,k)) ) * Te(i,je+1,k)
           energy(i,je+2,k)=( rho2*cp2*cvof(i,je+2,k) + rho1*cp1*(1.d0-cvof(i,je+2,k)) ) * Te(i,je+2,k) 
        enddo
     enddo
  endif

  if(coords(3)==nPz-1) then
     do i=imin,imax
        do j=jmin,jmax
           energy(i,j,ke+1)=( rho2*cp2*cvof(i,j,ke+1) + rho1*cp1*(1.d0-cvof(i,j,ke+1)) ) * Te(i,j,ke+1)
           energy(i,j,ke+2)=( rho2*cp2*cvof(i,j,ke+2) + rho1*cp1*(1.d0-cvof(i,j,ke+2)) ) * Te(i,j,ke+2) 
        enddo
     enddo
  endif
end subroutine SetBoundaryEnergy

!Boundary values for temperature
subroutine SetTempBc
  use module_grid
  use module_boil
  implicit none
  integer :: i,j,k

  Tumask=1.0d0; Tvmask=1.0d0; Twmask=1.0d0
  
  !BC on x-
  if (coords(1)==0) then
     if (BDRY_T(1)<2) Tumask(is-1,jmin:jmax,kmin:kmax)=0.0d0
     do j=jmin,jmax; do k=kmin,kmax
        if (BDRY_T(1)==0) then !Dirichlet
           Te(is-1,j,k)=2.0d0*BC_T(1)-Te(is,j,k)
           Te(is-2,j,k)=Te(is-1,j,k)
        endif
        if (BDRY_T(1)==1) then
           Te(is-1,j,k)=Te(is,j,k)-BC_T(1)*dxh(is-1)
           Te(is-2,j,k)=Te(is-1,j,k)-BC_T(1)*dxh(is-1)
        endif
     enddo; enddo
  endif

  !BC on x+
  if (coords(1)==nPx-1) then
     if (BDRY_T(4)<2) Tumask(ie,jmin:jmax,kmin:kmax)=0.0d0
     do j=jmin,jmax; do k=kmin,kmax
        if (BDRY_T(4)==0) then 
           Te(ie+1,j,k)=2.0d0*BC_T(4)-Te(ie,j,k)
           Te(ie+2,j,k)=Te(ie+1,j,k)
        endif
        if (BDRY_T(4)==1) then
           Te(ie+1,j,k)=Te(ie,j,k)+BC_T(4)*dxh(ie)
           Te(ie+2,j,k)=Te(ie+1,j,k)+BC_T(4)*dxh(ie)
        endif
     enddo; enddo
  endif

  !BC on y-
  if (coords(2)==0) then
     if (BDRY_T(2)<2) Tvmask(imin:imax,js-1,kmin:kmax)=0.0d0
     do i=imin,imax; do k=kmin,kmax
        if (BDRY_T(2)==0) then !Dirichlet
           Te(i,js-1,k)=2.0d0*BC_T(2)-Te(i,js,k)
           Te(i,js-2,k)=Te(i,js-1,k)
        endif
        if (BDRY_T(2)==1) then
           Te(i,js-1,k)=Te(i,js,k)-BC_T(2)*dyh(js-1)
           Te(i,js-2,k)=Te(i,js-1,k)-BC_T(2)*dyh(js-1)
        endif
     enddo; enddo
  endif

  !BC on y+
  if (coords(2)==nPy-1) then
     if (BDRY_T(5)<2) Tvmask(imin:imax,je,kmin:kmax)=0.0d0
     do i=imin,imax; do k=kmin,kmax
        if (BDRY_T(5)==0) then 
           Te(i,je+1,k)=2.0d0*BC_T(5)-Te(i,je,k)
           Te(i,je+2,k)=Te(i,je+1,k)
        endif
        if (BDRY_T(5)==1) then 
           Te(i,je+1,k)=Te(i,je,k)+BC_T(5)*dyh(je)
           Te(i,je+2,k)=Te(i,je+1,k)+BC_T(5)*dyh(je)
        endif
     enddo; enddo
  endif

  !BC on z-
  if (coords(3)==0) then
     if (BDRY_T(3)<2) Twmask(imin:imax,jmin:jmax,ks-1)=0.0d0
     do i=imin,imax; do j=jmin,jmax
        if (BDRY_T(3)==0) then !Dirichlet
           Te(i,j,ks-1)=2.0d0*BC_T(3)-Te(i,j,ks)
           Te(i,j,ks-2)=Te(i,j,ks-1)
        endif
        if (BDRY_T(3)==1) then
           Te(i,j,ks-1)=Te(i,j,ks)-BC_T(3)*dzh(ks-1)
           Te(i,j,ks-2)=Te(i,j,ks-1)-BC_T(3)*dzh(ks-1)
        endif
     enddo; enddo
  endif

  !BC on z+
  if (coords(3)==nPz-1) then
     if (BDRY_T(6)<2) Twmask(imin:imax,jmin:jmax,ke)=0.0d0
     do i=imin,imax; do j=jmin,jmax
        if (BDRY_T(6)==0) then 
           Te(i,j,ke+1)=2.0d0*BC_T(6)-Te(i,j,ke)
           Te(i,j,ke+2)=Te(i,j,ke+1)
        endif
        if (BDRY_T(6)==1) then 
           Te(i,j,ke+1)=Te(i,j,ke)+BC_T(6)*dzh(ke)
           Te(i,j,ke+2)=Te(i,j,ke+1)+BC_T(6)*dzh(ke)
        endif
     enddo; enddo
  endif

end subroutine SetTempBC

subroutine setuppoisson_phase_change(coeff)
  use module_grid
  use module_flow
  use module_boil
  implicit none
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: coeff
  integer :: i,j,k,i0,j0,k0,i1,j1,k1,d,sign
  integer :: ierr

  do k=ks,ke; do j=js,je; do i=is,ie
     coeff(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rho(i-1,j,k)+rho(i,j,k)))
     coeff(i,j,k,2) = 2d0*dt*umask(i  ,j,k)/(dx(i)*dxh(i  )*(rho(i+1,j,k)+rho(i,j,k)))
     coeff(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rho(i,j-1,k)+rho(i,j,k)))
     coeff(i,j,k,4) = 2d0*dt*vmask(i,j  ,k)/(dy(j)*dyh(j  )*(rho(i,j+1,k)+rho(i,j,k)))
     coeff(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rho(i,j,k-1)+rho(i,j,k)))
     coeff(i,j,k,6) = 2d0*dt*wmask(i,j,k  )/(dz(k)*dzh(k  )*(rho(i,j,k+1)+rho(i,j,k)))
     coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
     coeff(i,j,k,8) = s_v(i,j,k) - ( (u(i,j,k)-u(i-1,j,k))/dx(i) &
          +  (v(i,j,k)-v(i,j-1,k))/dy(j) &
          +  (w(i,j,k)-w(i,j,k-1))/dz(k) )
     !===============================================================================================================
  enddo; enddo; enddo
  call Poisson_BCs(coeff)
end subroutine setuppoisson_phase_change

subroutine setuppoisson_phase_change_ext(coeff)
  use module_grid
  use module_flow
  use module_boil
  implicit none
  real(8), dimension(is:ie,js:je,ks:ke,8), intent(out) :: coeff
  integer :: i,j,k,i0,j0,k0,i1,j1,k1,d,sign
  integer :: ierr

  coeff=0.0d0
  do k=ks,ke; do j=js,je; do i=is,ie
     if ( topo_mask(i,j,k)==-2 .or. topo_mask(i,j,k)==-1 ) then
        coeff(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rho(i-1,j,k)+rho(i,j,k)))
        coeff(i,j,k,2) = 2d0*dt*umask(i  ,j,k)/(dx(i)*dxh(i  )*(rho(i+1,j,k)+rho(i,j,k)))
        coeff(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rho(i,j-1,k)+rho(i,j,k)))
        coeff(i,j,k,4) = 2d0*dt*vmask(i,j  ,k)/(dy(j)*dyh(j  )*(rho(i,j+1,k)+rho(i,j,k)))
        coeff(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rho(i,j,k-1)+rho(i,j,k)))
        coeff(i,j,k,6) = 2d0*dt*wmask(i,j,k  )/(dz(k)*dzh(k  )*(rho(i,j,k+1)+rho(i,j,k)))
        coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
        coeff(i,j,k,8) = 0.0d0
     else if (topo_mask(i,j,k)==0) then
        if (topo_mask(i-1,j,k)==-1 .or. topo_mask(i-1,j,k)==0) then
           coeff(i,j,k,1) = 2d0*dt*umask(i-1,j,k)/(dx(i)*dxh(i-1)*(rho(i-1,j,k)+rho(i,j,k)))
        endif
        if (topo_mask(i+1,j,k)==-1 .or. topo_mask(i+1,j,k)==0) then
           coeff(i,j,k,2) = 2d0*dt*umask(i  ,j,k)/(dx(i)*dxh(i  )*(rho(i+1,j,k)+rho(i,j,k)))
        endif
        if (topo_mask(i,j-1,k)==-1 .or. topo_mask(i,j-1,k)==0) then
           coeff(i,j,k,3) = 2d0*dt*vmask(i,j-1,k)/(dy(j)*dyh(j-1)*(rho(i,j-1,k)+rho(i,j,k)))
        endif
        if (topo_mask(i,j+1,k)==-1 .or. topo_mask(i,j+1,k)==0) then
           coeff(i,j,k,4) = 2d0*dt*vmask(i,j  ,k)/(dy(j)*dyh(j  )*(rho(i,j+1,k)+rho(i,j,k)))
        endif
        if (topo_mask(i,j,k-1)==-1 .or. topo_mask(i,j,k-1)==0) then
           coeff(i,j,k,5) = 2d0*dt*wmask(i,j,k-1)/(dz(k)*dzh(k-1)*(rho(i,j,k-1)+rho(i,j,k)))
        endif
        if (topo_mask(i,j,k+1)==-1 .or. topo_mask(i,j,k+1)==0) then
           coeff(i,j,k,6) = 2d0*dt*wmask(i,j,k  )/(dz(k)*dzh(k  )*(rho(i,j,k+1)+rho(i,j,k)))
        endif
        coeff(i,j,k,7) = sum(coeff(i,j,k,1:6))
        coeff(i,j,k,8) = -s_v(i,j,k)
     else
        coeff(i,j,k,1) = 0.0d0
        coeff(i,j,k,2) = 0.0d0
        coeff(i,j,k,3) = 0.0d0
        coeff(i,j,k,4) = 0.0d0
        coeff(i,j,k,5) = 0.0d0
        coeff(i,j,k,6) = 0.0d0
        coeff(i,j,k,7) = 1.0d0
        coeff(i,j,k,8) = 0.0d0
     endif

  enddo; enddo; enddo
end subroutine setuppoisson_phase_change_ext

subroutine phase_change_vel_correction
  use module_grid
  use module_flow
  use module_boil
  implicit none
  integer :: i,j,k

  do k=ks,ke;  do j=js,je; do i=is,ieu    ! CORRECT THE u-velocity 
     if ( topo_mask(i,j,k)>-3 .and. topo_mask(i,j,k)<=0 ) then
        if ( correct(1,i,j,k) ) then
           u_s(i,j,k)=-dt*(2.0*umask(i,j,k)/dxh(i))*(p_s(i+1,j,k)-p_s(i,j,k))/(rho(i+1,j,k)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo

  do k=ks,ke;  do j=js,jev; do i=is,ie    ! CORRECT THE v-velocity
     if ( topo_mask(i,j,k)>-3 .and. topo_mask(i,j,k)<=0 ) then
        if ( correct(2,i,j,k) ) then
           v_s(i,j,k)=-dt*(2.0*vmask(i,j,k)/dyh(j))*(p_s(i,j+1,k)-p_s(i,j,k))/(rho(i,j+1,k)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo

  do k=ks,kew;  do j=js,je; do i=is,ie   ! CORRECT THE w-velocity
     if ( topo_mask(i,j,k)>-3 .and. topo_mask(i,j,k)<=0 ) then
        if ( correct(3,i,j,k) ) then
           w_s(i,j,k)=-dt*(2.0*wmask(i,j,k)/dzh(k))*(p_s(i,j,k+1)-p_s(i,j,k))/(rho(i,j,k+1)+rho(i,j,k))
        endif
     endif
  enddo; enddo; enddo

contains
  function correct(d,i,j,k)
    integer, intent(in) :: d,i,j,k
    integer :: i0,j0,k0
    logical :: correct
    call init_i0j0k0 (d,i0,j0,k0)
    correct = ( (topo_mask(i,j,k)==-2 .and. topo_mask(i+i0,j+j0,k+k0)==-1) &
         .or. (topo_mask(i,j,k)==-1) &
         .or. (topo_mask(i,j,k)==0 .and. &
         (topo_mask(i+i0,j+j0,k+k0)==0 .or. topo_mask(i+i0,j+j0,k+k0)==-1) ) )
  end function correct
end subroutine phase_change_vel_correction

subroutine get_ghost_masks
  use module_grid
  use module_VOF
  use module_BC
  use module_tmpvar
  use module_boil
  implicit none
  include 'mpif.h'
  integer :: phase
  integer :: i,j,k,d
  integer :: i0,j0,k0,i1,j1,k1,level
  integer :: req(12),sta(MPI_STATUS_SIZE,12),ierr,phase_sign

  !First loop to set level 0 velocities using staggered phase
  !initialize all masks to 3
  topo_mask = 3*(2*vof_phase - 1)
  
  do k=ks,ke; do j=js,je; do i=is,ie
     if (vof_flag(i,j,k) == 2) topo_mask(i,j,k) = 0
  enddo; enddo; enddo
  call do_all_ighost(topo_mask)

  !Set levels 1 to 2
  do level=1,2
     do k=ks,ke; do j=js,je; do i=is,ie
        do phase=0,1
           phase_sign=2*phase-1
           if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
              if ((topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                   (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                   (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                   (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                   (topo_mask(i,j,k+1)==phase_sign*(level-1)).or.&
                   (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                 topo_mask(i,j,k)=level*phase_sign
              endif
           endif
        enddo
     enddo; enddo; enddo
     call do_all_ighost(topo_mask)
  enddo

  ! Get boundary cell masks, assuming dc/dn=0
  if(coords(1)==0) then
     do level=1,2
        do k=ks,ke; do j=js,je; do i=imin,is-1 ! ng must be 2
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ((topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k+1)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif
  if(coords(1)==nPx-1) then
     do level=1,2
        do k=ks,ke; do j=js,je; do i=ie+1,imax
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ( (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k+1)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif
  if(coords(2)==0) then
     do level=1,2 
        do k=ks,ke; do i=is,ie; do j=jmin,js-1
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ( (topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k+1)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif
  if(coords(2)==nPy-1) then
     do level=1,2
        do k=ks,ke; do i=is,ie; do j=je+1,jmax
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ((topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k+1)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif
  if(coords(3)==0) then
     do level=1,2
        do j=js,je; do i=is,ie; do k=kmin,ks-1
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ((topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k+1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif
  if(coords(3)==nPz-1) then
     do level=1,2
        do j=js,je; do i=is,ie; do k=ke+1,kmax
           do phase=0,1
              phase_sign=2*phase-1
              if (topo_mask(i,j,k)==3*phase_sign .and. vof_phase(i,j,k)==phase) then
                 if ((topo_mask(i+1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i-1,j,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j+1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j-1,k)==phase_sign*(level-1)).or.&
                      (topo_mask(i,j,k-1)==phase_sign*(level-1))) then
                    topo_mask(i,j,k)=level*phase_sign
                 endif
              endif
           enddo
        enddo; enddo; enddo
     enddo
  endif

end subroutine get_ghost_masks

subroutine sweep_mass_transfer(rho_l,dt,dh)
  use module_grid
  use module_boil
  use module_VOF
  use module_BC
  implicit none
  real(8), intent(in) :: rho_l,dt,dh
  real(8), dimension(imin:imax,jmin:jmax,kmin:kmax) :: nshift
  integer :: i,j,k
  nshift=mdot/rho_l*dt/dh
  call LevelSetShift2VOFCHANGE(nshift,dc)

end subroutine sweep_mass_transfer

subroutine NuAvg(t)
  use module_grid
  use module_boil
  implicit none
  include 'mpif.h'
  real(8), intent(in) :: t
  real(8) :: count, Nu_loc, NuTot, total
  integer :: i,j,k,ierr

  count=0.0d0; Nu_loc=0.0d0; 
  if (coords(2)==0) then
     do k=ks,ke; do i=is,ie
        count=count+1.0d0
        Nu_loc=Nu_loc + 2.d0*(BC_T(2)-Te(i,js,k))/dy(js)
     enddo; enddo
  endif
  call MPI_ALLREDUCE(Nu_loc, NuTot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  call MPI_ALLREDUCE(count, total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_Domain, ierr)
  NuTot=NuTot/total
  if (rank==0) then
     open(unit=70,file='AvgNu',position='append')
     write(70,'(2e14.5)')t,NuTot
     close(70)
  endif
end subroutine NuAvg

#endif

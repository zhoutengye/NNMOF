Module ModMOF
  Use ModGlobal
  Use SussmanMOF
  Implicit None
  Integer, Parameter  :: MofItermax = 1
  Real(sp), Parameter :: tol = 1.0e-8
  Real(sp), Parameter :: local_tol = 1.0e-8
  Real(sp), Parameter :: CenTol = 1e-13
  Real(sp), Parameter :: GaussNewtonTol = 1e-13
  Real(sp), Parameter :: MOF_PI = 3.1415926535897932d0
Contains
  Subroutine MOF_Init
    Implicit None
    Call Initialize_SussmanMOF
  End Subroutine MOF_Init

  Subroutine NormSussmanMOF(f, c, norm, init_norm)
    Use geometry_intersect_module
    Use MOF_routines_module
    Implicit None
    Real(sp) , Intent(In)   :: f
    Real(sp) , Intent(In)   :: c(3)
    Real(sp) , Intent(out)  :: norm(3)
    Real(sp) , Intent(In), optional  :: init_norm(3)

    mofdata = 0.0_sp
    mofdata(1) = f
    mofdata(2) = c(1) - 0.5_sp
    mofdata(3) = c(2) - 0.5_sp
    mofdata(4) = c(3) - 0.5_sp
    mofdata(5) = 1.0_sp

    mofdata(10)  = 1.0-f
    mofdata(11) = - f * c(1) / mofdata(10)
    mofdata(12) = - f * c(2) / mofdata(10)
    mofdata(13) = - f * c(3) / mofdata(10)
    mofdata(14) = 2.0_sp

    Call multimaterial_MOF( &
        bfact,dl,xsten0,nhalf0, &
        mof_verbose, &
        use_ls_data, &
        LS_stencil, &
        xtetlist_vof, &
        xtetlist_cen, &
        nmax, &
        mofdata, &
        multi_centroidA, &
        continuous_mof, &
        levelrz,nmat,sdim, &
        ngeom_recon, &
        caller_id)

    norm(1) = -mofdata(6)
    norm(2) = -mofdata(7)
    norm(3) = -mofdata(8)

  End Subroutine NormSussmanMOF

  Subroutine MOFZY(f, c, norm)
    Implicit None
    Real(sp), Intent(In)  :: f
    Real(sp), Intent(In)  :: c(3)
    Real(sp), Intent(Out) :: norm(3)

    Real(sp), Dimension(3)   :: norm_2
    Real(sp) :: delta_theta
    Real(sp) :: delta_theta_max
    Real(sp), Dimension(2)   :: delangle, angle_init, angle_previous, new_angle
    Real(sp), Dimension(2)   :: angle_base, angle_plus, angle_minus
    Real(sp), Dimension(2)   :: err_plus, err_minus
    Real(sp), Dimension(3)   :: dbase, dopt, dp, dm
    Real(sp), Dimension(3,2) :: d_plus, d_minus
    Real(sp), Dimension(3)   :: cenopt, cenp, cenm, cen_init
    Real(sp), Dimension(3,2) :: cen_plus, cen_minus
    Real(sp), Dimension(3,2) :: dgrad
    Real(sp), Dimension(2,2) :: JTJ, JTJINV
    Real(sp), Dimension(2)   :: RHS
    Real(sp) :: DET
    Real(sp), Dimension(3,MOFITERMAX+1) :: d_array, cen_array
    Real(sp), Dimension(2,MOFITERMAX+1) :: angle_array
    Real(sp), Dimension(MOFITERMAX+1)   :: err_array

    Real(sp) :: err, err_local_min
    Integer :: singular_flag
    Integer :: i_angle, j_angle, dir, iter, i

    ! Setting the initial guess as 0
    angle_init = 0.0_sp
    do dir=1,2
      angle_array(dir,1)= angle_init(dir)
    enddo
    Call FindCentroid(angle_init,f,cen_init)
    do dir=1,3
      d_array(dir,1) = c(dir) - cen_init(dir)
      cen_array(dir,1)=cen_init(dir)
      err_array(1)=err_array(1) + d_array(dir,1) * d_array(dir,1)
    enddo

    delta_theta=MOF_Pi/180.0  ! 1 degree=pi/180
    delta_theta_max=10.0*MOF_Pi/180  ! 10 degrees

    iter = 0
    err = err_array(1)
    err_local_min = err_array(1)
    Do While ((iter.lt.MOFITERMAX).and. &
        (err.gt.tol).and. &
        (err_local_min.gt.local_tol))

      do dir=1,3
        dbase(dir)=d_array(dir,iter+1)
      enddo

      do i_angle=1,2
        angle_base(i_angle)=angle_array(i_angle,iter+1)
      enddo

      Do i_angle=1,2

        Do j_angle=1,2
          angle_plus(j_angle)=angle_base(j_angle)
          angle_minus(j_angle)=angle_base(j_angle)
        End Do
        angle_plus(i_angle)=angle_plus(i_angle)+delta_theta
        angle_minus(i_angle)=angle_minus(i_angle)-delta_theta
        Call FindCentroid(angle_plus,f,cenp)
        err_plus(i_angle)= 0.0_sp
        do dir=1,3
          dp(dir) = c(dir) - cenp(dir)
          err_plus(i_angle)=err_plus(i_angle)+dp(dir)**2
          d_plus(dir,i_angle)=dp(dir)
          cen_plus(dir,i_angle)=cenp(dir)
        end do
        err_plus(i_angle)=sqrt(err_plus(i_angle))
        Call FindCentroid(angle_minus,f,cenm)
        err_minus(i_angle)= 0.0_sp
        do dir=1,3
          dm(dir) = c(dir) - cenm(dir)
          err_minus(i_angle)=err_minus(i_angle)+dm(dir)**2
          d_minus(dir,i_angle)=dm(dir)
          cen_minus(dir,i_angle)=cenm(dir)
        enddo
        err_minus(i_angle)=sqrt(err_minus(i_angle))
        !!! jacobian matrix has:
        do dir=1,3
          dgrad(dir,i_angle)=(dp(dir)-dm(dir))/(2.0_sp*delta_theta)
        enddo

      End Do

      Do i_angle=1,2
        Do j_angle=1,2
          JTJ(i_angle,j_angle)= 0.0_sp
          Do dir=1,3
            JTJ(i_angle,j_angle)=JTJ(i_angle,j_angle)+ &
                dgrad(dir,i_angle)*dgrad(dir,j_angle)
          EndDo
        EndDo
      EndDo

      DET=JTJ(1,1)*JTJ(2,2)-JTJ(1,2)*JTJ(2,1)
      ! DET has dimensions of length squared
      if (abs(DET).ge.CENTOL) then 
        singular_flag=0
      else if (abs(DET).le.CENTOL) then
        singular_flag=1
      else
        print *,"DET bust"
        stop
      endif

      if (singular_flag.eq.0) then
        JTJINV(1,1)=JTJ(2,2)
        JTJINV(2,2)=JTJ(1,1)
        JTJINV(1,2)=-JTJ(1,2)
        JTJINV(2,1)=-JTJ(2,1)

        do i_angle=1,2
          do j_angle=1,2
            JTJINV(i_angle,j_angle)=JTJINV(i_angle,j_angle)/DET
          enddo
        enddo

        do i_angle=1,2  ! compute -JT * r
          RHS(i_angle)= 0.0_sp
          do dir=1,3
            RHS(i_angle)=RHS(i_angle)-dgrad(dir,i_angle)*dbase(dir)
          enddo
        enddo

        err_local_min= 0.0_sp
        do i_angle=1,2
          err_local_min=err_local_min+RHS(i_angle)**2
        enddo
        err_local_min=sqrt(err_local_min)

        do i_angle=1,2  ! compute JTJ^-1 (RHS)
          delangle(i_angle)= 0.0_sp
          do j_angle=1,2
            delangle(i_angle)=delangle(i_angle)+ &
                JTJINV(i_angle,j_angle)*RHS(j_angle)
          enddo
          if (delangle(i_angle).gt.delta_theta_max) then
            delangle(i_angle)=delta_theta_max
          else if (delangle(i_angle).lt.-delta_theta_max) then
            delangle(i_angle)=-delta_theta_max
          endif
        enddo
        ! -pi<angle<pi
        do i_angle=1,2
          angle_previous(i_angle)=angle_base(i_angle)
          call advance_angle(angle_base(i_angle),delangle(i_angle))
        enddo
        Call FindCentroid(angle_base,f, cenopt)
        do dir=1,3
          dopt(dir) = c(dir) - cenopt(dir)
        enddo

      else if (singular_flag.eq.1) then
        err_local_min= 0.0_sp

        do dir=1,3
          dopt(dir)=dbase(dir)
        enddo
        do i_angle=1,2
          angle_base(i_angle)=angle_array(i_angle,iter+1)
        enddo
      else
        print *,"singular_flag invalid"
        stop
      End If

      err= 0.0_sp
      do dir=1,3
        err=err+dopt(dir)**2
      enddo
      err=sqrt(err)

      if (singular_flag.eq.1) then
        ! do nothing
      else if (singular_flag.eq.0) then
         do i_angle=1,2
          do j_angle=1,2
           angle_plus(j_angle)=angle_previous(j_angle)
           angle_minus(j_angle)=angle_previous(j_angle)
          enddo
          angle_plus(i_angle)=angle_plus(i_angle)+delta_theta
          angle_minus(i_angle)=angle_minus(i_angle)-delta_theta

          if ((err.le.err_plus(i_angle)).and. &
              (err.le.err_minus(i_angle))) then
           ! do nothing
          else

           if (err.ge.err_plus(i_angle)) then
            err=err_plus(i_angle)
            do dir=1,3
             dopt(dir)=d_plus(dir,i_angle)
             cenopt(dir)=cen_plus(dir,i_angle)
            enddo 
            do j_angle=1,2
             angle_base(j_angle)=angle_plus(j_angle)
            enddo
           endif

           if (err.ge.err_minus(i_angle)) then
            err=err_minus(i_angle)
            do dir=1,3
             dopt(dir)=d_minus(dir,i_angle)
             cenopt(dir)=cen_minus(dir,i_angle)
            enddo
            do j_angle=1,2
             angle_base(j_angle)=angle_minus(j_angle)
            enddo
           endif

          endif ! use safe update

         enddo ! i_angle     if (singular_flag.eq.1) then
       else
         print *,"singular_flag invalid"
         stop
       endif

       do dir=1,3
         d_array(dir,iter+2)=dopt(dir)
         cen_array(dir,iter+2)=cenopt(dir)
       enddo

       do i_angle=1,2
         angle_array(i_angle,iter+2)=angle_base(i_angle)
       enddo
       err_array(iter+2)=err
       iter=iter+1
       print *, 'iter=',iter
    End Do

    do i = 1, mofitermax+1
      print *, '=====step',i,'=========='
      print *, cen_array(:,i)
      print *, angle_array(:,i)
      print *, err_array(i)
    end do

    do dir=1,2
      new_angle(dir)=angle_array(dir,iter+1)
    enddo

    call Angle2Norm(new_angle,norm_2)

    scale = 1.0_sp / abs(norm_2(1)+abs(norm_2(2)) + abs(norm_2(3)))
    norm = norm_2 * scale

  End Subroutine MOFZY

  Subroutine Angle2Norm(angle, norm)
    Implicit None
    Real(sp), Intent(In)  :: angle(2)
    Real(sp), Intent(Out) :: norm(3)
    norm(3)=cos(angle(2))
    norm(1)=sin(angle(2))*cos(angle(1))
    norm(2)=sin(angle(2))*sin(angle(1))
  End Subroutine Angle2Norm

  Subroutine MOF2SZ(norm1, alpha, norm2, dis)
    Implicit None
    Real(sp), Intent(In)  :: alpha
    Real(sp), Intent(In)  :: norm1(3)
    Real(sp), Intent(Out) :: norm2(3)
    Real(sp), Intent(Out) :: dis
    Real(sp) :: scale

    scale = 1.0_sp / sqrt( norm1(1)*norm1(1) + norm1(2)*norm1(2) + norm1(3)*norm1(3))
    norm2 = norm1 * scale

  End Subroutine MOF2SZ

  Subroutine FindCentroid(angle, vof, cen_mof)
    Implicit None
    Real(sp), Intent(In) :: angle(2)
    Real(sp), Intent(In) :: vof
    Real(sp), Intent(Out) :: cen_mof(3)
    Real(sp) :: cen_sz(3)
    Real(sp) :: norm_2(3)
    Real(sp) :: norm_1(3)
    Real(sp) :: alpha
    Real(sp) :: x0(3)
    Real(sp) :: deltax(3)
    Real(sp) :: FloodSZ_Backward

    Call Angle2Norm(angle, norm_2)
    scale = 1.0_sp / (abs(norm_2(1))+abs(norm_2(2)) + abs(norm_2(3)))
    norm_1 = norm_2 * scale
    alpha = FloodSZ_Backward(norm_1,vof)
    x0 = 0.0_sp
    deltax = 1.0_sp
    Call FloodSZ_ForwardC(norm_1,alpha,x0,deltax,vof,cen_sz)
    cen_mof = cen_sz - 0.5_sp

    print *, cen_mof

  End Subroutine FindCentroid

  subroutine advance_angle(angle,delangle)
    use global_utility_module
    IMPLICIT NONE
    REAL(sp) angle,delangle

    angle=angle+delangle
    if (angle.lt.-MOF_PI) then
      angle=angle+2.0_sp*MOF_PI
    endif
    if (angle.gt.MOF_PI) then
      angle=angle-2.0_sp*MOF_PI
    endif

    return
  end subroutine advance_angle

End Module MODMOF

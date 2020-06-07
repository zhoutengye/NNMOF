Program main
    Use Decision_Tree
    Use MOF_ML
    Use geometry_intersect_module
    Use MOF_routines_module


    Implicit None

    Real(8) :: time1, time2, time3
    Integer :: i

    Call Initialize_Data

    Call initmof( &
       order_algorithm_in, &
       nmat,MOFITERMAX_in, &
       MOF_DEBUG_RECON_in, &
       MOF_TURN_OFF_LS_in, &
       nthreads, &
       nmax_in)

        mofdata(1) = 0.125d0
        mofdata(2) = 1.d0/6.d0 - 0.5d0
        mofdata(3) = 1.d0/6.d0 - 0.5d0
        mofdata(4) = 1.0

    mofdata(8)  = 0.875d0
    mofdata(9)  = ( 0.5d0 - 1.d0 / 48.d0 ) / (7.d0 / 8.d0) - 0.5d0
    mofdata(10) = ( 0.5d0 - 1.d0 / 48.d0 ) / (7.d0 / 8.d0) - 0.5d0 
    mofdata(11) = 2.0
    ! 0.99995076606427225       0.99995076606427225       -2.4076418992368088E-005  
    ! -2.4076418992368088E-005   2.4569222751200787E-005   2.4569222751200787E-005


    ! Call Read_Decision_Tree

    ! mofdata(1) =  0.99995076606427225       
    ! mofdata(2) =  -2.4076418992368088E-005 
    ! mofdata(3) =  2.4569222751200787E-005
    ! mofdata(8) = 1 - mofdata(1)
    ! mofdata(9) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(2)+0.5d0) ) / mofdata(8) - 0.5d0
    ! mofdata(10) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(3)+0.5d0) ) / mofdata(8) - 0.5d0

    ! mofdata(1:3) = [0.68641093,  0.03268939, -0.15239356]
    ! mofdata(8) = 1 - mofdata(1)
    ! mofdata(9) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(2)+0.5d0) ) / mofdata(8) - 0.5d0
    ! mofdata(10) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(3)+0.5d0) ) / mofdata(8) - 0.5d0

    ! Call cpu_time(time1)
    ! DO i = 1,1000000
    !     mofdata(5:7) = Predict_Decision_Tree(mofdata(1:3))
    ! end do
    ! call cpu_time(time2)
    ! write(*,*) 'continuous_mof = 8'
    ! write(*,*) mofdata(5:7)
    ! write(*,*) time2 - time1

    ! DO i = 1,1000000
    ! Call multimaterial_MOF( &j
    !     bfact,dx,xsten0,nhalf0, &
    !     mof_verbose, &
    !     use_ls_data, &
    !     LS_stencil, &
    !     xtetlist_vof, &
    !     xtetlist_cen, &
    !     nmax, &
    !     mofdata, &
    !     multi_centroidA, &
    !     continuous_mof, &
    !     levelrz,nmat,sdim, &
    !     ngeom_recon, &
    !     caller_id)
    ! enddo
    ! call cpu_time(time3)
    ! write(*,*) time3 - time2
    ! write(*,*) 'continuous_mof = 0'
    ! write(*,*) mofdata(5:7)

    ! MOF_TURN_OFF_LS_in = 0
    ! use_ls_data    = 1
    ! continuous_mof = 6

    ! Call multimaterial_MOF( &
    !     bfact,dx,xsten0,nhalf0, &
    !     mof_verbose, &
    !     use_ls_data, &
    !     LS_stencil, &
    !     xtetlist_vof, &
    !     xtetlist_cen, &
    !     nmax, &
    !     mofdata, &
    !     multi_centroidA, &
    !     continuous_mof, &
    !     levelrz,nmat,sdim, &
    !     ngeom_recon, &
    !     caller_id)

    ! write(*,*) 'continuous_mof = 6'
    ! write(*,*) mofdata(5:7)

    ! MOF_TURN_OFF_LS_in = 0
    ! use_ls_data    = 1
    ! continuous_mof = 7

    ! Call multimaterial_MOF( &
    !     bfact,dx,xsten0,nhalf0, &2001 BELLEVUE WAY APT 62
    !     mof_verbose, &
    !     use_ls_data, &
    !     LS_stencil, &
    !     xtetlist_vof, &
    !     xtetlist_cen, &
    !     nmax, &
    !     mofdata, &
    !     multi_centroidA, &
    !     continuous_mof, &
    !     levelrz,nmat,sdim, &
    !     ngeom_recon, &
    !     caller_id)

    ! write(*,*) 'continuous_mof = 7'
    ! write(*,*) mofdata(5:7)

    ! mofdata(1:3)  = 0.d0
    ! mofdata(8:10) = 0.d0

    ! Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata)

    ! write(*,*) mofdata

    ! mofdata(5) = -3.0000000000000027E-002
    ! mofdata(6) = 0.99954989870441180
    ! mofdata(7) = 0.51000000000000001 
    ! mofdata(4) = 1.0

    ! Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata)

    ! write(*,*) mofdata

    ! Call Compare_MOF_LINEAR_2D

    ! Call Gen_Linear_training_2D(1)
!
    ! Call Gen_MOF_training_2D

    ! Call Gen_Linear_Test_2D(2)
! 
    ! Call system(adjustl(trim("python ML_Python/SK_REG_MOF.py")))

    ! mofdata(1:3) = [0.68641093,  0.03268939, -0.15239356]
    ! mofdata(8) = 1 - mofdata(1)
    ! mofdata(9) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(2)+0.5d0) ) / mofdata(8) - 0.5d0
    ! mofdata(10) = ( 0.5d0 - 1.d0 * mofdata(1) * (mofdata(3)+0.5d0) ) / mofdata(8) - 0.5d0
    
    Call Read_Decision_Tree

    ! Call multimaterial_MOF( &
    !     bfact,dx,xsten0,nhalf0, &
    !     mof_verbose, &
    !     use_ls_data, &
    !     LS_stencil, &
    !     xtetlist_vof, &
    !     xtetlist_cen, &
    !     nmax, &
    !     mofdata, &
    !     multi_centroidA, &
    !     continuous_mof, &
    !     levelrz,nmat,sdim, &
    !     ngeom_recon, &
    !     caller_id)

    ! write(*,*) mofdata(5:7)

    ! mofdata(5:7) = Predict_Decision_Tree(mofdata(1:3))

    ! write(*,*) mofdata(5:7)

    Call Compare_DT_MOF

    Return

  End Program main




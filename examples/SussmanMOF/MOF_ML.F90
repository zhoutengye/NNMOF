!================================================================================
!
! Machine Learning module for MOF
!
! 
!
!
!
!
!
!=================================================================================
! 
! By Zhouteng Ye at FSU
! Last update: August 24th, 2018
!
!=================================================================================

Module MOF_ML

    Use probcommon_module
    Implicit None
    Save

    Integer :: caller_id
    Integer :: bfact
    Integer :: nhalf0
    Integer :: use_ls_data 
    Integer :: mof_verbose
    Integer :: nmat
    Integer :: sdim
    Integer :: nmax
    Integer :: continuous_mof

! *** These two variables are defined in module probcommon_module
    ! Integer :: levelrz
    ! Integer :: ngeom_recon

    Integer :: order_algorithm_in     ! all components equal to nmat+1
    Integer :: MOFITERMAX_in 
    Integer :: MOF_DEBUG_RECON_in 
    Integer :: nthreads 
    Integer :: MOF_TURN_OFF_LS_in 
    Integer :: nmax_in

    Real(8), Allocatable :: LS_stencil(:,:,:) 

    Real(8), Allocatable :: dx(:)
             
    Real(8), Allocatable :: xtetlist_vof(:,:,:)
    Real(8), Allocatable :: xtetlist_cen(:,:,:)
    Real(8), Allocatable :: xsten0(:,:)
    Real(8), Allocatable :: mofdata(:)
    Real(8), Allocatable :: multi_centroidA(:,:)


    Integer :: Linear_recon_FLAG !  1: succesfull reconstruction 0: invalid line information 

    
Contains

!---------------------------------------------------------------------------------
!
! Subroutine for initializing data 
!
!    Most parameter can be obtained from the namelist file, then the array  
! varaibles  are allocated and initialized.
!
!---------------------------------------------------------------------------------
!
! Update: August 13th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Initialize_Data

        Implicit None

        Namelist /Inputs_1/ caller_id, bfact, nhalf0, use_ls_data , mof_verbose, & 
                      &  nmat, sdim, nmax, ngeom_recon, continuous_mof, levelrz

        Namelist /MOFINITIAL/ order_algorithm_in , MOFITERMAX_in , & 
                      & MOF_DEBUG_RECON_in, nthreads, MOF_TURN_OFF_LS_in , nmax_in  

        Open(1, file='parameters.namelist')
        Read(1, nml = Inputs_1)
        Read(1, nml = MOFINITIAL)
        ngeom_recon = 2 * sdim + 3
        num_materials = nmat
        nmax_in = nmax


        Allocate(LS_stencil(-1:1,-1:1,nmat))
        Allocate(dx(sdim))
        Allocate(xtetlist_vof(sdim+1,sdim,nmax))
        Allocate(xtetlist_cen(sdim+1,sdim,nmax)) 
        Allocate(xsten0(-nhalf0:nhalf0,sdim))
        Allocate(mofdata(nmat*ngeom_recon))
        Allocate(multi_centroidA(nmat,sdim))


        LS_stencil      = 0.d0
        dx              = 1.d0       
        xtetlist_vof    = 0.d0
        xtetlist_cen    = 0.d0
        xsten0          = 0.d0
        mofdata         = 1.d0
        multi_centroidA = 0.d0

        xsten0(-3,1) = -1.5d0
        xsten0(-2,1) = -1.0d0
        xsten0(-1,1) = -0.5d0
        xsten0(0,1)  = 0.d0
        xsten0(1,1)  = 0.5d0
        xsten0(2,1)  = 1.0d0
        xsten0(3,1)  = 1.5d0

        xsten0(-3,2) = -1.5d0
        xsten0(-2,2) = -1.0d0
        xsten0(-1,2) = -0.5d0
        xsten0(0,2)  = 0.d0
        xsten0(1,2)  = 0.5d0
        xsten0(2,2)  = 1.0d0
        xsten0(3,2)  = 1.5d0

    End Subroutine Initialize_Data

!---------------------------------------------------------------------------------
!
! Subroutine for calculating centroid and volume fraction from intercepts and slopes
!
!     Calculate the centroid and volume fraction from the intercepts and slopes of a 
! straight line. 
!
!     When nmat = 1, only one fraction will be calculated
!     When nmat = 2, both of the fractions will be calculated
!
!     Only the volume fraction and centroid will be updated.
!
!---------------------------------------------------------------------------------
! Six possible positions
!
! ---------------    ---------------    ---------------   ---------------
! |             |    |             |    |             |   |             |
! |             |    |          ****    ***           |   |             |
! |             *    |      ****   |    |  ****       |   *             |
! |           * |    |  ****       |    |      ****   |   | *           |
! |         *   |    ***           |    |          ****   |   *         |
! |       *     |    |             |    |             |   |     *       |
! ------*--------    ---------------    ---------------   --------*------
!     case 1              case 2             case 2           case 4     
!
! --------*------    ----*----------    ----------*----   ------*--------
! |     *       |    |    *        |    |        *    |   |       *     |
! |   *         |    |     *       |    |       *     |   |         *   |
! | *           |    |      *      |    |      *      |   |           * |
! *             |    |       *     |    |     *       |   |             *
! |             |    |        *    |    |    *        |   |             |
! |             |              *   |    |   *         |   |             |
! ---------------    -----------*---    ---*-----------   ---------------
!    case 5               case 3             case 3           case 6     
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 13th, 2018 
!
!---------------------------------------------------------------------------------
    Subroutine Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata)

        Implicit None

        Integer , Intent(In) :: nmat
        Integer , Intent(In) :: ngeom_recon
        Real(8) , Dimension(2*ngeom_recon), Intent(InOut) :: mofdata

        Real(8) :: a,b,c
        Real(8) :: fx0, fx1, f0y,f1y
        Real(8) :: area, area2
        Real(8) , dimension(2) :: centroid
        Integer :: n
        
        ! Those are parameters used for the equation of the straight line
        a = mofdata(5)
        b = mofdata(6)
        c = mofdata(7) * dsqrt(a*a+b*b)

        fx0 =   ( 0.5d0*b - c ) / a
        fx1 = - ( c + 0.5d0*b ) / a
        f0y =   ( 0.5d0*a - c ) / b
        f1y = - ( c + 0.5d0*a ) / b

        Linear_recon_FLAG = 1


        ! Six cases corresponds with 6 possible position of the line.
        If ( fx0.ge.-0.5d0 .and. fx0 .le. 0.5d0 .and. f1y .ge. -0.5d0 .and. f1y .le. 0.5d0) Then
            n = 1
        Elseif (f0y .ge. -0.5d0 .and. f0y .le. 0.5d0 .and. f1y .ge. -0.5d0 .and. f1y .le. 0.5d0) Then
            n = 2
        Elseif (fx0 .ge. -0.5d0 .and. fx0 .le. 0.5d0 .and. fx1 .ge. -0.5d0 .and. fx1 .le. 0.5d0) Then
            n = 3
        Elseif (f0y .ge. -0.5d0 .and. f0y .le. 0.5d0 .and. fx1 .ge. -0.5d0 .and. fx1 .le. 0.5d0) Then
            n = 4
        Elseif (fx0 .ge. -0.5d0 .and. fx0 .le. 0.5d0 .and. f0y .ge. -0.5d0 .and. f0y .le. 0.5d0) Then
            n = 5
        Elseif (fx1 .ge. -0.5d0 .and. fx1 .le. 0.5d0 .and. f1y .ge. -0.5d0 .and. f1y .le. 0.5d0) Then
            n = 6
        Else
            Linear_recon_FLAG = 0
            return
        EndIf

        ! Calculate the centroid and area of the polygons devided by the line 
        ! The position of the line in each case corresponds shown before this subroutine
        Select Case(n) 

        Case (1) 

            area = 0.5d0 * (0.5d0-fx0) * (f1y+0.5d0)
            centroid(1) = 0.5d0 - (0.5d0-fx0+1d-15) / 3.d0
            centroid(2) = (f1y+0.5d0) / 3.d0 - 0.5d0
            If (mofdata(5)<0) Then
                area2 = 1.d0 - area
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf
      
        Case (2)

            area = 0.5d0 * (f0y+f1y+1.d0)
            centroid(1) = 1.d0 / 3.d0 * (2*(f1y+0.5d0)+(f0y+0.5d0)) / (f1y+f0y+1.d0+1.d-15)-0.5d0
            centroid(2) = ((f0y+0.5d0)**2.d0+(f1y+0.5d0)**2.d0+(f0y+0.5d0)*(f1y+0.5d0)) / (f0y+f1y+1.d0+1d-15) / 3.d0 - 0.5d0
            If (mofdata(6)>0) Then
                area2 = 1 - area
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf

        Case (3)

            area = 0.5d0 * (fx0+fx1+1.d0)
            centroid(1) = ((fx0+0.5d0)**2.d0+(fx1+0.5d0)**2.d0+(fx0+0.5d0)*(fx1+0.5d0)) / (fx0+fx1+1.d0+1.d-15) / 3.d0 - 0.5d0
            centroid(2) = 1.d0 / 3.d0 * (2.d0*(fx1+0.5d0)+fx0+0.5d0) / (fx0+fx1+1+1d-15)-0.5d0
            If (mofdata(5) .gt. 0) Then
                area2 = 1.d0 - area;
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf

        Case (4)

            area = 0.5d0 * (0.5d0-f0y) * (fx1+0.5d0)
            centroid(1) = (fx1+0.5d0) / 3.d0 - 0.5d0
            centroid(2) = 0.5d0 - (0.5d0-f0y) / 3.d0
            If (mofdata(5) .gt. 0) Then
                area2 = 1.d0 - area;
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf

        Case (5)

            area = 0.5d0 * (f0y+0.5d0) * (fx0+0.5d0)
            centroid(1) = (fx0+0.5d0) / 3.d0 - 0.5d0
            centroid(2) = (f0y+0.5d0) / 3.d0 - 0.5d0
            If (mofdata(5) .gt. 0) Then
                area2= 1.d0 - area
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf

        Case (6)

            area = 0.5d0 * (0.5d0-f1y) * (0.5d0-fx1)
            centroid(1) = 0.5d0 - (0.5d0-fx1) / 3.d0
            centroid(2) = 0.5d0 - (0.5d0-f1y) / 3.d0
            If (mofdata(5) .lt. 0) Then
                area2 = 1 - area
                centroid(1) = -centroid(1)*area / (area2+1d-15)
                centroid(2) = -centroid(2)*area / (area2+1d-15)
                area = area2
            EndIf
         
        End Select

        ! Assign the first area (volume fraction) and centroid to mofdata
        mofdata(1) = area
        mofdata(2) = centroid(1)
        mofdata(3) = centroid(2)

        ! If 2 materials, assign the second area (volume fraction) and centroid to mofdata
        If(nmat .eq. 2) Then
            mofdata(8)  = 1.d0 - area
            mofdata(9)  = - mofdata(2) * area / (mofdata(8)+1d-15)
            mofdata(10) = - mofdata(3) * area / (mofdata(8)+1d-15)
        EndIf

        ! When the volume fraction is very small, reconstruction skipped
        If (mofdata(1).lt.1.d-4 .or. mofdata(8).lt.1.d-4) Linear_recon_FLAG = 0
        ! If (dabs(mofdata(5)).lt.3.d-1 .or. dabs(mofdata(6)).lt.3.d-1) Linear_recon_FLAG = 0
      
    End Subroutine Linear_Intercepts_Centroid_2D

!---------------------------------------------------------------------------------
!
! Subroutine of the generation of the trainning data from 
!          subroutine "Linear_Intercepts_Centroid_2D"
!
! Method:
!       1: Brute-force generation
!       2: Random generation
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 19th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Gen_Linear_training_2D(Method)

        Implicit None

        Integer :: norm_number = 499
        Integer :: intercepts_number = 201
        Real(8) :: intercepts_U   =  0.75d0
        Real(8) :: intercepts_L   = -0.75d0
        Real(8) :: norm_interval
        Real(8) :: intercepts_interval

        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_linear
        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_MOF

        Integer :: i,j,k

        Integer :: data_amouont = 0
        Integer :: rand_data_amount = 1d6

        Integer :: method

        Real(8) :: norm_x, norm_y, intercept

        Character(len=80) :: File_name

        Real(8) :: PI, angle

        File_name = "DATA/Linear_Training_2D.dat"
        201 Format(7F12.8)

        Open(21,file=file_name,status="Unknown")

        norm_interval = 2.d0 / dble(norm_number+1)
        intercepts_interval = (intercepts_U - intercepts_L) / dble(intercepts_number-1)

        PI = dacos(0.d0) * 2.d0

        write(*,*)

        If (method .eq. 1) Then
        ! Brute-force generation
            Do i = 1, norm_number
            Do k = 1, intercepts_number

                mofdata_linear = 0.d0

                angle = - PI + 2.d0*PI*dble(i)/dble(norm_number+1)
                norm_x = dsin(angle)
                norm_y = dcos(angle)

                mofdata_linear(5) = norm_x
                mofdata_linear(6) = norm_y
                mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(21,201) mofdata_linear(1:7)

            EndDo
            EndDo

        ElseIf (method .eq. 2) Then
        ! random_generation
            Do i = 1, rand_data_amount / 2

                Call random_number(norm_x)
                norm_x = ( norm_x - 0.5d0 ) * 2.d0 
                Call random_number(intercept)
                intercept = ( intercept - 0.5d0 ) * 20.d0

                mofdata_linear(5) = norm_x 
                mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
                mofdata_linear(7) = intercept

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(21,201) mofdata_linear(1:7)

            EndDo

            Do i = 1, rand_data_amount / 2

                Call random_number(norm_x)
                norm_x = ( norm_x - 0.5d0 ) * 2.d0 
                Call random_number(intercept)
                intercept = ( intercept - 0.5d0 ) * 20.d0

                mofdata_linear(6) = norm_x 
                mofdata_linear(5) = dsqrt( 1.d0 -  mofdata_linear(6) * mofdata_linear(6) )
                mofdata_linear(7) = intercept

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(21,201) mofdata_linear(1:7)

            EndDo

        EndIf


    End Subroutine Gen_Linear_training_2D

!---------------------------------------------------------------------------------
!
! Subroutine of the generation of the trainning data from 
!          subroutine "multimaterial_MOF"
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 24th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Gen_MOF_training_2D(method)

        Use geometry_intersect_module
        Use MOF_routines_module

        Implicit None

        Integer :: norm_number = 99
        Integer :: intercepts_number = 101
        Real(8) :: intercepts_U   =  10.d0
        Real(8) :: intercepts_L   = -10.d0
        Real(8) :: norm_interval
        Real(8) :: intercepts_interval
        Integer :: method

        Real(8) :: norm_x, norm_y, intercept
        Integer :: rand_data_amount = 1d6
        Integer :: data_amouont

        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_linear
        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_MOF

        Integer :: i,j,k

        Character(len=80) :: file_name

        File_name = "DATA/MOF_Training_2D.dat"

        Open(21,file=file_name,status="Unknown")
        201 Format(7F12.8)

        norm_interval = 1.d0 / dble(norm_number+1)
        intercepts_interval = (intercepts_U - intercepts_L) / dble(intercepts_number-1)

        If (method .eq. 1) Then

            Do i = 1, norm_number
            Do k = 1, intercepts_number

                mofdata_linear = 0.d0

                mofdata_linear(5) = -1.d0 + norm_interval*dble(i) 
                mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
                mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                mofdata_MOF = mofdata_linear

                Call multimaterial_MOF( &
                    bfact,dx,xsten0,nhalf0, &
                    mof_verbose, &
                    use_ls_data, &
                    LS_stencil, &
                    xtetlist_vof, &
                    xtetlist_cen, &
                    nmax, &
                    mofdata_MOF, &
                    multi_centroidA, &
                    continuous_mof, &
                    levelrz,nmat,sdim, &
                    ngeom_recon, &
                    caller_id)

                Write(21,201) mofdata_MOF(1:7)

            EndDo
            EndDo

            Do i = 1, norm_number
            Do k = 1, intercepts_number

                mofdata_linear = 0.d0

                mofdata_linear(5) = - 1.d0 + norm_interval*dble(i) 
                mofdata_linear(6) = - dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
                mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                mofdata_MOF = mofdata_linear

                Call multimaterial_MOF( &
                    bfact,dx,xsten0,nhalf0, &
                    mof_verbose, &
                    use_ls_data, &
                    LS_stencil, &
                    xtetlist_vof, &
                    xtetlist_cen, &
                    nmax, &
                    mofdata_MOF, &
                    multi_centroidA, &
                    continuous_mof, &
                    levelrz,nmat,sdim, &
                    ngeom_recon, &
                    caller_id)

                Write(21,201) mofdata_MOF(1:7)

            EndDo
            EndDo

        ElseIf (method .eq. 2) Then
        ! random_generation
            Do i = 1, rand_data_amount

                Call random_number(norm_x)
                norm_x = ( norm_x - 0.5d0 ) * 2.d0 
                Call random_number(intercept)
                intercept = ( intercept - 0.5d0 ) * 20.d0

                mofdata_linear(5) = norm_x 
                mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
                mofdata_linear(7) = intercept

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                mofdata_MOF = mofdata_linear

                Call multimaterial_MOF( &
                    bfact,dx,xsten0,nhalf0, &
                    mof_verbose, &
                    use_ls_data, &
                    LS_stencil, &
                    xtetlist_vof, &
                    xtetlist_cen, &
                    nmax, &
                    mofdata_MOF, &
                    multi_centroidA, &
                    continuous_mof, &
                    levelrz,nmat,sdim, &
                    ngeom_recon, &
                    caller_id)

                Write(21,201) mofdata_MOF(1:7)

            EndDo

        EndIf

    End Subroutine Gen_MOF_training_2D

!---------------------------------------------------------------------------------
!
! Compare between MOF reconstruction and linear reconstruction
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 14th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Compare_MOF_LINEAR_2D

        Use geometry_intersect_module
        Use MOF_routines_module

        Implicit None

        Integer :: norm_number = 99
        Integer :: intercepts_number = 101
        Real(8) :: intercepts_U   =  0.75d0
        Real(8) :: intercepts_L   = -0.75d0
        Real(8) :: norm_interval
        Real(8) :: intercepts_interval

        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_linear
        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_MOF

        Real(8) , Dimension(3) :: L1_Error = 0.d0
        Real(8) , Dimension(3) :: L2_Error = 0.d0

        Real(8) , Dimension(3,7) :: L1_Error_max = 0.d0
        Real(8) , Dimension(3,7) :: L2_Error_max = 0.d0

        Real(8) :: error1, error2, error3, error4, error5, error6

        Integer :: i,j,k

        Integer :: data_amouont = 0

        norm_interval = 1.d0 / dble(norm_number+1)
        intercepts_interval = (intercepts_U - intercepts_L) / dble(intercepts_number-1)

        Do i = 1, norm_number
        Do k = 1, intercepts_number

            mofdata_linear = 0.d0

            mofdata_linear(5) = -1.d0 + norm_interval*dble(i) 
            mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
            mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

            If(nmat .eq.2) Then
                mofdata_linear(12) =   mofdata_linear(5)
                mofdata_linear(13) =   mofdata_linear(6)
                mofdata_linear(14) = - mofdata_linear(7)
            EndIf

            Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

            If(Linear_recon_FLAG .eq. 0) cycle

            data_amouont = data_amouont + 1

            mofdata_MOF = mofdata_linear

            Call multimaterial_MOF( &
                bfact,dx,xsten0,nhalf0, &
                mof_verbose, &
                use_ls_data, &
                LS_stencil, &
                xtetlist_vof, &
                xtetlist_cen, &
                nmax, &
                mofdata_MOF, &
                multi_centroidA, &
                continuous_mof, &
                levelrz,nmat,sdim, &
                ngeom_recon, &
                caller_id)

            error1 = dabs(mofdata_MOF(5)-mofdata_linear(5))
            error2 = dabs(mofdata_MOF(6)-mofdata_linear(6))
            error3 = dabs(mofdata_MOF(7)-mofdata_linear(7))
            error4 = (mofdata_MOF(5)-mofdata_linear(5)) * (mofdata_MOF(5)-mofdata_linear(5))
            error5 = (mofdata_MOF(6)-mofdata_linear(6)) * (mofdata_MOF(6)-mofdata_linear(6))
            error6 = (mofdata_MOF(7)-mofdata_linear(7)) * (mofdata_MOF(7)-mofdata_linear(7))

            L1_Error(1) = L1_Error(1) + error1
            L1_Error(2) = L1_Error(2) + error2
            L1_Error(3) = L1_Error(3) + error3

            L2_Error(1) = L2_Error(1) + error4
            L2_Error(2) = L2_Error(2) + error5
            L2_Error(3) = L2_Error(3) + error6

            If ( L1_Error_max(1,1) .lt. error1 ) Then
                L1_Error_max(1,1) = error1
                L1_Error_max(1,2) = mofdata_MOF(1)
                L1_Error_max(1,3) = mofdata_Linear(1)
                L1_Error_max(1,4) = mofdata_MOF(2)
                L1_Error_max(1,5) = mofdata_Linear(2)
                L1_Error_max(1,6) = mofdata_MOF(3)
                L1_Error_max(1,7) = mofdata_Linear(3)
            EndIf

            If ( L1_Error_max(2,1) .lt. error2 ) Then
                L1_Error_max(2,1) = error2
                L1_Error_max(2,2) = mofdata_MOF(5)
                L1_Error_max(2,3) = mofdata_linear(5)
                L1_Error_max(2,4) = mofdata_MOF(6)
                L1_Error_max(2,5) = mofdata_linear(6)
                L1_Error_max(2,6) = mofdata_MOF(7)
                L1_Error_max(2,7) = mofdata_linear(7)
            EndIf

            If ( L1_Error_max(3,1) .lt. error3 ) Then
                L1_Error_max(3,1) = error3
                L1_Error_max(3,2) = mofdata_MOF(1)
                L1_Error_max(3,3) = mofdata_linear(1)
                L1_Error_max(3,4) = mofdata_MOF(2)
                L1_Error_max(3,5) = mofdata_linear(2)
                L1_Error_max(3,6) = mofdata_MOF(3)
                L1_Error_max(3,7) = mofdata_linear(3)
            EndIf

            If ( L2_Error_max(1,1) .lt. error4 ) Then
                L2_Error_max(1,1) = error4
                L2_Error_max(1,2) = mofdata_MOF(5)
                L2_Error_max(1,3) = mofdata_linear(5)
                L2_Error_max(1,4) = mofdata_MOF(6)
                L2_Error_max(1,5) = mofdata_linear(6)
                L2_Error_max(1,6) = mofdata_MOF(7)
                L2_Error_max(1,7) = mofdata_linear(7)
            EndIf

            If ( L2_Error_max(2,1) .lt. error5 ) Then
                L2_Error_max(2,1) = error5
                L2_Error_max(2,2) = mofdata_MOF(5)
                L2_Error_max(2,3) = mofdata_linear(5)
                L2_Error_max(2,4) = mofdata_MOF(6)
                L2_Error_max(2,5) = mofdata_linear(6)
                L2_Error_max(2,6) = mofdata_MOF(7)
                L2_Error_max(2,7) = mofdata_linear(7)
            EndIf

            If ( L2_Error_max(3,1) .lt. error6 ) Then
                L2_Error_max(3,1) = error6
                L2_Error_max(3,2) = mofdata_MOF(5)
                L2_Error_max(3,3) = mofdata_linear(5)
                L2_Error_max(3,4) = mofdata_MOF(6)
                L2_Error_max(3,5) = mofdata_linear(6)
                L2_Error_max(3,6) = mofdata_MOF(7)
                L2_Error_max(3,7) = mofdata_linear(7)
            EndIf



        EndDo
        EndDo

        L1_Error = L1_Error / dble(data_amouont)
        L2_Error = L2_Error / dble(data_amouont)

        write(*,*) L1_Error

        write(*,*) L2_Error

        write(*,*) L1_Error_max(1,:)
        ! write(*,*) L1_Error_max(2,:)
        write(*,*) L1_Error_max(3,:)

        ! write(*,*) L2_Error_max(1,:)
        ! write(*,*) L2_Error_max(2,:)
        ! write(*,*) L2_Error_max(3,:)


    End Subroutine Compare_MOF_LINEAR_2D

! ---------------------------------------------------------------------------------

! Trainning the data

! Input: slope, intercepts 

! Output: 


! ---------------------------------------------------------------------------------

! By Zhouteng Ye
! Update: August 17th, 2018 
! ---------------------------------------------------------------------------------
    Subroutine Trainning_2D(Method)

        Implicit None

        Integer, Intent(in) :: Method

        Select Case(method)

        Case (1)
            Write(*,*) 'Linear regression trainning'

        Case (2)
            Write(*,*) 'Neuron network training'
            Call system(adjustl(trim("python ML_Python/NN_MOF.py")))

        Case (3)
            Write(*,*) 'Nearest neighbour training'

        Case (4)
            Write(*,*) 'Decision Tree training'

        Case (5)
            Write(*,*) 'Random forest reaining'

        End Select 

    End Subroutine Trainning_2D

!---------------------------------------------------------------------------------
!
! Subroutine of the generation of the trainning data from 
!          subroutine "Linear_Intercepts_Centroid_2D"
!
! When 
!      method = 1: Brute-force generation with norm_y and intercept
!      method = 2: Random data
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 24th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Gen_Linear_Test_2D(method)

        Implicit None

        Integer :: norm_number = 199
        Integer :: intercepts_number = 101
        Real(8) :: intercepts_U   =  0.75d0
        Real(8) :: intercepts_L   = -0.75d0
        Real(8) :: norm_interval
        Real(8) :: intercepts_interval

        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_linear
        Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_MOF

        Integer :: i,j,k

        Integer :: data_amouont = 0
        Integer :: rand_data_amount = 1d6

        Real(8) :: norm_x, norm_y, intercept
        Integer :: method

        Character(len=80) :: file_name

        File_name = "DATA/Linear_Test_2D.dat"
        201 Format(7F12.8)

        Open(23,file=file_name,status="Unknown")

        norm_interval = 2.d0 / dble(norm_number+1)
        intercepts_interval = (intercepts_U - intercepts_L) / dble(intercepts_number-1)


        If (method.eq.1) Then
        ! brute-forth generation
            Do i = 1, norm_number
            Do k = 1, intercepts_number

                mofdata_linear = 0.d0

                mofdata_linear(6) = -1.d0 + norm_interval*dble(i) 
                mofdata_linear(5) = dsqrt( 1.d0 -  mofdata_linear(6) * mofdata_linear(6) )
                mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(23,201) mofdata_linear(1:7)
            EndDo
            EndDo

            ! Do i = 1, norm_number
            ! Do k = 1, intercepts_number

            !     mofdata_linear = 0.d0

            !     mofdata_linear(6) = -1.d0 + norm_interval*dble(i) 
            !     mofdata_linear(5) = - dsqrt( 1.d0 -  mofdata_linear(6) * mofdata_linear(6) )
            !     mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

            !     If(nmat .eq.2) Then
            !         mofdata_linear(12) =   mofdata_linear(5)
            !         mofdata_linear(13) =   mofdata_linear(6)
            !         mofdata_linear(14) = - mofdata_linear(7)
            !     EndIf

            !     Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

            !     If(Linear_recon_FLAG .eq. 0) cycle

            !     data_amouont = data_amouont + 1

            !     Write(23,201) mofdata_linear(1:7)
            ! EndDo
            ! EndDo

        ElseIf (method.eq.2) Then
        ! random_generation
            Do i = 1, rand_data_amount / 2

                Call random_number(norm_x)
                norm_x = ( norm_x - 0.5d0 ) * 2.d0 
                Call random_number(intercept)
                intercept = ( intercept - 0.5d0 ) * 1.5d0

                mofdata_linear(6) = norm_x 
                mofdata_linear(5) = dsqrt( 1.d0 -  mofdata_linear(6) * mofdata_linear(6) )
                mofdata_linear(7) = intercept

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(23,201) mofdata_linear(1:7)

            EndDo

            Do i = 1, rand_data_amount / 2

                Call random_number(norm_x)
                norm_x = ( norm_x - 0.5d0 ) * 2.d0 
                Call random_number(intercept)
                intercept = ( intercept - 0.5d0 ) * 20.d0

                mofdata_linear(5) = norm_x 
                mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(6) )
                mofdata_linear(7) = intercept

                If(nmat .eq.2) Then
                    mofdata_linear(12) =   mofdata_linear(5)
                    mofdata_linear(13) =   mofdata_linear(6)
                    mofdata_linear(14) = - mofdata_linear(7)
                EndIf

                Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

                If(Linear_recon_FLAG .eq. 0) cycle

                data_amouont = data_amouont + 1

                Write(23,201) mofdata_linear(1:7)

            EndDo

        EndIf


    End Subroutine Gen_Linear_Test_2D

!---------------------------------------------------------------------------------
!
! Subroutine for comparing between decision tree results and MOF results
!
! There are 4 cases:
!
! 1: traditional MOF
!    use_ls_data = 0,
!    continuous_mof
!
! 2: Use decision tree result without conserving mass
! 
!
!
! 3:
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: September 4th, 2018 
!---------------------------------------------------------------------------------
    Subroutine Compare_DT_MOF

        Use Decision_Tree
        Use geometry_intersect_module
        Use MOF_routines_module

        Implicit None
        real(8) :: time1, time2, time3
        Real(8), allocatable :: DT_results(:,:)
        Real(8), allocatable :: MOF_results(:,:)
        Real(8), allocatable :: Input(:,:)
        Real(8), allocatable :: Exact(:,:)
        Real(8), allocatable :: MOF_0(:,:)
        Real(8), allocatable :: MOF_6(:,:)
        Real(8), allocatable :: MOF_7(:,:)
        Real(8), allocatable :: MOF_8(:,:)
        Real(8), allocatable :: MOF_81(:,:)
        Real(8) :: skip_number
        integer i, n
        real(8) :: line
        Real(8) :: time01, time02, time61, time62
        Real(8) :: time71, time72, time81, time82
        Real(8) :: time811,time812
        Real(8) :: L1_0(3) = 0 
        Real(8) :: L1_6(3) = 0
        Real(8) :: L1_7(3) = 0
        Real(8) :: L1_8(3) = 0
        Real(8) :: L1_81(3) = 0
        Real(8) :: L2_0(3) = 0
        Real(8) :: L2_6(3) = 0
        Real(8) :: L2_7(3) = 0
        Real(8) :: L2_8(3) = 0
        Real(8) :: L2_81(3) = 0
        
        ! get the line number of training data
        n=0
        Open(71,file='DATA/Linear_Test_2D.dat',status='old')
        do while (.true.)
            read(71,*,end=100) line
            n=n+1
        enddo
        100 continue
        close(71)

        Allocate(Input(n,3))
        Allocate(Exact(n,3))
        Allocate(DT_results(n,3))
        Allocate(MOF_results(n,3))
        Allocate(MOF_0(n,14))
        Allocate(MOF_6(n,14))
        Allocate(MOF_7(n,14))
        Allocate(MOF_8(n,14))
        Allocate(MOF_81(n,14))

        
        Open(71,file='DATA/Linear_Test_2D.dat',status='old')
        do i = 1,n
            read(71,701) Input(i,:), skip_number, Exact(i,:)
            MOF_0(i,1:3) = Input(i,1:3)
            MOF_0(i,8) = 1 - MOF_0(i,1)
            MOF_0(i,9) = ( 0.5d0 - 1.d0 * MOF_0(i,1) * (MOF_0(i,2)+0.5d0) ) / MOF_0(i,8) - 0.5d0
            MOF_0(i,10) = ( 0.5d0 - 1.d0 * MOF_0(i,1) * (MOF_0(i,3)+0.5d0) ) / MOF_0(i,8) - 0.5d0
        enddo
        701 Format(7F12.8)
        MOF_6 = MOF_0
        MOF_7 = MOF_0
        MOF_8 = MOF_0
        MOF_81 = MOF_0

        write(*,*) "test data length:", n

        MOF_TURN_OFF_LS_in = 1
        use_ls_data    = 0
        continuous_mof = 0

        Call CPU_TIME(time01)
        Do i = 1,n
            mofdata = MOF_0(i,:)
            Call multimaterial_MOF( &
            bfact,dx,xsten0,nhalf0, &
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
            MOF_0(i,5:7) = mofdata(5:7)
        EndDo
        Call CPU_TIME(time02)

        write(*,*) "CPU time for traditional MOF", time02-time01

        MOF_TURN_OFF_LS_in = 0
        use_ls_data    = 1
        continuous_mof = 6

        Call CPU_TIME(time61)
        Do i = 1,n
            mofdata = MOF_6(i,:)
            Call multimaterial_MOF( &
            bfact,dx,xsten0,nhalf0, &
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
            MOF_6(i,5:7) = mofdata(5:7)
        EndDo
        Call CPU_TIME(time62)

        write(*,*) "CPU time for MOF with LS initialization", time62-time61

        MOF_TURN_OFF_LS_in = 0
        use_ls_data    = 1
        continuous_mof = 7

        Call CPU_TIME(time71)
        Do i = 1,n
            mofdata = MOF_7(i,:)
            Call multimaterial_MOF( &
            bfact,dx,xsten0,nhalf0, &
            mof_verbose,       &
            use_ls_data,       &
            LS_stencil,        &
            xtetlist_vof,      &
            xtetlist_cen,      &
            nmax,              &
            mofdata,           &
            multi_centroidA,   &
            continuous_mof,    &
            levelrz,nmat,sdim, &
            ngeom_recon,       &
            caller_id)
            MOF_7(i,5:7) = mofdata(5:7)
        EndDo
        Call CPU_TIME(time72)

        write(*,*) "CPU time for MOF with continuous MOF=7", time72-time71

        MOF_TURN_OFF_LS_in = 0
        use_ls_data    = 1
        continuous_mof = 8

        Call CPU_TIME(time81)
        Do i = 1,n
            mofdata = MOF_8(i,:)
            Call multimaterial_MOF( &
            bfact,dx,xsten0,nhalf0, &
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
            MOF_8(i,5:7) = mofdata(5:7)
        EndDo

        Call CPU_TIME(time82)

        write(*,*) "CPU time for MOF with continuous MOF=8", time82-time81

        MOF_TURN_OFF_LS_in = 0
        use_ls_data    = 1
        continuous_mof = 8

        Call CPU_TIME(time811)
        Do i = 1,n
            MOF_81(i,5:7) = Predict_Decision_Tree(MOF_81(i,1:3))
        EndDo
        Call CPU_TIME(time812)

        Write(*,*) "CPU time for MOF using DT prediction", time812-time811

        Do i = 1,n
            L1_0(1)  = L1_0(1)  + dabs( MOF_0(i,5)  - exact(i,1) )
            L1_0(2)  = L1_0(2)  + dabs( MOF_0(i,6)  - exact(i,2) )
            L1_0(3)  = L1_0(3)  + dabs( MOF_0(i,7)  - exact(i,3) )
            L1_6(1)  = L1_6(1)  + dabs( MOF_6(i,5)  - exact(i,1) )
            L1_6(2)  = L1_6(2)  + dabs( MOF_6(i,6)  - exact(i,2) )
            L1_6(3)  = L1_6(3)  + dabs( MOF_6(i,7)  - exact(i,3) )
            L1_7(1)  = L1_7(1)  + dabs( MOF_7(i,5)  - exact(i,1) )
            L1_7(2)  = L1_7(2)  + dabs( MOF_7(i,6)  - exact(i,2) )
            L1_7(3)  = L1_7(3)  + dabs( MOF_7(i,7)  - exact(i,3) )
            L1_8(1)  = L1_8(1)  + dabs( MOF_8(i,5)  - exact(i,1) )
            L1_8(2)  = L1_8(2)  + dabs( MOF_8(i,6)  - exact(i,2) )
            L1_8(3)  = L1_8(3)  + dabs( MOF_8(i,7)  - exact(i,3) )
            L1_81(1) = L1_81(1) + dabs( MOF_81(i,5) - exact(i,1) )
            L1_81(2) = L1_81(2) + dabs( MOF_81(i,6) - exact(i,2) )
            L1_81(3) = L1_81(3) + dabs( MOF_81(i,7) - exact(i,3) )
        EndDo

        Write(*,*) "Continuous_MOF = 0, L1 Error for x_norm",       L1_0(1) / float(n)
        Write(*,*) "Continuous_MOF = 0, L1 Error for y_norm",       L1_0(2) / float(n)
        Write(*,*) "Continuous_MOF = 0, L1 Error for intercepts",   L1_0(3) / float(n)
        Write(*,*) "Continuous_MOF = 6, L1 Error for x_norm",       L1_6(1) / float(n)
        Write(*,*) "Continuous_MOF = 6, L1 Error for y_norm",       L1_6(2) / float(n)
        Write(*,*) "Continuous_MOF = 6, L1 Error for intercepts",   L1_6(3) / float(n)
        Write(*,*) "Continuous_MOF = 7, L1 Error for x_norm",       L1_7(1) / float(n)
        Write(*,*) "Continuous_MOF = 7, L1 Error for y_norm",       L1_7(2) / float(n)
        Write(*,*) "Continuous_MOF = 7, L1 Error for intercepts",   L1_7(3) / float(n)
        Write(*,*) "Continuous_MOF = 8, L1 Error for x_norm",       L1_8(1) / float(n)
        Write(*,*) "Continuous_MOF = 8, L1 Error for y_norm",       L1_8(2) / float(n)
        Write(*,*) "Continuous_MOF = 8, L1 Error for intercepts",   L1_8(3) / float(n)
        Write(*,*) "Direct decision tree, L1 Error for x_norm",     L1_81(1) / float(n)
        Write(*,*) "Direct decision tree, L1 Error for y_norm",     L1_81(2) / float(n)
        Write(*,*) "Direct decision tree, L1 Error for intercepts", L1_81(3) / float(n)

    End Subroutine Compare_DT_MOF

! !---------------------------------------------------------------------------------
! !
! ! Subroutine of the generation of the trainning data from 
! !          subroutine "multimaterial_MOF"
! !
! !---------------------------------------------------------------------------------
! !
! ! By Zhouteng Ye
! ! Update: August 19th, 2018 
! !---------------------------------------------------------------------------------
!     Subroutine Gen_MOF_Test_2D

!         Use geometry_intersect_module
!         Use MOF_routines_module

!         Implicit None

!         Integer :: norm_number = 199
!         Integer :: intercepts_number = 101
!         Real(8) :: intercepts_U   =  0.75d0
!         Real(8) :: intercepts_L   = -0.75d0
!         Real(8) :: norm_interval
!         Real(8) :: intercepts_interval

!         Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_linear
!         Real(8) , Dimension(nmat*ngeom_recon) :: mofdata_MOF

!         Integer :: i,j,k

!         Character(len=80) :: file_name

!         File_name = "DATA/MOF_Test_2D.dat"

!         Open(24,file=file_name,status="Unknown")
!         201 Format(7F12.8)

!         norm_interval = 1.d0 / dble(norm_number+1)
!         intercepts_interval = (intercepts_U - intercepts_L) / dble(intercepts_number-1)

!         Do i = 1, norm_number
!         Do k = 1, intercepts_number

!             mofdata_linear = 0.d0

!             mofdata_linear(5) = -1.d0 + norm_interval*dble(i) 
!             mofdata_linear(6) = dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
!             mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

!             If(nmat .eq.2) Then
!                 mofdata_linear(12) =   mofdata_linear(5)
!                 mofdata_linear(13) =   mofdata_linear(6)
!                 mofdata_linear(14) = - mofdata_linear(7)
!             EndIf

!             Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

!             If(Linear_recon_FLAG .eq. 0) cycle

!             mofdata_MOF = mofdata_linear

!             Call multimaterial_MOF( &
!                 bfact,dx,xsten0,nhalf0, &
!                 mof_verbose, &
!                 use_ls_data, &
!                 LS_stencil, &
!                 xtetlist_vof, &
!                 xtetlist_cen, &
!                 nmax, &
!                 mofdata_MOF, &
!                 multi_centroidA, &
!                 continuous_mof, &
!                 levelrz,nmat,sdim, &
!                 ngeom_recon, &
!                 caller_id)

!             Write(24,201) mofdata_MOF(1:7)
!         EndDo
!         EndDo

!         Do i = 1, norm_number
!         Do k = 1, intercepts_number

!             mofdata_linear = 0.d0

!             mofdata_linear(5) = - 1.d0 + norm_interval*dble(i) 
!             mofdata_linear(6) = - dsqrt( 1.d0 -  mofdata_linear(5) * mofdata_linear(5) )
!             mofdata_linear(7) = intercepts_L + intercepts_interval*dble(k-1) 

!             If(nmat .eq.2) Then
!                 mofdata_linear(12) =   mofdata_linear(5)
!                 mofdata_linear(13) =   mofdata_linear(6)
!                 mofdata_linear(14) = - mofdata_linear(7)
!             EndIf

!             Call Linear_Intercepts_Centroid_2D(nmat,ngeom_recon,mofdata_linear)

!             If(Linear_recon_FLAG .eq. 0) cycle

!             mofdata_MOF = mofdata_linear

!             Call multimaterial_MOF( &
!                 bfact,dx,xsten0,nhalf0, &
!                 mof_verbose, &
!                 use_ls_data, &
!                 LS_stencil, &
!                 xtetlist_vof, &
!                 xtetlist_cen, &
!                 nmax, &
!                 mofdata_MOF, &
!                 multi_centroidA, &
!                 continuous_mof, &
!                 levelrz,nmat,sdim, &
!                 ngeom_recon, &
!                 caller_id)

!             Write(24,201) mofdata_MOF(1:7)

!         EndDo
!         EndDo

!     End Subroutine Gen_MOF_Test_2D

!---------------------------------------------------------------------------------
!
! Test the data
! Input
!
!
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 13th, 2018 
!---------------------------------------------------------------------------------
    ! Subroutine Test_2D

    !     Implicit None

    ! End Subroutine Test_2D


!---------------------------------------------------------------------------------
!
! Call Python code for 
!
!---------------------------------------------------------------------------------
! 
! Call
! 
!---------------------------------------------------------------------------------
!
! By Zhouteng Ye
! Update: August 24th, 2018 
!---------------------------------------------------------------------------------
    ! Subroutine Decision_Tree_Trainging

    !     Implicit None

    !     Call system(adjustl(trim("python ML_Python/DT_MOF.py")))

    ! End Subroutine Decision_Tree_Training

End Module MOF_ML

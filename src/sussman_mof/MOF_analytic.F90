Program main

  Implicit None

  Integer :: sdim

  Real(8), Allocatable :: delta1(:), delta2(:)
  Real(8) :: volume
  Real(8), Allocatable :: centroid(:)
  Real(8) :: intercepts
  Real(8), Allocatable :: normal(:) 
  Real(8), Allocatable :: X0(:), Xk(:), Xk12(:), Xk1(:)

  Real(8), Allocatable :: dp(:), dx(:)
  Real(8) :: p, ddp

  Integer :: i

  sdim = 2


  Allocate(delta1(sdim))
  Allocate(delta2(sdim))
  Allocate(centroid(sdim))
  Allocate(normal(sdim))
  Allocate(X0(sdim))
  Allocate(Xk(sdim))
  Allocate(Xk12(sdim))
  Allocate(Xk1(sdim))

  Allocate(dp(sdim))
  Allocate(dx(sdim))

  i = 0

  Volume = 1.0d0/8.0d0
  X0(1)  = 1.0d0/18.0d0
  X0(2)  = 1.0d0/4.0d0
  Xk(1)  = dsqrt(Volume*2.0d0/9.0d0)
  Xk(1)  = dsqrt(10.0d0/13.0d0/36.0d0)
  Xk(2)  = Volume * 2.0d0 / 9.0d0 /Xk(1)
  xk1    = 0.0d0

  Volume = 1.0d0/4.0d0
  X0(1)  = 0.4d0
  X0(2)  = 0.1d0
  Xk(1)  = 0.5d0
  ! Xk(1)  = dsqrt(10.0d0/13.0d0/36.0d0)
  Xk(2)  = 1.0/8.0d0
  xk1    = 0.0d0

  ! Volume = 3.0d0/64.0d0
  ! X0(1)  = 1.0d0/6.0d0
  ! X0(2)  = 1.0d0/6.0d0
  ! x0(3)  = 1.d00/6.0d0
  ! Xk(1)  = (Volume*2.0d0/9.0d0)**(1.d0/3.d0)
  ! ! Xk(1)  = dsqrt(10.0d0/13.0d0/36.0d0)
  ! Xk(2)  = xk(1)
  ! xk(3)  = xk(1)
  ! xk1    = 0.0d0



    ! !! 2D Triangle
    ! Do While( (Xk(1) - Xk1(1)) .gt. 1.0d-8 .or. (Xk(2) - Xk1(2)).gt. 1.0d-8 )
    !   p = Xk(1) * Xk(2) - 2.d0/9.d0*Volume
    !   dx = x0 - xk
    !   dp(1) = Xk(2)
    !   dp(2) = Xk(1)
    !   ddp = dp(1) * dp(1) + dp(2) * dp(2)

    !   delta1 = - p * dp / ddp
    !   xk12 = xk + delta1
    !   delta2 = dx - (dx(1)*dp(1)+dx(2)*dp(2))/ddp*dp
    !   xk1  = xk12 + delta2

    !   xk12 = xk1
    !   xk1  = xk
    !   xk   = xk12

    !   i = i+1

    ! End Do

    !! 2D Triangle
    Do While( (Xk(1) - Xk1(1)) .gt. 1.0d-8 .or. (Xk(2) - Xk1(2)).gt. 1.0d-8 )
      p = (0.5d0-Xk1(1))*(0.5d0-Xk1(1))*Volume*6 +Volume / 2.0d0 - Xk(2)
      dx = x0 - Xk
      dp(1) = - 12 * Xk(1) * Volume * (0.5d0-Xk(1))
      dp(2) = - 1.0d0
      ddp = dp(1) * dp(1) + dp(2) * dp(2)

      delta1 = - p * dp / ddp
      xk12 = xk + delta1
      delta2 = dx - (dx(1)*dp(1)+dx(2)*dp(2))/ddp*dp
      xk1  = xk12 + delta2

      xk12 = xk1
      xk1  = xk
      xk   = xk12

      i = i+1

    End Do

    ! !! 3D Triangle
    ! Do While( (Xk(1) - Xk1(1)) .gt. 1.0d-8 .or. (Xk(2) - Xk1(2)).gt. 1.0d-8 .or. (Xk(3) - Xk1(3)).gt. 1.0d-8)
    !   p = Xk(1) * Xk(2) * Xk(3) - 2.d0/9.d0*Volume
    !   dx = x0 - xk
    !   dp(1) = Xk(2)*Xk(3)
    !   dp(2) = Xk(1)*Xk(3)
    !   dp(3) = Xk(1)*Xk(2)
    !   ddp = dp(1) * dp(1) + dp(2) * dp(2) + dp(3) * dp(3)

    ! !! 3D Triangle
    ! Do While( (Xk(1) - Xk1(1)) .gt. 1.0d-8 .or. (Xk(2) - Xk1(2)).gt. 1.0d-8 .or. (Xk(3) - Xk1(3)).gt. 1.0d-8)
    !   p = Xk(1) * Xk(2) * Xk(3) - 2.d0/9.d0*Volume
    !   dx = x0 - xk
    !   dp(1) = Xk(2)*Xk(3)
    !   dp(2) = Xk(1)*Xk(3)
    !   dp(3) = Xk(1)*Xk(2)
    !   ddp = dp(1) * dp(1) + dp(2) * dp(2) + dp(3) * dp(3)

    !   delta1 = - p * dp / ddp
    !   xk12 = xk + delta1
    !   delta2 = dx - (dx(1)*dp(1)+dx(2)*dp(2)+dx(3)*dp(3))/ddp*dp
    !   xk1  = xk12 + delta2

    !   xk12 = xk1
    !   xk1  = xk
    !   xk   = xk12

    !   i = i+1

    ! EndDo
    !   delta1 = - p * dp / ddp
    !   xk12 = xk + delta1
    !   delta2 = dx - (dx(1)*dp(1)+dx(2)*dp(2)+dx(3)*dp(3))/ddp*dp
    !   xk1  = xk12 + delta2

    !   xk12 = xk1
    !   xk1  = xk
    !   xk   = xk12

    !   i = i+1

    ! EndDo


  write(*,*) i, xk, xk(1) * xk(2) * xk(3) *9.d0/2.d0*64.d0/3.d0
  ! write(*,*) x0(1:2)
  ! x12(1:2) = x0(1:2) * x0(1:2)
  ! write(*,*) sum(x0(1:2))


End Program

 
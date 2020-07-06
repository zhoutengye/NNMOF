Module InitSH
  Use ModGlobal
  use ModTools
Contains
  Real(8) Function ShapeZalesak(x,y,z)
    Implicit None
    Real(sp) :: x, y, z
    Real(sp) :: lv1, lv2, lv3
    ! lv1 = 0.2_sp - sqrt(x**2+y**2+z**2)
    ! lv2 = min(0.04_sp-y, 0.04_sp+y)
    ! lv3 = 0.1_sp-x
    ! ShapeZalesak = min( -min(lv3, lv2), lv1)
    Shapezalesak = 0.2_sp - sqrt(x**2+y**2+z**2)

  End Function ShapeZalesak

  Subroutine InitVolume
    Use ModVOFFunc
    Implicit None
    Real(sp), Allocatable :: xc(:)
    Real(sp), Allocatable :: xl(:)
    Real(sp), Allocatable :: ls(:,:,:)
    Type(Octree) :: gr
    Integer :: i,j,k
    Integer :: flag
    Real(sp) :: dx
    Integer :: Octreelevel = 5
    Real(sp) :: xx, yy, zz
    Real(sp) :: n3(3), c3(3)
    Integer :: ng

    ng = nl(1)
    Allocate(ls(0 :ng+1,0:ng+1,0:ng/2+1))

    dx = 1.0_sp / dble(ng)
    Allocate(xc(ng))
    Allocate(xl(ng))
    Do k = 1, ng
      xc(k) = dble(k)/ dble(ng) - 0.5_sp * dx
      xl(k) = dble(k)/ dble(ng) - 1.0_sp * dx
    End Do
    dx = 1.0_sp / dble(ng)
    ShapeLevelSet => ShapeZalesak

    Do k = 1, ng / 2
      Do j = 1, ng
        Do i = 1, ng
          xx = xc(i) - 0.5_sp
          yy = xc(j) - 0.5_sp
          zz = xc(k) - 0.25_sp
          ls(i,j,k) = ShapeZalesak(xx,yy,zz)
        End Do
      End Do
    End Do

    Do k = 1, ng / 2
      Do j = 1, ng
        Do i = 1, ng
          xx = xc(i) - 0.5_sp
          yy = xc(j) - 0.5_sp
          zz = xc(k) - 0.25_sp
          flag = Inout(xx, yy, zz, dx, dx, dx)
          If (flag .eq. 1) Then
            phi(i,j,k) = 1.0
            cx(i,j,k)  = 0.5_sp
            cy(i,j,k)  = 0.5_sp
            cz(i,j,k)  = 0.5_sp
          Elseif (flag .eq. 2) Then
            phi(i,j,k) = 0.0
            cx(i,j,k)  = 0.0_sp
            cy(i,j,k)  = 0.0_sp
            cz(i,j,k)  = 0.0_sp
          Elseif (flag .eq. 3) Then
            gr%level = 1
            gr%maxlevel = Octreelevel
            gr%xc(1) = xx
            gr%xc(2) = yy
            gr%xc(3) = zz
            gr%dx(1) = dx
            gr%dx(2) = dx
            gr%dx(3) = dx
            Call VolumeOctree(gr)
            phi(i,j,k) = gr%vof
            ! Call VolumeCentroidOctree(gr)
            ! if ( phi(i,j,k) .gt. epsc ) cx(i,j,k)  = gr%centroid(1) + 0.5_sp
            ! if ( phi(i,j,k) .gt. epsc ) cy(i,j,k)  = gr%centroid(2) + 0.5_sp
            ! if ( phi(i,j,k) .gt. epsc ) cz(i,j,k)  = gr%centroid(3) + 0.5_sp
            n3(1) = ls(i-1,  j,  k) - ls(i+1,  j,  k)
            n3(2) = ls(  i,j-1,  k) - ls(  i,j+1,  k)
            n3(3) = ls(  i,  j,k-1) - ls(  i,  j,k+1)
            Call Normalization1(n3)
            Call FloodSZ_BackwardC(n3,phi(i,j,k),c3)
            cx(i,j,k) = c3(1)
            cy(i,j,k) = c3(2)
            cz(i,j,k) = c3(3)
          EndIf
        End Do
      End Do
    End Do

  End Subroutine InitVolume

End Module InitSH

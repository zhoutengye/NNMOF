module mod_bound
  use mpi
  use ModGlobal
  implicit none
  private
  Integer :: rp = sp
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,n,bc,isoutflow,dl,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in), dimension(0:1,3) :: isoutflow
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: u,v,w
    integer :: q,idir,sgn,ioutflowdir
    !
    call updthalo((/n(1),n(2)/),1,u)
    call updthalo((/n(1),n(2)/),2,u)
    call updthalo((/n(1),n(2)/),1,v)
    call updthalo((/n(1),n(2)/),2,v)
    call updthalo((/n(1),n(2)/),1,w)
    call updthalo((/n(1),n(2)/),2,w)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1,1),0,n(1),1,.false.,bc(0,1,1),dl(1),u)
      call set_bc(cbc(0,1,2),0,n(1),1,.true. ,bc(0,1,2),dl(1),v)
      call set_bc(cbc(0,1,3),0,n(1),1,.true. ,bc(0,1,3),dl(1),w)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1,1),1,n(1),1,.false.,bc(1,1,1),dl(1),u)
      call set_bc(cbc(1,1,2),1,n(1),1,.true. ,bc(1,1,2),dl(1),v)
      call set_bc(cbc(1,1,3),1,n(1),1,.true. ,bc(1,1,3),dl(1),w)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2,1),0,n(2),2,.true. ,bc(0,2,1),dl(2),u)
      call set_bc(cbc(0,2,2),0,n(2),2,.false.,bc(0,2,2),dl(2),v)
      call set_bc(cbc(0,2,3),0,n(2),2,.true. ,bc(0,2,3),dl(2),w)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2,1),1,n(2),2,.true. ,bc(1,2,1),dl(2),u)
      call set_bc(cbc(1,2,2),1,n(2),2,.false.,bc(1,2,2),dl(2),v)
      call set_bc(cbc(1,2,3),1,n(2),2,.true. ,bc(1,2,3),dl(2),w)
    endif
    call set_bc(cbc(0,3,1),0,n(3),3,.true. ,bc(0,3,1),dzc(0)   ,u)
    call set_bc(cbc(0,3,2),0,n(3),3,.true. ,bc(0,3,2),dzc(0)   ,v)
    call set_bc(cbc(0,3,3),0,n(3),3,.false.,bc(0,3,3),dzf(0)   ,w)
    call set_bc(cbc(1,3,1),1,n(3),3,.true. ,bc(1,3,1),dzc(n(3)),u)
    call set_bc(cbc(1,3,2),1,n(3),3,.true. ,bc(1,3,2),dzc(n(3)),v)
    call set_bc(cbc(1,3,3),1,n(3),3,.false.,bc(1,3,3),dzf(n(3)),w)
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(n,ioutflowdir,dl,dzf,u,v,w)
        endif
      enddo
    enddo
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,n,bc,dl,dzc,dzf,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n 
    real(rp)         , intent(in), dimension(0:1,3) :: bc
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzc,dzf
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    !
    call updthalo((/n(1),n(2)/),1,p)
    call updthalo((/n(1),n(2)/),2,p)
    !
    if(left .eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,1),0,n(1),1,.true.,bc(0,1),dl(1),p)
    endif
    if(right.eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,1),1,n(1),1,.true.,bc(1,1),dl(1),p)
    endif
    if(front.eq.MPI_PROC_NULL) then
      call set_bc(cbc(0,2),0,n(2),2,.true.,bc(0,2),dl(2),p)
     endif
    if(back .eq.MPI_PROC_NULL) then
      call set_bc(cbc(1,2),1,n(2),2,.true.,bc(1,2),dl(2),p)
    endif
    call set_bc(cbc(0,3),0,n(3),3,.true.,bc(0,3),dzc(0)   ,p)
    call set_bc(cbc(1,3),1,n(3),3,.true.,bc(1,3),dzc(n(3)),p)
    return
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,n,idir,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound,n,idir
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and.centered) then
      factor = 2.*factor
      sgn    = -1.
    endif
    if(ctype.eq.'N'.and.centered) then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1.
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !p(0  ,:,:) = p(n,:,:)
        !p(n+1,:,:) = p(1,:,:)
      case(2)
        !p(:,0  ,:) = p(:,n,:)
        !p(:,n+1,:) = p(:,1,:)
      case(3)
        !$OMP WORKSHARE
        p(:,:,0  ) = p(:,:,n)
        p(:,:,n+1) = p(:,:,1)
        !$OMP END WORKSHARE
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0  ,:,:) = factor+sgn*p(1,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n+1,:,:) = factor+sgn*p(n,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0  ,:) = factor+sgn*p(:,1,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n+1,:) = factor+sgn*p(:,n,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0  ) = factor+sgn*p(:,:,1)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n+1) = factor+sgn*p(:,:,n)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(0,:,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(n  ,:,:) = factor
            p(n+1,:,:) = p(n-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,0,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,n  ,:) = factor
            p(:,n+1,:) = p(:,n-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,0) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,n  ) = factor
            p(:,:,n+1) = p(:,:,n-1)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(0,:,:) = 1.*factor + p(1  ,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(n,:,:) = 1.*factor + p(n-1,:,:)
            p(n+1,:,:) = p(n,:,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,0,:) = 1.*factor + p(:,1  ,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,n,:) = 1.*factor + p(:,n-1,:)
            p(:,n+1,:) = p(:,n,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,0) = 1.*factor + p(:,:,1  )
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,n) = 1.*factor + p(:,:,n-1)
            p(:,:,n+1) = p(:,:,n) ! not needed
            !$OMP END WORKSHARE
          endif
        end select
      endif
    end select
    return
  end subroutine set_bc
  !
  subroutine outflow(n,idir,dl,dzf,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(rp) :: dx,dy,dxi,dyi
    real(rp), dimension(0:n(3)+1) :: dzfi
    integer :: i,j,k
    !
    dx   = dl(1)     
    dxi  = dl(1)**(-1)
    dy   = dl(2)     
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    !
    ! determine face velocity from zero divergence
    !
    select case(idir)
    case(1) ! x direction, right
      if(right.eq.MPI_PROC_NULL) then
        i = n(1) + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k) &
        !$OMP SHARED(n,i,u,v,w,dx,dyi,dzfi)
        do k=1,n(3)
          do j=1,n(2)
            u(i  ,j,k) = u(i-1,j,k) - dx*((v(i,j,k)-v(i,j-1,k))*dyi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
            u(i+1,j,k) = u(i  ,j,k) ! not needed
          enddo
        enddo
        !$OMP END PARALLEL DO
      endif
    case(2) ! y direction, back
      if(back.eq.MPI_PROC_NULL) then
        j = n(2) + 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k) &
        !$OMP SHARED(n,j,u,v,w,dy,dxi,dzfi)
        do k=1,n(3)
          do i=1,n(1)
            v(i,j  ,k) = v(i,j-1,k) - dy*((u(i,j,k)-u(i-1,j,k))*dxi+(w(i,j,k)-w(i,j,k-1))*dzfi(k))
            v(i,j+1,k) = v(i,j  ,k) ! not needed
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(3) ! z direction, top
      k = n(3) + 0
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j) &
      !$OMP SHARED(n,k,u,v,w,dzf,dxi,dyi)
      do j=1,n(2)
        do i=1,n(1)
          w(i,j,  k) = w(i,j,k-1) - dzf(k)*((u(i,j,k)-u(i-1,j,k))*dxi+(v(i,j,k)-v(i,j-1,k))*dyi)
          w(i,j,k+1) = w(i,j,k  ) ! not needed
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(-1) ! x direction, left
      if(left.eq.MPI_PROC_NULL) then
        i = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(j,k) &
        !$OMP SHARED(n,i,u,v,w,dx,dyi,dzfi)
        do k=1,n(3)
          do j=1,n(2)
            u(i,j,k) = u(i+1,j,k) + dx*((v(i+1,j,k)-v(i+1,j-1,k))*dyi+(w(i+1,j,k)-w(i+1,j,k-1))*dzfi(k))
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-2) ! y direction, front
      if(front.eq.MPI_PROC_NULL) then
        j = 0
        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(i,k) &
        !$OMP SHARED(n,j,u,v,w,dy,dxi,dzfi)
        do k=1,n(3)
          do i=1,n(1)
            v(i,j,k) = v(i,j+1,k) + dy*((u(i,j+1,k)-u(i-1,j+1,k))*dxi+(w(i,j+1,k)-w(i,j+1,k-1))*dzfi(k))
          enddo
        enddo 
        !$OMP END PARALLEL DO
      endif
    case(-3) ! z direction, bottom
      k = 0
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j) &
      !$OMP SHARED(n,k,u,v,w,dzf,dxi,dyi)
      do j=1,n(2)
        do i=1,n(1)
          w(i,j,k) = w(i,j,k+1) + dzf(k)*((u(i,j,k+1)-u(i-1,j,k+1))*dxi+(v(i,j,k+1)-v(i,j-1,k+1))*dyi)
        enddo
      enddo 
      !$OMP END PARALLEL DO
    end select
    return
  end subroutine outflow
  !
  subroutine inflow(n,idir,dl,dzf,vel2d,u,v,w)
    implicit none
    integer, intent(in), dimension(3) :: n
    integer, intent(in) :: idir
    real(rp), intent(in), dimension(3) :: dl
    real(rp), intent(in), dimension(0:) :: dzf
    real(rp), dimension(0:,0:), intent(in) :: vel2d
    real(rp), dimension(0:,0:,0:), intent(inout) :: u,v,w
    real(rp) :: dx,dy,dxi,dyi
    real(rp), dimension(0:n(3)+1) :: dzfi
    integer :: i,j,k
    !
    dx   = dl(1)
    dxi  = dl(1)**(-1)
    dy   = dl(2)
    dyi  = dl(2)**(-1)
    dzfi = dzf**(-1)
    select case(idir)
      case(1) ! x direction
        if(left.eq.MPI_PROC_NULL) then
          i = 0
          do k=1,n(3)
            do j=1,n(2)
              u(i,j,k) = vel2d(j,k)
            enddo
          enddo 
        endif
      case(2) ! y direction
        j = 0
        if(front.eq.MPI_PROC_NULL) then
          do k=1,n(3)
            do i=1,n(1)
              v(i,j,k) = vel2d(i,k)
            enddo
          enddo 
        endif
      case(3) ! z direction
        k = 0
        do j=1,n(2)
          do i=1,n(1)
            w(i,j,k) = vel2d(i,j)
          enddo
        enddo 
    end select
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,n,rhsbx,rhsby,rhsbz,p)
    implicit none
    character, intent(in), dimension(3) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(:,:,0:) :: rhsbx,rhsby,rhsbz
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    integer , dimension(3) :: q
    integer :: idir
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(left.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1        ,1:n(2),1:n(3)) = p(1        ,1:n(2),1:n(3)) + rhsbx(:,:,0)
      !$OMP END WORKSHARE
    endif  
    if(right.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(n(1)-q(1),1:n(2),1:n(3)) = p(n(1)-q(1),1:n(2),1:n(3)) + rhsbx(:,:,1)
      !$OMP END WORKSHARE
    endif
    if(front.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),1        ,1:n(3)) = p(1:n(1),1        ,1:n(3)) + rhsby(:,:,0)
      !$OMP END WORKSHARE
    endif
    if(back.eq.MPI_PROC_NULL) then
      !$OMP WORKSHARE
      p(1:n(1),n(2)-q(2),1:n(3)) = p(1:n(1),n(2)-q(2),1:n(3)) + rhsby(:,:,1)
      !$OMP END WORKSHARE
    endif
    !$OMP WORKSHARE
    p(1:n(1),1:n(2),1        ) = p(1:n(1),1:n(2),1        ) + rhsbz(:,:,0)
    p(1:n(1),1:n(2),n(3)-q(3)) = p(1:n(1),1:n(2),n(3)-q(3)) + rhsbz(:,:,1)
    !$OMP END WORKSHARE
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(n,idir,p)
    implicit none
    integer , dimension(2), intent(in) :: n
    integer , intent(in) :: idir
    real(rp), dimension(0:,0:,0:), intent(inout) :: p
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(1     ,0,0),1,xhalo,left ,0, &
                        p(n(1)+1,0,0),1,xhalo,right,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(n(1),0,0),1,xhalo,right,0, &
                        p(0   ,0,0),1,xhalo,left ,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0     ,0,0),1,xhalo,left ,1, &
         !               comm_cart,requests(2),error)
         !call MPI_IRECV(p(n(1)+1,0,0),1,xhalo,right,0, &
         !               comm_cart,requests(1),error)
         !call MPI_ISSEND(p(n(1),0,0),1,xhalo,right,1, &
         !               comm_cart,requests(4),error)
         !call MPI_ISSEND(p(1   ,0,0),1,xhalo,left ,0, &
         !               comm_cart,requests(3),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(p(0,1     ,0),1,yhalo,front,0, &
                        p(0,n(2)+1,0),1,yhalo,back ,0, &
                        comm_cart,status,ierr)
      call MPI_SENDRECV(p(0,n(2),0),1,yhalo,back ,0, &
                        p(0,0   ,0),1,yhalo,front,0, &
                        comm_cart,status,ierr)
         !call MPI_IRECV(p(0,n(2)+1,0),1,yhalo,back ,0, &
         !               comm_cart,requests(1),error)
         !call MPI_IRECV(p(0,0     ,0),1,yhalo,front,1, &
         !               comm_cart,requests(2),error)
         !call MPI_ISSEND(p(0,1   ,0),1,yhalo,front,0, &
         !               comm_cart,requests(3),error)
         !call MPI_ISSEND(p(0,n(2),0),1,yhalo,back ,1, &
         !               comm_cart,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    end select
    return
  end subroutine updthalo
end module mod_bound

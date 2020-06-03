c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
c GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION AL3D(b1,b2,b3,cc)
c***
C      INCLUDE   'param.h'
      use param
      DOUBLE PRECISION m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
      DOUBLE PRECISION p,p12,q,teta,cs
      DOUBLE PRECISION UNTIER,V1,V2,V3
      DOUBLE PRECISION alpha, w, vm1, vm2, vm3, vm12 
      DOUBLE PRECISION a0, a1, a2, q0, sp, th ,CONST_TINY,CONST_PI
      
      DOUBLE PRECISION xi, invp, vma, vmb, vmc
      DOUBLE PRECISION, parameter :: ONE = 1.0d0, PB = 1.49d0, 
     &   PC2 = 0.239d0, PC1 = 0.132d0, 
     &   PC0 = (PB * (PB * PC2 + 4d0 * PC1 - 8d0) / 16d0), 
     &   PA = (PB * PB * (PB - 1d0))
      
      PARAMETER (UNTIER=1.d0/3.d0)
      INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS
      
      CONST_TINY = 1D-25
      CONST_PI = 3.14159265358979323846d0
c***  
c     (1) order coefficients: m1<m2<m3; (2) get ranges: V1<V2<v3;
c     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
c*(1)* 
      m1 = DMIN1(b1,b2)
      m3 = DMAX1(b1,b2)
      m2 = b3
      if (m2 .LT. m1) then
         tmp = m1
         m1 = m2
         m2 = tmp
      else if (m2 .GT. m3) then
         tmp = m3
         m3 = m2
         m2 = tmp
      endif
c*(2)*
      m12 = m1 + m2 
      pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
      V1  = m1*m1*m1/pr
      V2  = V1 + 0.5d0*(m2-m1)/m3
      if (m3 .LT. m12) then
         mm = m3
         V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +
     %        m2*m2*(m2-3.d0*m3))/pr
      else
         mm = m12
         V3 = 0.5d0*mm/m3
      endif
c*(3)*
      ch = DMIN1(cc,1.d0-cc)
c*(4)*      
      
      
      
! Aoki APPLIC method 
        vma = ABS(b1)   
        vmb = ABS(b2)   
        vmc = ABS(b3)    
    
        w = min(cc, 1.0D0 - cc)
        xi = (PB - vma) * (PB - vmb) * (PB - vmc) - PA
        invp = (xi + PC0) / ((PC2 * xi + PC1) * xi + PC0)
        alpha = 0.5d0 * exp(log(w + w) * invp) 
        if (cc > 0.5d0) alpha = 1d0 - alpha 
        AL3D = alpha
c***
      return
      end

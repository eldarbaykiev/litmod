 
      SUBROUTINE isoloveC_fund
      USE M_surf_waves

CC    outputs fundamental-mode phase velocities at fixed frequencies
CC    no eigenfunctions or derivatives are written
CC    option of no outfile (say none)

c     this is an adaptation of Guy Master's program tor1.f

c     differences are:
c (1) modes are computed at fixed frequency for non-integer angular
c     order
c (2) output is conform that given by surfc.f - at this stage anisotropy
c     is ignored (although it is still used in the computations, the
c     model output and partial derivatives are for vph and vsv)
c     the eigenvector scales to produce unit kinetic energy as
c     defined by Dahlen, GJ, 59:19-42,1979, and unlike Aki&Richards,
c     who omit a factor l(l+1) under the integral. Units are MKS on
c     output (8/92 G.N.).

c compile: g77 -O -fzeros -Wall -fbounds-checks -fno-automatic -finit-local-zero isoloveC_fund.f isoltableC_fund.f -oisoloveC_fund_lit 
c source text of several subroutines is ratfor, in isoltableC_fund.r
      REAL*8 rn4
      character*256 filnam

       
!      print *,'enter input model file: '
!      read(*,'(a256)') filnam
!       filnam='ME01_litmod_love_newver' 
!       rn4=6371D3
!      open(7,file=filnam,status='old',form='formatted')
!      print *,'enter output file: '
!      read(*,'(a256)') filnam
!      if (filnam(1:4) .eq. 'none') then
        nfile = -1
!      else
!        nfile = 8
!      endif
      if (nfile .gt. 0) open(nfile,file=filnam,form='formatted')
!      print *,'phase velocity output file: '
!      read(*,'(a256)') filnam
!       filnam='out_phase_love.dat'  
!!!      open(9,file='dummyeig',form='unformatted')
!!!      call model(7,nfile,9)
      call model(7,nfile,rn4)
      close(7)            
!!!      close(9)
!      open(9,file=filnam,form='formatted')
      call ltableC(nfile,9,2,rn4)
      if (nfile .gt. 0) close(8)
!      close(9)
!!!      close(19)
!!!      call system("rm dummyeig")
       
      END SUBROUTINE

      subroutine detqn_love(wdim,knt,det,ifeif)
c**** supervises the integration of the equations,it returns the value
c**** of the secular determinant as det and the count of zero crossings.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl

!!!!      real*4 dcda,dcdb,dcdrh,y1,dy1dr,dum,y2,dy2dr

      common/eifx/a(5,223),dcda(223),dcdb(223),dcdrh(223),y1(223),
     +   dy1dr(223),y2(223),dy2dr(223),dum(892)
!!!!     +   dy1dr(223),y2(223),dy2dr(223),dum(3345)

      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     *  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/stdep/ls,maxlayr
      dimension ass(6)
      w=wdim/wn
      wsq=w*w
      kount=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi 
      ass(1)=1.d0
      do 5 ind=2,6
    5   ass(ind)=0.d0
      q=0.d0
      ls=nocp1
      call startl_love(nocp1,nsl,fmu,ls,q)
      if(ls.ne.nocp1) call tps(ls,ass,ifeif)
      if(ifeif.eq.1) goto 10
      call tprop(ls,nsl,ass)
      det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      if(knsw.ne.1) return
      knt=kount-2             
      if(l.eq.1) return
      knt=knt+1
      irem=mod(knt,2)
      if(irem.eq.0.and.det.lt.0.d0) return
      if(irem.ne.0.and.det.gt.0.d0) return
      knt=knt+1
      return
   10 call etprop(ls,nsl,ass)
      det=ass(2)/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      rnrm=1.d0/(w*dsqrt(ass(3)))
      cg=(fl+0.5d0)*ass(4)/(w*ass(3))
      qinv=ass(5)/(wsq*ass(3))
      wray=dsqrt(ass(6)/ass(3))
c  compute partials (isotropy assumed, follows Takeuchi and Saito 1972)
!!!      do 11 i=1,n
!!!      y1(i)=0.
!!!      dy1dr(i)=0.
!!!      dcdrh(i)=0.
!!!      dcda(i)=0.
!!!11    dcdb(i)=0.
!!!      ekin=wsq*ass(3)
      vf=w/(fl+0.5d0)
!!!      x1=vf*vf/(ekin*cg)
!!!      do 12 i=ls,nsl
!!!      z1=a(1,i)**2
!!!      z2=a(2,i)**2
!!!      rr=r(i)**2
!!!      dcdrh(i)=0.5d0*x1*(-wsq*rho(i)*rr*z1+rr*z2/lcon(i)+
!!!     +  (fl-1.0d0)*(fl+2.0d0)*ncon(i)*z1)/rho(i)
!!!      dcda(i)=0.
!!!      dcdb(i)=x1*(rr*z2/lcon(i)+(fl-1.0d0)*(fl+2.0d0)*ncon(i)*z1)*
!!!     +  sqrt(rho(i)/lcon(i))
!!!12    continue
!!!c   transform a2 to da1/dr
!!!      do 15 i=ls,nsl
!!!   15   a(2,i)=a(1,i)/r(i)+a(2,i)/(lcon(i)*(1.d0+qshear(i)*fct))
c   normalize to unite kinetic energy (but in funny units)
!!!      do 20 i=ls,nsl
!!!        do 20 j=1,2
!!!   20     a(j,i)=a(j,i)*rnrm
c   normalize to MKS, such that omega**2 x I1 = 1.0 Nm
c   note I1 includes factor fl3 (as in Dahlen, GJ, 59:19-43, 1979 and
c   unlike Aki & Richards). Remove division by sfl3 to conform to A&R.
c   note dcdb is in units 1/rn since it is an integrand
c   the factors below have been obtained as follows:
c   7.78d-10=1/sqrt(wn*wn*rhobar*rn)
c   1.22d-16=7.78d-10/rn
c   1.94d-7=vn/(rhobar*rn)
c   1.56d-7=1/rn
!!!      do 25 i=ls,nsl
!!!      y1(i)=a(1,i)*7.78808387d-10/sfl3
!!!      dy1dr(i)=a(2,i)*1.22242723d-16/sfl3
!!!      dcdrh(i)=dcdrh(i)*1.949574d-7
!!!      dcdb(i)=dcdb(i)*1.569612d-7
!!!25    continue
      
      return
      end

      subroutine tps(i,a,ifeif)
c*** toroidal mode start soln using sph bessel fns.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension a(2)
      fu=fmu(i)*(1.d0+qshear(i)*fct)
      sqk=wsq*rho(i)/fu
      zsq=r(i)*r(i)*sqk
      call bfs(fl,zsq,eps,fp)
      a(1)=r(i)
      a(2)=fu*(fp-1.d0)       
c      if(ifeif.ne.1) return
c      qrmu=fu*qshear(i)
c      z=dsqrt(zsq)
c      fnor=r(i)*r(i)/dsqrt(sqk)
c      fi1=bfint(l,z)
c      fi2=z*(z*z+fp*(fp+1)-fl3)/2.d0
c      fi3=z*fp+fi2-fl3*fi1
c      fi4=fi3-z+fi1
c      a(3)=fi2*rho(i)*fnor/sqk
c      a(4)=fi1*ncon(i)*fnor
c      a(5)=qrmu*(fi4+(fl3-1)*fi1)*fnor
c      a(6)=(lcon(i)*(fi4+fi1)+(fl3-2)*ncon(i)*fi1)*fnor
      return
      end

      subroutine etprop(jf,jl,f)
c*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,223),dum(11,223)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(6,10),s(6),f(6)
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      qy=1.d0/y+dsqrt(dabs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2) 
      if(i.eq.jl) return
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+dsqrt(dabs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=dmax1(qx,qy)
      del=step(maxo)/q
      dxs=0.d0   
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        do 20 nf=1,6
   20     s(nf)=f(nf)
        do 40 ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +        t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          h(3,ni)=ro*f(1)*f(1)*z*z
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
          t1=(h(1,ni)/z-f(1))**2
          t2=fl3m2*f(1)*f(1)
          h(4,ni)=nn*f(1)*f(1)
          h(5,ni)=(t1+t2)*ll*qshear(iq)
          h(6,ni)=ll*t1+t2*nn
   40     call rkdot(f,s,h,6,ni)
        x=y
        if(y.ne.r(i)) goto 15
      goto 10
      end          

      subroutine tprop(jf,jl,f)
c*** propagates f from jf to jl - toroidal modes ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,223),dum(11,223)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension h(2,10),s(2),f(2)
      fl3m2=fl3-2.d0
      maxo1=maxo-1
      y=r(jf)
      qy=1.d0/y+dsqrt(dabs(rho(jf)*wsq/fmu(jf)-fl3/(y*y)))
      i=jf
   10 a(1,i)=f(1)
      a(2,i)=f(2)
      if(i.eq.jl) return                                     
      iq=i
      i=i+1
      x=y
      y=r(i)
      if(y.eq.x) goto 10
      qll=1.d0+qshear(iq)*fct
      qx=qy
      qy=1.d0/y+dsqrt(dabs(rho(i)*wsq/fmu(i)-fl3/(y*y)))
      q=dmax1(qx,qy)
      del=step(maxo)/q
      dxs=0.d0   
   15   y=x+del
        if(y.gt.r(i)) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        s(1)=f(1)
        s(2)=f(2)
        do 40 ni=1,in
          z=x+c(ni)
          t=z-r(iq)
          ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
          ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
          nn=ll
          if(ifanis.ne.0) nn=(ncon(iq)+
     +      t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          z=1.d0/z
          h(1,ni)=z*f(1)+f(2)/ll
          h(2,ni)=(nn*fl3m2*z*z-ro*wsq)*f(1)-3.d0*z*f(2)
   40     call rkdot(f,s,h,2,ni)
        if(knsw.eq.1) then
          if(s(2)*f(2).le.0.d0) then
            if(f(2).eq.0.d0) then
              tes=-s(2)*s(1)
            else
              tes=f(2)*s(1)-s(2)*f(1)
            endif
            if(tes.lt.0.d0) kount=kount+1
            if(tes.gt.0.d0) kount=kount-1
          end if
        end if
        x=y
        if(y.ne.r(i)) goto 15
      goto 10
      end          

      subroutine model(iin,iout,rn4)
      USE M_surf_waves
      implicit real*8(a-h,o-z)
      integer*4 ititle(20),n4,in4

!!!!      real*4 rn4,rh4,vp4,vs4,q4

      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/vpv(223),vph(223),vsv(223),vsh(223),eta(223),wrk(2453)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n

      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
      pi=3.14159265358979d0
!!!J      read(iin,100) (ititle(i),i=1,20)
!!!J  100 format(20a4)
!!!J      read(iin,*) ifanis,tref,ifdeck
!      ifanis=0 
       ifanis=I_anis 
      tref=50 ! Reference period is 50 s by default
      ifdeck=1 !We always use the deck of layer model (rather than polynomial, ifdeck=0)
      if(ifdeck.eq.0) go to 1000
c*** card deck model ***
!!!J      read(iin,*) n,nic,noc
       n=N_L2M   !number of layers
       nic=13 ! from defalut reference model file ME01
       noc=37 ! from defalut reference model file ME01
c added if-block to allow for isotropic models too. SL Aug 93
      if (ifanis.eq.0) then
!!!J          do 101 i=1,n
!!!J101      read(iin,*) r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i)
       r(1:n)=R_L2M(1:n)
       rho(1:n)=rho_L2M(1:n)
       vpv(1:n)=vpv_L2M(1:n)
       vsv(1:n)=vsv_L2M(1:n)
       qkappa(1:n)=Q_kappa_L2M(1:n)
       qshear(1:n)=Q_s_L2M(1:n)
      else
!           do 102 i=1,n
!102      read(iin,*) r(i),rho(i),vpv(i),vsv(i),
!     +     qkappa(i),qshear(i),vph(i),vsh(i),eta(i)
       r(1:n)=R_L2M(1:n)
       rho(1:n)=rho_L2M(1:n)
       vpv(1:n)=vpv_L2M(1:n)
       vsv(1:n)=vsv_L2M(1:n)
       vph(1:n)=vph_L2M(1:n)
       vsh(1:n)=vsh_L2M(1:n)
       qkappa(1:n)=Q_kappa_L2M(1:n)
       qshear(1:n)=Q_s_L2M(1:n)
       eta(1:n)=1D0
      endif
      go to 2000
c*** polynomial model ***
 1000 read(iin,*) nreg,nic,noc,rx
      rx=rx*tau
      n=0
      knt=0
      jj=5
      if(ifanis.ne.0) jj=8
      do 10 nn=1,nreg
      read(iin,*) nlay,r1,r2
      r1=r1*tau
      r2=r2*tau
      dr=(r2-r1)/float(nlay-1)
      do 15 i=1,nlay
      n=n+1
   15 r(n)=r1+dr*float(i-1)
      do 20 j=1,jj
      read(iin,*) (wrk(i),i=1,5)
      do 20 i=1,nlay
      ind=knt+i
      rt=r(ind)/rx
      val=wrk(1)+rt*(wrk(2)+rt*(wrk(3)+rt*(wrk(4)+rt*wrk(5))))
      if(j.eq.1) rho(ind)=val*tau
      if(j.eq.2) vpv(ind)=val*tau
      if(j.eq.3) vsv(ind)=val*tau
      if(j.eq.4) qkappa(ind)=val
      if(j.eq.5) qshear(ind)=val
      if(ifanis.eq.0) goto 20
      if(j.eq.6) vph(ind)=val*tau
      if(j.eq.7) vsh(ind)=val*tau
      if(j.eq.8) eta(ind)=val
   20 continue
   10 knt=knt+nlay
 2000 if(ifanis.ne.0) go to 3000
      do 25 i=1,n
      vph(i)=vpv(i)
      vsh(i)=vsv(i)
   25 eta(i)=1.d0
 3000 if(iout.lt.0) goto 30
c*** write out model ***
      write(iout,900) (ititle(k),k=1,20),tref
  900 format(1x,20a4,' ref per =',f6.1,' secs',///,2x,'level',
     1 4x,'radius(m)',5x,'rho(kg/m3)',2x,'vpv(m/s)',4x,'vph(m/s)',4x,
     2 'vsv(m/s)',4x,'vsh(m/s)',4x,'eta',9x,'qmu ',8x,'qkap',/)
      write(iout,905) (i,r(i),rho(i),vpv(i),vph(i),vsv(i),vsh(i),
     1 eta(i),qshear(i),qkappa(i),i=1,n)
  905 format(3x,i3,f12.1,5f12.2,f12.5,2f12.2)
  30  continue
      if(r(n).lt.6.e6.or.rho(n).lt.500.0.or.vpv(n).lt.100.0) then
       write(6,*) 'WARNING!!!! Is your input in MKS units???'
!!! write(iout,*) 'WARNING!!!! Is your input in MKS units???'
c if not MKS units --> make MKS units:  (SL Aug 93)
        do 32 i=1,n
         r(i) = r(i)*tau
         rho(i) = rho(i)*tau
         vpv(i) = vpv(i)*tau
         vsv(i) = vsv(i)*tau
         vph(i) = vph(i)*tau
         vsh(i) = vsh(i)*tau
32      continue
       write(6,*) 'corrected for WARNING as follows:'
       write(6,*) 'r,rho,vpv,vsv,vph,vsh multiplied by 1000 to'
       write(6,*) 'convert to MKS units.'
!!! write(iout,*) 'corrected for WARNING as follows:'
!!! write(iout,*) 'r,rho,vpv,vsv,vph,vsh multiplied by 1000 to'
!!! write(iout,*) 'convert to MKS units.'
      endif
   60 nsl=n
      if(vsv(nsl).gt.0.d0) go to 70
   65 nsl=nsl-1
      if(vsv(nsl).le.0.d0) go to 65
   70 nicp1=nic+1
      nocp1=noc+1
      nslp1=nsl+1
      n4=n
      in4=nsl
      rn4=r(n)
!!!      write(ioeig) rn4,in4,n4
!!!      do 75 i=1,n
!!!      rn4=r(i)
!!!      vp4=vph(i)
!!!      vs4=vsv(i)
!!!      rh4=rho(i)
!!!      q4=qshear(i)
!!!75    write(ioeig) rn4,vp4,vs4,rh4,q4
c*** normalise and spline ***
      rn=r(n)
      gn=pi*bigg*rhobar*rn
      vn2=gn*rn
      vn=dsqrt(vn2)
      wn=vn/rn
      do 45 i=1,n
      r(i)=r(i)/rn
      if(i.gt.1.and.dabs(r(i)-r(i-1)).lt.1.d-7) r(i)=r(i-1)
      if(qshear(i).gt.0.d0) qshear(i)=1.d0/qshear(i)
      if(qkappa(i).gt.0.d0) qkappa(i)=1.d0/qkappa(i)
      rho(i)=rho(i)/rhobar
      lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
      ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
   45 fmu(i)=lcon(i)
      call drspln(1,n,r,rho,qro,wrk)
      call drspln(1,n,r,lcon,lspl,wrk)
      if(ifanis.eq.0) goto 80
      call drspln(1,n,r,ncon,nspl,wrk)
80      tref=0.5d0*tref/pi
      return
      end        

      subroutine startl_love(jf,jl,v,ls,q)
c*** finds start level between jf and jl using velocityv and ang. ord. l.
c*** upon entry q is the value of the exponent at r(jf) or at the turning
c*** point(q=0) depending on previous calls to startl_love. upon exit q is the
c*** value of the exponent at the starting level ls.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),qshear(223),qkappa(223),
     + rho(223),qro(3,223),lcon(223),lspl(3,223),ncon(223),nspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension rrlog(223),p(223),v(223)
      save rrlog,p
      data ifirst/1/
      if(ifirst.ne.1) goto 5
      ifirst=0
      vertno=-dlog(eps)
      do 1 i=3,n
    1 rrlog(i)=.5d0*dlog(r(i)/r(i-1))
    5 do 10 j=jf,jl
      pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
      if(pp.le.0.d0) goto 15
   10 p(j)=dsqrt(pp)
   15 p(j)=0.d0
   20 k=j
      j=j-1
      if(j.le.jf) go to 25
      q=q+rrlog(k)*(p(j)+p(k))
      if(q.lt.vertno) go to 20
      ls=j
      return
   25 ls=jf
      return
      end



      subroutine dsplin(n,x,y,q,f)
      implicit real*8(a-h,o-z)
      dimension x(3),y(3),q(3,3),f(3,3),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      a0=0.d0
      j2=n-2
      h=x(2)-x(1)
      h2=x(3)-x(1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
      b0=(y(1)*(h-h2)+y(2)*h2-y(3)*h)/y0
      b1=b0
      do 5 i=1,j2
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h-2.d0*a0
      h3a=2.d0*h-3.*a0
      h2b=h2*b0
      q(1,i)=h2/ha
      q(2,i)=-ha/(h2a*h2)
      q(3,i)=-h*h2a/h3a
      f(1,i)=(y0-h*b0)/(h*ha)
      f(2,i)=(h2b-y0*(2.d0*h-a0))/(h*h2*h2a)
      f(3,i)=-(h2b-3.d0*y0*ha)/(h*h3a)
      a0=q(3,i)
    5 b0=f(3,i)
      i=j2+1
      h=x(i+1)-x(i)
      y0=y(i+1)-y(i)
      h2=h*h
      ha=h-a0
      h2a=h*ha
      h2b=h2*b0-y0*(2.d0*h-a0)
      q(1,i)=h2/ha
      f(1,i)=(y0-h*b0)/h2a
      ha=x(j2)-x(i+1)
      y0=-h*ha*(ha+h)
      ha=ha*ha
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
      do 10 j=1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
   10 i=k
      q(1,i)=b1
      do 15 j=1,3
   15 q(j,n)=yy(j)
      return
      end



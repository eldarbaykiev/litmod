      subroutine ltablec(iout,iophase,mtype,rn4)
      USE M_surf_waves
      implicit real*8(a-h,o-z)
      integer*4 ls0,m1,kind,null
      real*4 zero4
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,
     &     fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      common/eifx/a(5,223),dcda(223),dcdb(223),dcdrh(223),y1(223), 
     &     dy1dr(223),y2(223),dy2dr(223),dum(892)
      common/stdep/ls,maxlayr
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension xl1(0:600),xl2(0:600),d1(0:600),d2(0:600),k1(0:600),k2(
     &     0:600)
      dimension xlpred(0:600),dlpred(0:600)
!      dimension freq(50), period(50)
      dimension freq(500), period(500)
      character*2 ichar(2)
      data nmx/600/,inss/5/,ichar/' R',' L'/,zero/0./,zero4/0./,null/0/
      chz=2.0d0*pi
      stepf=1.d0
      eps = 1.0d-7
      eps1=eps
      eps2=eps
      modtot=0
      call steps(eps)
      nbran = 0
      cmin = 1.5
      cmax = 7.
      i_cont=1
!      print *,'enter number of periods:'
       IF(mtype==1) THEN
        nfreq=N_ray
        period(1:nfreq)=v_ray_dat(:,1) 
       ELSE
        nfreq=N_love
        period(1:nfreq)=v_love_dat(:,1)
       ENDIF
!!      read(*,*) nfreq
!      print *,'enter values of the periods:'
!!      read(*,*) (period(i), i = 1, nfreq)
!      period(1:nfreq)=v_ray_dat(:,1)
!      write(*,*)'N_ray',N_ray
!      write(*,*)'period',period(1:nfreq)
      do 23000 i = 1, nfreq
      freq(i) = 1/period(i)
23000 continue
      maxlayr = 0
      if(.not.(iout .gt. 0))goto 23002
      write(iout,80) cmin,cmax
23002 continue
80    format(/,'Phase velocity bounds:',2g12.4,' km/s')
      if(.not.(iout .gt. 0))goto 23004
      write(iout,100) eps,eps1
23004 continue
100   format(/,'Integration precision =',g12.4,'  root precision =', 
     &     g12.4,///,3x,'mode',5x,'order', 7x,'f(hz)',10x,'t(secs)',4x,
     &     'phase vel(km/s)',6x,'grp vel(km/s)', 8x,'q',13x,'raylquo',/)
      wmin=df*chz
      wmax=nfreq*wmin
      if(.not.(nbran.lt.0.or.nbran.gt.nmx))goto 23006
      nbran=nmx
23006 continue
      df4=df
      cmax4=cmax
!!!      rewind iophase
      radius = rn4/1000.
!      write(*,*)'radius',radius
      do 23008 i=1,nfreq
      xlpred(i)=0.0
23008 continue
      do 23010 i=1,nfreq
      dlpred(i)=0.0
23010 continue
      do 23012 ifreq=1,nfreq 
      wrad=freq(ifreq) * 2.0d0 * pi
      fhz=freq(ifreq)
      flmin=radius*wrad/cmax - 0.5d0
      if(.not.(flmin.lt.10.0d0))goto 23014
      flmin=10.0d0
23014 continue
      flmax=radius*wrad/cmin - 0.5d0
      if(.not.(flmax.lt.10.0d0))goto 23016
      flmax=10.0d0
23016 continue
      if(.not.(flmax.le.flmin))goto 23018
      goto 23012
23018 continue
      knsw=1
      maxo=inss
!      IF(mtype/=1) write(*,*)'love antes de  ldetqn'
!       write(*,*)'antes ldetqn kmax'
      call ldetqn(flmin,wrad,kmax,dmax,0,mtype)
!       write(*,*)'desp ldetqn kmax'
!       write(*,*)'flmax,wrad,kmin,dmin',flmax,wrad,kmin,dmin
!       pause
      call ldetqn(flmax,wrad,kmin,dmin,0,mtype)
!        write(*,*)'desp ldetqn kmin'
!      IF(mtype/=1) write(*,*)'love tras de  ldetqn'
      mmax=kmax
      mmin=kmin+1
!      write(*,*)'mmin,mmax antes',mmin,mmax
      if(.not.(mmax.lt.mmin))goto 23020
      goto 23012
23020 continue
      if(.not.(mmax.gt.nbran))goto 23022
      mmax=nbran
23022 continue
!       write(*,*)'mmin,mmax',mmin,mmax
      do 23024 m=mmin,mmax 
      xl1(m)=flmin
      k1(m)=kmax
      d1(m)=dmax
      xl2(m)=flmax
      k2(m)=kmin
      d2(m)=dmin
23024 continue
      do 23026 m=mmin,mmax 
23028 if(.not.(k1(m).ne.k2(m)+1))goto 23029
      fltry=0.5d0*(xl1(m)+xl2(m))
      call ldetqn(fltry,wrad,ktry,dtry,0,mtype)
      do 23030 mm=ktry+1,mmax 
      if(.not.(xl2(mm).gt.fltry))goto 23032
      xl2(mm)=fltry
      k2(mm)=ktry
      d2(mm)=dtry 
23032 continue
23030 continue
      do 23034 mm=mmin,ktry 
      if(.not.(xl1(mm).lt.fltry))goto 23036
      xl1(mm)=fltry
      k1(mm)=ktry
      d1(mm)=dtry 
23036 continue
23034 continue
!      write(*,*)'keep ongoing endless loop'
      goto 23028
23029 continue
23026 continue
      do 23038 m=mmin,mmax 
      knsw=0
      maxo=8
      if(.not.(xl2(m).lt.10.0d0))goto 23040
      goto 23038
23040 continue
      call rootf(wrad,flroot,detroot,xl1(m),xl2(m),d1(m),d2(m),eps,
     *          mtype)
      if(.not.(flroot.lt.10.0d0))goto 23042
      goto 23038
23042 continue
      IF(mtype==1) THEN
!       write(*,*)'llamado por ray'
       call detqn(wrad,kroot,droot,1)
      ELSE
!       write(*,*)'llamado por love'
       call detqn_love(wrad,kroot,droot,1)
      ENDIF  
      wdiff=(wrad-wray*wn)/wrad
      gcom=vn*cg/1000.d0
      qmod=0.0
      if(.not.(qinv .gt. 0.0))goto 23044
      qmod=1.0/qinv
23044 continue
      vphase=radius*wrad/(fl+0.5d0)
      w4=wrad
      ekin=1.0
      fl4=fl
      sqll=sfl3
      ls00=max(ls,nsl-maxlayr)
      ls0=ls00
      m1=m+1
      kind=3-mtype
      if(.not.(wdiff.gt.0.2))goto 23046
      if(.not.(iout .gt. 0))goto 23048
      write(iout,190) m,ichar(mtype),flroot,fhz,1.0/fhz,vphase,gcom,
     &     qmod,wdiff
23048 continue
190   format(i5,a2,f10.2,6g16.7,' SKIP')
      goto 23038
23046 continue
      isig=+1
      if(.not.(y1(ls).lt.0.))goto 23050
      isig=-1
      do 23052 i=ls00,n 
      y1(i)=-y1(i)
      dy1dr(i)=-dy1dr(i)
      if(.not.(kind.eq.1))goto 23054
      goto 23052
23054 continue
      y2(i)=-y2(i)
      dy2dr(i)=-dy2dr(i)
23052 continue
23050 continue
      if(.not.(iout .gt. 0))goto 23056
      write(iout,200) m,ichar(mtype),flroot,fhz,1.0/fhz,vphase,gcom,
     &     qmod,wdiff,isig
23056 continue
200   format(i5,a2,f10.2,6g16.7,i5)
      vphase=vphase*1000.0
      if(.not.(qmod.gt.0.))goto 23058
      qmod=1.0/qmod
23058 continue
!      write(iophase,'(f8.3,f10.3)') REAL(1./fhz),REAL(vphase)
!      write(*,*)'vphase,i',vphase,i_cont
!Modificacion LITMOD     
      IF(mtype==1) THEN  
       v_ray_calc(i_cont)=REAL(vphase)
      ELSE
       v_love_calc(i_cont)=REAL(vphase)
      ENDIF
      i_cont=i_cont+1 
!      write(*,*) 'sub ltable, i_cont',i_cont

23038 continue
23012 continue
      return
      end
      subroutine ldetqn(fldum,wdum,kdum,ddum,ifeif,mtype)
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,
     &     fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      data nmx/600/,inss/5/
      l=fl
      fl=fldum
      fl1=fl+1.d0
      fl2=fl+fl1
      fl3=fl*fl1
      sfl3=dsqrt(fl3)
      IF(mtype==1) THEN 
       call detqn(wdum,kdum,ddum,ifeif)
      ELSE
       call detqn_love(wdum,kdum,ddum,ifeif)
      ENDIF
      f4=fl
      d4=ddum
      return
      end
      subroutine rootf(wrad,flroot,detroot,x1,x2,d1,d2,tol,mtype)
      implicit real*8 (a-h,o-z)
      parameter (itmax=100,eps=3.e-8)
      a=x1
      b=x2
      fa=d1
      fb=d2
      if(.not.(fb*fa.gt.0.))goto 23060
      stop 'root must be bracketed for rootf.'
23060 continue
      fc=fb
      do 23062 iter=1,itmax 
      if(.not.(fb*fc.gt.0.))goto 23064
      c=a
      fc=fa
      d=b-a
      e=d
23064 continue
      if(.not.(abs(fc).lt.abs(fb)))goto 23066
      a=b
      b=c
      c=a
      fa=fb
      fb=fc
      fc=fa
23066 continue
      tol1=2.*eps*abs(b)+0.5*tol
      xm=.5*(c-b)
      if(.not.(abs(xm).le.tol1 .or. fb.eq.0.))goto 23068
      flroot=b
      return
23068 continue
      if(.not.(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)))goto 23070
      s=fb/fa
      if(.not.(a.eq.c))goto 23072
      p=2.*xm*s
      q=1.-s
      goto 23073
23072 continue
      q=fa/fc
      r=fb/fc
      p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
      q=(q-1.)*(r-1.)*(s-1.)
23073 continue
      if(.not.(p.gt.0.))goto 23074
      q=-q
23074 continue
      p=abs(p)
      if(.not.(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))))goto 23076
      e=d
      d=p/q
      goto 23077
23076 continue
      d=xm
      e=d
23077 continue
      goto 23071
23070 continue
      d=xm
      e=d
23071 continue
      a=b
      fa=fb
      if(.not.(abs(d) .gt. tol1))goto 23078
      b=b+d
      goto 23079
23078 continue
      b=b+sign(tol1,xm)
23079 continue
       call ldetqn(b,wrad,kroot,detroot,0,mtype)
      fb=detroot
23062 continue
      print *, 'rootf exceeding maximum iterations.'
      flroot=b
      return
      end
      subroutine partials(y,dw,i)
      implicit real*8 (a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl, fl1,
     &     fl2,fl3,sfl3,nord,ll,kount,knsw,ifanis,iback
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223), xa2(223)
     &,xlam(223),rho(223),qro(3,223),g(223),qg(3,223), fcon(223),fspl(3,
     &223),lcon(223),lspl(3,223),ncon(223), nspl(3,223),ccon(223),cspl(
     &3,223),acon(223),aspl(3,223)
      dimension y(4),dw(3)
      a=acon(i)
      f=fcon(i)
      c=ccon(i)
      xl=lcon(i)
      xn=ncon(i)
      d=rho(i)
      z=r(i)
      z2=z*z
      y1=y(1)
      y3=y(3)
      y1r=y1/z
      y3r=y3/z
      t1=2.0d0*y1-fl3*y3
      y2=c*y(2)+f*t1/z
      y4=xl*(y(4)+y1r-y3r)
      vsv=xl/d
      if(.not.(xl.gt.0.0d0))goto 23080
      vsv=dsqrt(vsv)
23080 continue
      vph=sqrt(a/d)
      t2=t1*t1
      t3=a-xl-xl
      dw(3)=-wsq*d*z2*(y1*y1+fl3*y3*y3)+z2*y2*y2/c+ (a-xn-f*f/c)*t2+(fl-
     &     1.0d0)*fl3*(fl+2.0d0)*xn*y3*y3
      if(.not.(xl.gt.0.0d0))goto 23082
      dw(3)= fl3*z2*y4*y4/xl+dw(3)
23082 continue
      dw(3)=dw(3)-2.0d0*d*z*y1*g(i)*t1
      dw(3)=0.5d0*w*dw(3)/d
      dw(1)=((z*y2+2.0d0*f*xl*t1/t3)**2)/c+a*t2*(1.0d0-f*f*a/(c*t3*t3))
      dw(1)=w*dw(1)/vph
      if(.not.(xl.gt.0.0d0))goto 23084
      dw(2)=fl3*z2*y4*y4/xl-4.0d0*f*xl*z*y(2)*t1/t3+2.0d0*xn*(fl3*y3*(
     &     y1- y3)-y1*t1)
      dw(2)=w*dw(2)/vsv
      goto 23085
23084 continue
      dw(2)=0.0d0
23085 continue
      return
      end

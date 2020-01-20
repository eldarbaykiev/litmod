 

       SUBROUTINE rayleighC_fund
       USE M_surf_waves
!    outputs fundamental-mode phase velocities at fixed frequencies
!CC    no eigenfunctions or derivatives are written
!CC    option of no outfile (say none)

!c     this is an adaptation of Guy Masters program sph1.f

!c     differences are:
!c (1) modes are computed at fixed frequency for non-integer angular
!c     order
!c (2) output is conform that given by surfc.f - at this stage anisotropy
!c     is ignored (although it is still used in the computations, the
!c     model output and partial derivatives are for vph and vsv)
!c     the eigenvector scales to produce unit kinetic energy as
!c     defined by Dahlen, GJ, 59:19-42,1979, and unlike Aki&Richards,
!c     who omit a factor l(l+1) under the integral. Units are MKS on
!c     output (8/92 G.N.).

!c compile:  g77 -O -fzeros -Wall -fbounds-checks -fno-automatic -finit-local-zero isoraylC_fund.f isoltableC_fund.f -oisoraylC_fund_lit
! gfortran -fno-automatic isoraylC_fund.f isoltableC_fund.f -oisoraylC_fund_lit

!c source text of several subroutines is ratfor, in ltable1.r
      REAL*8 rn4
      character*256 filnam
       write(*,*)'1'
!      print *,'enter input model file: '
!!      read(*,'(a256)') filnam
!       filnam='ME01_litmod_ray_newver'
!       rn4=6371D3
!      open(7,file=filnam,status='old',form='formatted')
!      print *,'enter output file: '
!!      read(*,'(a256)') filnam
!!       if (filnam(1:4) .eq. 'none') then
        nfile = -1
!!       else
!!         nfile = 8
!!       endif
      write(*,*)'2'
      if (nfile .gt. 0) open(nfile,file=filnam,form='formatted')
       write(*,*)'3'
!      print *,'phase velocity output file: '
!!      read(*,'(a256)') filnam
! Output file
!         filnam='out_phase_ray.dat'      
!!!      open(9,file='dummyeig',form='unformatted')
!!!      call modeln(7,nfile,9,rn4)
!      write(*,*) 'isoray antes model'
      call modeln(7,nfile,rn4)
      write(*,*)'4'
      close(7)
!!!      close(9)
!      open(9,file=filnam,form='formatted')
!       write(*,*) 'isoray desp model'
      call ltableC(nfile,9,1,rn4)
      write(*,*)'5'
!       write(*,*) 'isoray tras ltable'
      if (nfile .gt. 0) close(8)
      write(*,*)'6'
!      close(9)
!!!      close(19)
!!!      call system("rm dummyeig")
      END SUBROUTINE 

      subroutine detqn(wdim,knt,det,ifeif)
!c**** supervises the integration of the equations,it returns the value
!c**** of the secular determinant as det and the count of zero crossings.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl

!!!!      real*4 dcda,dcdb,dcdrh,y1,dy1dr,dum,y2,dy2dr

      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback

      common/eifx/a(5,223),dcda(223),dcdb(223),dcdrh(223),y1(223),
     +   dy1dr(223),y2(223),dy2dr(223),dum(892)
!!!!     +   dy1dr(223),y2(223),dy2dr(223),dum(3345)

      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/stdep/ls,maxlayr
      dimension ass(5),vf(223)!!! ,y(4),dwdp(3)
      iback=0
      w=wdim/wn
      wsq=w*w
      iexp=0
      kount=0
      fct=0.d0
      if(tref.gt.0.d0) fct=2.d0*dlog(tref*wdim)/pi
      call sdepth(wdim,ls)
      if(ls.gt.nocp1) goto 25
      if(ls.gt.nicp1) goto 20
      if(ls.le.2) then
        r10=4.5d-4*(fl+.5d0)/wdim
        if(r10.lt.r(2)) then
          r(1)=r10
          g(1)=rho(1)*r(1)*1.333333333333333d0
          ls=1
        endif
      endif
      call spsm(ls,ass)
!c*** propagate through inner core ***
      call sprpmn(ls,nic,ass,vf,5,iexp)
      r(1)=0.d0
      g(1)=0.d0
      call sfbm(ass,iback)
   20 is=max0(ls,nicp1)
      if(is.eq.ls) call fpsm(ls,ass)
!c*** propagate through outer core ***
      call fprpmn(is,noc,ass,vf,2,iexp)
      call fsbm(ass,iback)
   25 is=max0(ls,nocp1)
      if(is.eq.ls) call spsm(ls,ass)
!c*** propagate through mantle ***
      call sprpmn(is,nsl,ass,vf,5,iexp)
      if(nsl.eq.n) then
        dnorm=a(1,nsl)*a(1,nsl)
        do 30 i=2,5
   30     dnorm=dnorm+a(i,nsl)*a(i,nsl)
        det=a(5,nsl)/dsqrt(dnorm)
      else
        call sfbm(ass,iback)
!c*** propagate through ocean ***
        call fprpmn(nslp1,n,ass,vf,2,iexp)
        det=a(2,n)/dsqrt(a(1,n)*a(1,n)+a(2,n)*a(2,n))
      endif
      if(ls.gt.noc) det=-det
      if(knsw.eq.1) then
        if(ls.gt.noc) kount=kount-2
        irem=mod(kount,2)
        if(irem.eq.0.and.det.lt.0.d0) kount=kount+1
        if(irem.ne.0.and.det.gt.0.d0) kount=kount+1
        knt=kount                       
        if(knt.le.-2) then
          knt=-1
          if(det.gt.0.) det=-det
        endif
      endif
      if(ifeif.eq.0) return
!c*** this does eigenfunction calculation for spheroidal modes ***
      iback=1
      jexp=0
      do 55 i=1,4
   55   ass(i)=0.d0
      if(n.ne.nsl) then
        ass(1)=dsign(1.d0,a(1,n))
        call fprpmn(n,nslp1,ass,vf,1,jexp)
        call fsbm(ass,iback)
      else
        asi1=a(3,n)*a(3,n)
        asi2=a(4,n)*a(4,n)
        if(asi2.le.asi1) ass(1)=dsign(1.d0,a(3,n))
        if(asi2.gt.asi1) ass(2)=dsign(1.d0,a(2,n))
      end if
      nto=max0(ls,nocp1)
      call sprpmn(nsl,nto,ass,vf,4,jexp)
      if(nto.eq.ls) goto 90
      call sfbm(ass,iback)
      nto=max0(ls,nicp1)
      call fprpmn(noc,nto,ass,vf,1,jexp)
      if(nto.eq.ls) goto 90
      call fsbm(ass,iback)
      nto=max0(ls,2)
      call sprpmn(nic,nto,ass,vf,4,jexp)
   90 call eifout(ls)
!c  compute partials and copy to single precision
      cn=w/(fl+0.5d0)
!c  tt converts d(omega)/dp to d(vphase)/dp
!c  ttr scales funny units to mks for dc/da and dc/db, ttd for dc/drho.
      tt=cn*cn/(w*cg)
      ttr=tt/rn
      ttd=tt*1.949575d-7
!c  rnrm1 scales y1 to unit mks kinetic energy, rnrm2 does that for y2.
      rnrm1=7.78808387d-10
      rnrm2=rnrm1/rn
!c  component 3 and 4 have already been scaled by sfl3, so scaling is:
      rnrm3=rnrm1
      rnrm4=rnrm2
!!!      do 100 i=ls,n
!!!      do 95 j=1,4
!!!95    y(j)=a(j,i)
!!!      call partials(y,dwdp,i)
!!!      dcda(i)=dwdp(1)*ttr
!!!      dcdb(i)=dwdp(2)*ttr
!!!      dcdrh(i)=dwdp(3)*ttd
!!!      y1(i)=a(1,i)*rnrm1
!!!      dy1dr(i)=a(2,i)*rnrm2
!!!      y2(i)=a(3,i)*rnrm3
!!!      dy2dr(i)=a(4,i)*rnrm4
!!!100   continue
      return
      end

      subroutine norm(f,iexp,nvec)
      implicit real*8(a-h,o-z)
      dimension f(8)
      data econst/1048576.d0/
      size=dabs(f(1))
      do 10 j=2,nvec
   10   size=dmax1(size,dabs(f(j)))
   20  if(size.lt.1024.d0) return
        do 30 j=1,nvec
   30     f(j)=f(j)/econst
        size=size/econst
        iexp=iexp+20
        goto 20
      end

      subroutine derms(iq,z,f,fp,iknt,qff,qll,qaa)
!c*** calculates minor vector derivative (fp) in a solid ***
      implicit real*8(a-h,o-z)
      real*8 nn,ll,lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      dimension f(5),fp(5)
      if(iknt.ne.0) goto 20
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        gr=g(iq)
        ff=fcon(iq)*qff
        ll=lcon(iq)*qll
        nn=ncon(iq)*qll
        cc=ccon(iq)*qaa
        aa=acon(iq)*qaa
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
      endif
      if(ifanis.eq.0) then
        nn=ll
        cc=ff+ll+ll
        aa=cc
      else
        nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
        cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
        aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
      endif
      zr=1.d0/z
      sfl3z=sfl3*zr
      rogr=ro*gr
      c11=1.d0/cc
      c22=1.d0/ll
      dmg=aa-nn-ff*ff*c11
      zdmg=zr*dmg
      t11=-2.d0*ff*zr*c11+zr
      t12=sfl3z*ff*c11
      t21=-sfl3z
      t22=zr+zr
      s22=-ro*wsq
      s11=s22+4.d0*zr*(zdmg-rogr)
      s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
      s12=sfl3z*(rogr-zdmg-zdmg)
      s11=s11+4.d0*ro*ro
      if(iback.eq.1) goto 25
      b11=t11+t22
      b33=t11-t22
      fp(1)=b11*f(1)+c22*f(3)-c11*f(4)
      fp(2)=s12*f(1)-t21*f(3)+t12*f(4)
   20 fp(3)=s22*f(1)-2.d0*t12*f(2)+b33*f(3)+c11*f(5)
      fp(4)=-s11*f(1)+2.d0*t21*f(2)-b33*f(4)-c22*f(5)
      fp(5)=-2.d0*s12*f(2)+s11*f(3)-s22*f(4)-b11*f(5)
      return
   25 fp(1)=t22*f(1)-t21*f(2)-c22*f(3)
      fp(2)=-t12*f(1)+t11*f(2)-c11*f(4)
      fp(3)=-s22*f(1)+s12*f(2)-t22*f(3)+t12*f(4)
      fp(4)=s12*f(1)-s11*f(2)+t21*f(3)-t11*f(4)
      return
      end

      subroutine dermf(iq,z,f,fp,iknt,qff)
c*** calculates minor vector derivative (fp) in a fluid ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      dimension f(5),fp(5)
      if(iknt.ne.0) goto 10
      t=z-r(iq)
      if(t.eq.0.d0) then
        ro=rho(iq)
        flu=fcon(iq)*qff
        gr=g(iq)
      else
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        flu=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
      endif
      t21=-4.d0*ro
      zr=1.d0/z
      t12=fl3*zr*zr/wsq
      t11=gr*t12-zr
      s11=ro*(gr*gr*t12-wsq)+t21*gr*zr
      c11=-t12/ro+1.d0/flu
   10 fp(1)=t11*f(1)+c11*f(2)
      fp(2)=(s11-t21*ro)*f(1)-t11*f(2)
      return
      end

      subroutine sfbm(ass,iback)
c*** convert minor vector at a solid/fluid boundary ***
      implicit real*8(a-h,o-z)
      dimension ass(5),as(5)
      do 10 j=1,5
        as(j)=ass(j)
   10   ass(j)=0.d0
      if(iback.ne.1) then
        ass(1)=as(3)
        ass(2)=as(5)
      else
        ass(1)=-as(3)
      end if
      return
      end

      subroutine fsbm(ass,iback)
c*** convert minor vector at a fluid/solid boundary ***
      implicit real*8(a-h,o-z)
      dimension ass(5),as(5)
      do 10 j=1,5
        as(j)=ass(j)
   10   ass(j)=0.d0
      if(iback.ne.1) then
        ass(1)=as(1)
        ass(4)=-as(2)
      else
        ass(1)=-as(1)
      end if
      return
      end

      subroutine fpsm(ls,ass)
c*** spheroidal mode start solution in a fluid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      dimension ass(5)
      x=r(ls)
      fla=flam(ls)*(1.d0+xlam(ls)*fct)
      vpsq=fla/rho(ls)
      xi=g(ls)/x
      qsq=(wsq+xi-fl3*xi*xi/wsq)/vpsq
      zsq=qsq*x*x
      call bfs(fl,zsq,eps,fp)
      ass(1)=-(fl3*xi/wsq+fp)/qsq
      ass(2)=x*fla
      sum=1.d0/dsqrt(ass(1)*ass(1)+ass(2)*ass(2))
      if(ass(2).lt.0.d0) sum=-sum
      do 5 i=1,2
    5   ass(i)=ass(i)*sum
      return
      end

      subroutine spsm(ls,ass)
c*** spheroidal mode start solution in a solid region using sph. bessel fns.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      dimension a(6,2),e(15),ass(5)
      x=r(ls)
      ro=rho(ls)
      fu=fmu(ls)*(1.d0+qshear(ls)*fct)
      flu=flam(ls)*(1.d0+xlam(ls)*fct)+2.d0*fu
      vssq=fu/ro
      vpsq=flu/ro
      zeta=4.d0*ro
      xi=g(ls)/x
      alfsq=(wsq+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vpsq*vssq))
      fksq=.5d0*(alfsq+betasq+delsq)
      qsq=fksq-delsq
      zsq=qsq*x*x
      b=xi/(vssq*(betasq-qsq))
      k=1
    5   call bfs(fl,zsq,eps,fp)
        a(1,k)=fl3*b+fp
        a(2,k)=1.d0+b+b*fp
        a(3,k)=-zsq
        a(4,k)=b*a(3,k)
        a(5,k)=1.d0
        a(6,k)=fl1-fl3*b
        if(k.eq.2) goto 10
        zsq=fksq*x*x
        b=-flu/(fu*fl3*b)
        k=2
      goto 5
   10 ll=0
      do 15 i=1,3
        i1=i+1
        do 15 j=i1,4
          ll=ll+1
   15     e(ll)=a(i,1)*a(j,2)-a(j,1)*a(i,2)
      ass(1)=x*x*e(1)
      ass(2)=fu*x*sfl3*(2.d0*e(1)-e(5))
      ass(3)=fu*x*(e(3)-2.d0*e(1))
      ass(4)=x*(flu*e(4)+4.d0*fu*e(1))
      ass(5)=fu*(flu*(e(6)+2.d0*e(4))+4.d0*fu*(fl3*(e(5)-e(1))
     +     -e(3)+2.d0*e(1)))
      sum=ass(1)*ass(1)
      do 30 i=2,5
   30   sum=sum+ass(i)*ass(i)
      sum=1.d0/dsqrt(sum)
      if(ass(5).lt.0.d0) sum=-sum
      do 35 i=1,5
   35   ass(i)=ass(i)*sum
      return
      end

      subroutine sprpmn(jf,jl,f,h,nvesm,iexp)
c*** propagate a minor vector in a solid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/ar(5,223),inorm(223),kdum(4683)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(5),h(nvesm,10),s(5),fp(5),rne(4)
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvesm)
      if(iback.ne.0) then
        inorm(i)=inorm(i)+iexp
        rne(1)=-ar(1,i)*f(3)+ar(2,i)*f(2)-ar(3,i)*f(1)
        rne(2)=-ar(1,i)*f(4)+ar(4,i)*f(2)+ar(2,i)*f(1)
        rne(3)=-ar(2,i)*f(4)+ar(4,i)*f(3)-ar(5,i)*f(1)
        rne(4)=-ar(3,i)*f(4)-ar(2,i)*f(3)-ar(5,i)*f(2)
        do 10 jj=1,4
   10     ar(jj,i)=rne(jj)
      else
        inorm(i)=iexp
        do 15 j=1,nvesm
   15     ar(j,i)=f(j)
      endif
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)-fl3/(x*x)
      q=(dsqrt(dabs(fksq))+dsqrt(dabs(fksq-delsq))+2.d0/zs)/stepf
      del=jud*step(maxo)/q
      dxs=0.d0
   20   y=x+del
        if(jud*(y-r(i)).gt.0.d0) y=r(i)
        dx=y-x
        if(dx.ne.dxs) call baylis(q,maxo1)
        dxs=dx
        do 30 j=1,nvesm
   30     s(j)=f(j)
   31   do 35 ni=1,in
          z=x+c(ni)
          call derms(iq,z,f,h(1,ni),0,qff,qll,qaa)
   35     call rkdot(f,s,h,nvesm,ni)
        if(knsw.eq.1) then
          call derms(iq,y,f,fp,1,qff,qll,qaa)
          call zknt1(s,h,f,fp,x,y,1,isplit)
          if(isplit.eq.1) then
            if(dabs(dx).lt.1.d-6) goto 200
            do 100 j=1,nvesm
  100         f(j)=s(j)
            dx=dx*0.5d0
            y=x+dx
            if(jud*(y-r(i)).gt.0.d0) y=r(i)
            dx=y-x
            call baylis(q,max01)
            dxs=dx
            goto 31
          end if
        end if
  200   x=y
      if(y.ne.r(i)) goto 20
      goto 5
      end

      subroutine fprpmn(jf,jl,f,h,nvefm,iexp)
c*** propagate the minor vector in a fluid region from level jf to jl ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/ar(5,223),inorm(223),kdum(4683)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(5),h(nvefm,10),s(2),fp(2)
      data econst/1048576.d0/
      if(nvefm.eq.1) then
        do 90 i=jl,jf
          inorm(i)=inorm(i)+iexp
          do 90 j=1,2
   90       ar(j,i)=ar(j,i)*f(1)
        return
      end if
      maxo1=maxo-1
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 call norm(f,iexp,nvefm)
      inorm(i)=iexp
      do 10 j=1,nvefm
   10   ar(j,i)=f(j)
      if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min0(i,i-jud)
      qff=1.d0+xlam(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=(dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs)/stepf
      del=jud*step(maxo)/q
      dxs=0.d0
   15 y=x+del
      if(jud*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,maxo1)
      dxs=dx
      do 30 j=1,nvefm
   30   s(j)=f(j)
   31 do 35 ni=1,in
        z=x+c(ni)
        call dermf(iq,z,f,h(1,ni),0,qff)
   35   call rkdot(f,s,h,nvefm,ni)
      if(knsw.eq.1) then
        call dermf(iq,y,f,fp,1,qff)
        call zknt1(s,h,f,fp,x,y,0,isplit)
        if(isplit.eq.1) then
          if(dabs(dx).lt.1.d-6) goto 200
          do 100 j=1,nvefm
  100     f(j)=s(j)
          dx=dx*0.5d0
          y=x+dx
          if(jud*(y-r(i)).gt.0.d0) y=r(i)
          dx=y-x
          call baylis(q,max01)
          dxs=dx
          goto 31
       end if
      end if
 200  x=y
      if(y.ne.r(i)) goto 15
      goto 5
      end


      subroutine zknt1(s,sp,f,fp,x,y,ifsol,isplit)
      implicit real*8(a-h,o-z)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      dimension s(*),f(*),sp(*),fp(*)
      isplit=0
      if(ifsol.ne.0) then
        y1=s(5)
        y2=f(5)
        t1=s(3)-s(4)
        t2=f(3)-f(4)    
        y2p=s(5)+sp(5)*(y-x)
      else
        y1=s(2)
        y2=f(2)
        t1=s(1)
        t2=f(1)             
        y2p=s(2)+sp(2)*(y-x)
      endif
      if(y1*y2.le.0.d0) goto 5
      if(y2*y2p.lt.0.d0) isplit=1
      return
    5 if(t1*t2.le.0.d0) then
        isplit=1
        return
      end if
      if(y2.eq.0.d0) then
        tes=-y1*t1
      else
        tes=y2*t1-y1*t2
      endif
      if(tes.lt.0.d0) kount=kount+1
      if(tes.gt.0.d0) kount=kount-1
      return
      end

      subroutine eifout(lsmin)
c*** massages spheroidal mode eigenfunctions before output ***
      implicit real*8(a-h,o-z)
      real*8 ll,lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,223),inorm(223),kdum(4683)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension zi(4)
c*** take out normalization here
      imax=0
      do 60 iq=lsmin,n
   60   imax=max0(inorm(iq),imax)
      do 65 iq=lsmin,n
        iexp=inorm(iq)-imax
        al=0.d0
        if(iexp.ge.-80) al=2.d0**iexp
        do 65 j=1,4
   65     a(j,iq)=a(j,iq)*al
      lsm1=max0(1,lsmin-1)
      do 70 i=1,lsm1
        do 70 j=1,4
   70     a(j,i)=0.d0
      if(l.eq.1.and.lsmin.le.2) then
        a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
        a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
      end if
c*** do integrals here
      call quod(lsmin,zi)
      cg=zi(2)/(w*zi(1))
      wray=dsqrt(2.d0*zi(4)/zi(1))
      qinv=2.d0*zi(3)/(wsq*zi(1))
      rnorm=1.d0/(w*dsqrt(zi(1)))        
      rnorm1=rnorm/sfl3
c*********
      i1=min0(nic,max0(2,lsmin))
      i2=nic
    5 if(i1.ne.i2) then
        do 10 iq=i1,i2
          ff=fcon(iq)*(1.d0+xlam(iq)*fct)
          ll=lcon(iq)*(1.d0+qshear(iq)*fct)
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          d=1.d0/(ccon(iq)*(1.d0+xa2(iq)*fct))
          v=a(2,iq)
          a(2,iq)=(zr-2.d0*ff*d*zr)*a(1,iq)+sfl3z*ff*d*v+d*a(3,iq)
          a(4,iq)=-sfl3z*a(1,iq)+(zr+zr)*v+a(4,iq)/ll
   10     a(3,iq)=v
      end if
      if(i2.eq.nsl) goto 25
      i1=min0(nsl,max0(lsmin,nocp1))
      i2=nsl
      goto 5
   25 i1=min0(noc,max0(lsmin,nicp1))
      i2=noc
   30 if(i1.ne.i2) then
        do 35 iq=i1,i2
          zr=1.d0/r(iq)
          sfl3z=sfl3*zr
          ffi=1.d0/(flam(iq)*(1.d0+xlam(iq)*fct))
          p=a(2,iq)
          a(3,iq)=sfl3z*(g(iq)*a(1,iq)-p/rho(iq))/wsq
          a(2,iq)=sfl3z*a(3,iq)-a(1,iq)*zr+p*ffi
   35   a(4,iq)=sfl3z*(a(1,iq)+p*(qro(1,iq)/(rho(iq)**2)+g(iq)*ffi)/wsq)
      end if
      if(n.eq.nsl.or.i2.eq.n) goto 55
      i1=nslp1
      i2=n
      goto 30
c
55    continue
c     transform to U,dU/dr,V,dV/dr, scale to unit kinetic energy
c     (in funny units) and scale V, dV/dr by sfl3.
      do 90 iq=lsmin,n
        zr=1.d0/r(iq)
        a(1,iq)=a(1,iq)*zr
        a(2,iq)=(a(2,iq)-a(1,iq))*zr
        a(3,iq)=a(3,iq)*zr
        a(4,iq)=(a(4,iq)-a(3,iq))*zr
        a(1,iq)=a(1,iq)*rnorm
        a(2,iq)=a(2,iq)*rnorm
        a(3,iq)=a(3,iq)*rnorm1
   90   a(4,iq)=a(4,iq)*rnorm1
      if(lsmin.gt.2.or.l.gt.2) return
      if(l.eq.2) goto 95
      a(1,1)=a(1,2)-.5d0*a(2,2)*r(2)
      a(2,1)=0.d0
      a(3,1)=a(3,2)-.5d0*a(4,2)*r(2)
      a(4,1)=0.d0
      return
   95 a(2,1)=1.5d0*a(1,2)/r(2)-.5d0*a(2,2)
      a(4,1)=1.5d0*a(3,2)/r(2)-.5d0*a(4,2)
      return
      end

      subroutine quod(ls,zi)
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/a(5,223),dum(2453)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      dimension f(8),zi(4)
      do 5 i=1,8
    5   f(i)=0.
      if(ls.gt.nocp1) goto 25
      if(ls.gt.nicp1) goto 15
      if(ls.lt.2) ls=2
c*** propagate through inner core ***
      call esprop(ls,nic,f)
      do 10 i=1,4
   10   f(i+2)=f(i+4)
   15 is=max0(ls,nicp1) 
c*** propagate through outer core ***
      call efprop(is,noc,f) 
      do 20 i=1,4
   20 f(9-i)=f(7-i)
   25 is=max0(ls,nocp1)      
c*** propagate through mantle ***
      call esprop(is,nsl,f)   
      if(nsl.ne.n) then
c*** propagate through ocean ***
        do 35 i=1,4
   35     f(i+2)=f(i+4)
        call efprop(nslp1,n,f)
        do 45 i=1,4
   45     zi(i)=f(i+2)
      else                
        do 50 i=1,4
   50     zi(i)=f(i+4)
      end if
      return
      end

      subroutine efprop(jf,jl,f)
c    fprop propagates the fundamental matrix f from jf to jl (a fluid region)
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/ar(5,223),inorm(223),kdum(4683)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(6),s(6),h(6,10)
      fac=(fl+.5d0)/sfl3
      d=fl3/wsq
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)
      if(y.eq.x) goto 5
      iq=min0(i,i-jud)
      f(1)=ar(1,iq)
      f(2)=ar(2,iq)                 
      qff=1.d0+xlam(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      alfsq=(wsq+4.d0*rho(i)+xi-fl3*xi*xi/wsq)*rho(i)/flam(i)
      q=dmax1(sfl3/x,dsqrt(dabs(alfsq-fl3/(x*x)))+1.d0/zs)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(float(jud)*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx 
      do 20 j=1,6            
   20    s(j)=f(j)
      do 40 ni=1,in
        z=x+c(ni)
        t=z-zs
        zr=1.d0/z
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        gr=(g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq))))*zr
        qrka=ff*qkappa(iq)
        rogr=ro*gr
        t21=-4.d0*ro
        t12=d*zr*zr
        t11=(gr*d-1.d0)*zr
        s11=-ro*(wsq+4.d0*gr-gr*gr*d)
        c11=-t12/ro+1.d0/ff
        s11=s11-t21*ro
        h(1,ni)=t11*f(1)+c11*f(2)
        h(2,ni)=s11*f(1)-t11*f(2)
        f2=sfl3*(gr*f(1)-zr*f(2)/ro)/(wsq*z)
        f1=f(1)/z
        f1sq=f1*f1
        f2sq=f2*f2
        h(3,ni)=z*z*ro*(f1sq+f2sq)
        h(4,ni)=(sfl3*ff*f2sq+f2*((z*z*rogr-ff)*f1-ff*h(1,ni)))*fac
        h(5,ni)=(.5d0*(f1sq+fl3*f2sq)-f1*sfl3*f2+
     +    f1*h(1,ni)-f2*sfl3*h(1,ni)+.5d0*h(1,ni)*h(1,ni))*qrka
        h(6,ni)=.5d0*((4.d0*z*z*ro*(ro-gr)+ff)*f1sq+
     +    ff*(fl3*f2sq+h(1,ni)*h(1,ni)))+
     +    f1*(sfl3*(z*z*rogr-ff)*f2+ff*h(1,ni))-f2*sfl3*ff*h(1,ni)
  40    call rkdot(f,s,h,6,ni)
      x=y
      if(y.ne.r(i)) goto 15
      goto 5
      end

      subroutine esprop(jf,jl,f)
c    sprop propagates the fundamental matrix f from jf to jl (a solid region)
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl,nn,ll
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/ar(5,223),inorm(223),kdum(4683)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension f(8),s(8),h(8,10)
      data d1,d2,d3,d4,d5,d6,d7/.111111111111111d0,
     + 0.066666666666667d0,0.666666666666667d0,1.333333333333333d0,
     + 2.666666666666667d0,3.333333333333333d0,5.333333333333333d0/
      fac=(fl+.5d0)/sfl3
      jud=1
      if(jl.lt.jf) jud=-1
      y=r(jf)
      i=jf
    5 if(i.eq.jl) return
      i=i+jud
      x=y
      y=r(i)   
      if(x.eq.y) goto 5
      iq=min0(i,i-jud)
      do 10 j=1,4
   10   f(j)=ar(j,iq)
      qff=1.d0+xlam(iq)*fct
      qll=1.d0+qshear(iq)*fct
      qaa=1.d0+xa2(iq)*fct
      zs=dmin1(x,y)
      xi=g(i)/y
      vpsq=(flam(i)+2.d0*fmu(i))/rho(i)
      vssq=fmu(i)/rho(i)
      alfsq=(wsq+4.d0*rho(i)+xi)/vpsq
      betasq=wsq/vssq
      delsq=dsqrt((betasq-alfsq)**2+4.d0*fl3*xi*xi/(vssq*vpsq))
      fksq=.5d0*(alfsq+betasq+delsq)
      al=fl3/(x*x)
      aq=fksq-delsq-al
      qs=dsqrt(dabs(fksq-al))+1.d0/zs
      qf=dsqrt(dabs(aq))+1.d0/zs
      q=dmax1(sfl3/x,qs,qf)
      del=jud*step(8)/q
      dxs=0.d0
   15 y=x+del
      if(jud*(y-r(i)).gt.0.d0) y=r(i)
      dx=y-x
      if(dx.ne.dxs) call baylis(q,7)
      dxs=dx
      do 20 j=1,8
   20   s(j)=f(j)
      do 40 ni=1,in
        z=x+c(ni)
        t=z-zs
        ro=rho(iq)+t*(qro(1,iq)+t*(qro(2,iq)+t*qro(3,iq)))
        gr=g(iq)+t*(qg(1,iq)+t*(qg(2,iq)+t*qg(3,iq)))
        ff=(fcon(iq)+t*(fspl(1,iq)+t*(fspl(2,iq)+t*fspl(3,iq))))*qff
        ll=(lcon(iq)+t*(lspl(1,iq)+t*(lspl(2,iq)+t*lspl(3,iq))))*qll
        if(ifanis.eq.0)then
          nn=ll
          cc=ff+ll+ll
          aa=cc
        else
          nn=(ncon(iq)+t*(nspl(1,iq)+t*(nspl(2,iq)+t*nspl(3,iq))))*qll
          cc=(ccon(iq)+t*(cspl(1,iq)+t*(cspl(2,iq)+t*cspl(3,iq))))*qaa
          aa=(acon(iq)+t*(aspl(1,iq)+t*(aspl(2,iq)+t*aspl(3,iq))))*qaa
        endif
        qrka=d1*(4.d0*(aa+ff-nn)+cc)*qkappa(iq)
        qrmu=d2*(aa+cc-2.d0*ff+5.d0*nn+6.d0*ll)*qshear(iq)
        zr=1.d0/z
        sfl3z=sfl3*zr
        rogr=ro*gr
        c11=1.d0/cc
        c22=1.d0/ll
        dmg=aa-nn-ff*ff*c11
        zdmg=zr*dmg
        t11=-2.d0*ff*zr*c11+zr
        t12=sfl3z*ff*c11
        t21=-sfl3z
        t22=zr+zr
        s22=-ro*wsq
        s11=s22+4.d0*zr*(zdmg-rogr)+4.d0*ro*ro
        s22=s22+zr*zr*(fl3*(dmg+nn)-nn-nn)
        s12=sfl3z*(rogr-zdmg-zdmg)
        h(1,ni)=t11*f(1)+t12*f(2)+c11*f(3)
        h(2,ni)=t21*f(1)+t22*f(2)+c22*f(4)
        h(3,ni)=s11*f(1)+s12*f(2)-t11*f(3)-t21*f(4)
        h(4,ni)=s12*f(1)+s22*f(2)-t12*f(3)-t22*f(4)
        f1=f(1)*zr
        f2=f(2)*zr
        f1sq=f1*f1
        f2sq=f2*f2
        g1=h(1,ni)
        g2=h(2,ni)
        h(5,ni)=z*z*ro*(f1sq+f2sq)
        h(6,ni)=(sfl3*(ll*f1sq+aa*f2sq)+f2*((z*rogr+2.d0*(nn-aa-ll)+ff)
     +     *f1-ff*g1)+ll*g2*f1)*fac
        t2=qrka+d7*qrmu
        t3=qrka+d4*qrmu
        t4=qrka+d6*qrmu
        t5=qrka-d5*qrmu
        t6=qrka-d3*qrmu
        h(7,ni)=.5d0*((fl3*qrmu+t2)*f1sq+(2.d0*qrmu+fl3*t3)*f2sq)-
     +    f1*sfl3*t4*f2+f1*(t5*g1+sfl3*qrmu*g2)+
     +    f2*(-2.d0*qrmu*g2-sfl3*t6*g1)+.5d0*(t3*g1*g1+qrmu*g2*g2)
        h(8,ni)=.5d0*((fl3*ll+4.d0*(z*ro*(z*ro-gr)+aa-nn-ff)+cc)*f1sq+
     +    (4.d0*ll-nn-nn+fl3*aa)*f2sq+cc*g1*g1+
     +    ll*g2*g2) +f1*(sfl3*(z*ro*gr+2.d0*(nn-aa-ll)+ff)*f2+
     +    (ff+ff-cc)*g1+sfl3*ll*g2)-f2*(sfl3*ff*g1+(ll+ll)*g2)
   40   call rkdot(f,s,h,8,ni)
      x=y
      if(y.ne.r(i)) goto 15
      goto 5
      end

      subroutine grav(g,rho,qro,r,n)
c*** given rho and spline coeffs,computes gravity ***
      implicit real*8(a-h,o-z)
      dimension g(223),rho(223),qro(3,223),r(223)
      g(1)=0.d0
      do 10 i=2,n
        im1=i-1
        del=r(i)-r(im1)
        rn2=r(im1)*r(im1)
        trn=2.d0*r(im1)
        c1=rho(im1)*rn2
        c2=(qro(1,im1)*rn2+trn*rho(im1))*0.5d0
        c3=(qro(2,im1)*rn2+trn*qro(1,im1)+rho(im1))/3.d0
        c4=(qro(3,im1)*rn2+trn*qro(2,im1)+qro(1,im1))*.25d0
        c5=(trn*qro(3,im1)+qro(2,im1))*0.2d0
   10   g(i)=(g(im1)*rn2+4.d0*del*(c1+del*(c2+del*(c3+del*(c4+del*
     +    (c5+del*qro(3,im1)/6.d0))))))/(r(i)*r(i))
      return
      end

      subroutine sdepth(wdim,ls)
c*** finds starting level,ls, for a given l and w ***
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      data aw,bw,dw/-2.d-3,2.25d-3,1.28d-3/
      q=0.d0
      w=wdim/wn
      wsoc=aw+dw*fl
      if(wdim.lt.wsoc) then
        call startl(nocp1,nsl,fmu,ls,q)
        if(ls.eq.nsl) ls=ls-1
        if(ls.gt.nocp1) return
      end if
      wsic=aw+bw*fl
      if(wdim.le.wsic) then
        call startl(nicp1,noc,flam,ls,q)
        if(ls.eq.noc) ls=ls-1
        if(ls.gt.nicp1) return
      end if
      call startl(2,nic,fmu,ls,q)
      if(ls.eq.nic) ls=ls-1
      return
      end

      subroutine startl(jf,jl,v,ls,q)
c*** finds start level between jf and jl using velocityv and ang. ord. l.
c*** upon entry q is the value of the exponent at r(jf) or at the turning
c*** point(q=0) depending on previous calls to startl. upon exit q is the
c*** value of the exponent at the starting level ls.
      implicit real*8(a-h,o-z)
      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
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
    1   rrlog(i)=.5d0*dlog(r(i)/r(i-1))
    5 do 10 j=jf,jl
        pp=fl3-wsq*r(j)*r(j)*rho(j)/v(j)
        if(pp.le.0.d0) goto 15
   10   p(j)=dsqrt(pp)
   15 p(j)=0.d0
   20   k=j
        j=j-1
        if(j.le.jf) go to 25
        q=q+rrlog(k)*(p(j)+p(k))
        if(q.lt.vertno) go to 20
      ls=j
      return
   25 ls=jf
      return
      end

      subroutine modeln(iin,iout,rn4)
      USE M_surf_waves
      implicit real*8(a-h,o-z)
      integer*4 ititle(20)

!!!!      real*4 rn4,rh4,vp4,vs4,q4

      real*8 lcon,ncon,lspl,nspl
      common r(223),fmu(223),flam(223),qshear(223),qkappa(223),
     + xa2(223),xlam(223),rho(223),qro(3,223),g(223),qg(3,223),
     + fcon(223),fspl(3,223),lcon(223),lspl(3,223),ncon(223),
     + nspl(3,223),ccon(223),cspl(3,223),acon(223),aspl(3,223)
      common/bits/pi,rn,vn,wn,w,wsq,wray,qinv,cg,tref,fct,eps,fl,
     +  fl1,fl2,fl3,sfl3,nord,l,kount,knsw,ifanis,iback
      common/eifx/vpv(223),vph(223),vsv(223),vsh(223),eta(223),dum(2453)
      common/rindx/nic,noc,nsl,nicp1,nocp1,nslp1,n
      common/arem/wrk(1784)

      data bigg,tau/6.6723d-11,1.d3/,rhobar/5515.d0/
      pi=3.14159265358979d0
!!!J      read(iin,100) (ititle(i),i=1,20)
!!!J  100 format(20a4)
!!!J      read(iin,*) ifanis,tref,ifdeck
        ifanis=I_anis
!!J       ifanis=0 
       tref=50 ! Reference period is 50 s by default
       ifdeck=1 !We always use the deck of layer model (rather than polynomial, ifdeck=0)
      if(ifdeck.eq.0) go to 1000
c*** card deck model ***
!      read(iin,*) n,nic,noc
      n=N_L2M   !number of layers
      nic=13 ! from defalut reference model file ME01
      noc=37 ! from defalut reference model file ME01
c added if-block to allow for isotropic models too. SL Aug 93
      if (ifanis.eq.0) then
!!!J       do 101 i=1,n
!!!J 101      read(iin,*) r(i),rho(i),vpv(i),vsv(i),qkappa(i),qshear(i)
       r(1:n)=R_L2M(1:n)
       rho(1:n)=rho_L2M(1:n)
       vpv(1:n)=vpv_L2M(1:n)
       vsv(1:n)=vsv_L2M(1:n)
       qkappa(1:n)=Q_kappa_L2M(1:n)
       qshear(1:n)=Q_s_L2M(1:n)
      else
!!      do 102 i=1,n
!!102      read(iin,*) r(i),rho(i),vpv(i),vsv(i),
!!     +     qkappa(i),qshear(i),vph(i),vsh(i),eta(i)
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
      nsl=n
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
      acon(i)=rho(i)*vph(i)*vph(i)/vn2
      ccon(i)=rho(i)*vpv(i)*vpv(i)/vn2
      lcon(i)=rho(i)*vsv(i)*vsv(i)/vn2
      ncon(i)=rho(i)*vsh(i)*vsh(i)/vn2
      fcon(i)=eta(i)*(acon(i)-2.d0*lcon(i))
      fmu(i)=(acon(i)+ccon(i)-2.d0*fcon(i)+5.d0*ncon(i)+
     1 6.d0*lcon(i))/15.d0
      flam(i)=(4.d0*(acon(i)+fcon(i)-ncon(i))+ccon(i))/9.d0
     +    -2.d0*fmu(i)/3.d0
      rat=4.d0*fmu(i)/(3.d0*(flam(i)+2.d0*fmu(i)))
      xlam(i)=((1.d0-rat)*qkappa(i)-.5d0*rat*qshear(i))/(1.d0-1.5d0*rat)
   45 xa2(i)=(1.d0-rat)*qkappa(i)+rat*qshear(i)
      call drspln(1,n,r,rho,qro,wrk)
c*** compute g *****
      call grav(g,rho,qro,r,n)
      call drspln(1,n,r,g,qg,wrk)
      call drspln(1,n,r,fcon,fspl,wrk)
      call drspln(1,n,r,lcon,lspl,wrk)
      if(ifanis.eq.0) goto 60
      call drspln(1,n,r,acon,aspl,wrk)
      call drspln(1,n,r,ccon,cspl,wrk)
      call drspln(1,n,r,ncon,nspl,wrk)
60    tref=0.5d0*tref/pi
      return
      end

      subroutine baylis(q,maxo1)
c    baylis returns the coefficients for rks integration.
c    see e. baylis shanks(1966 a. m. s.) and references therein for the
c    coefficients. the eight runge-kutta-shanks formulae are (1-1) (2-2)
c    (3-3) (4-4) (5-5) (6-6) (7-7) (8-10). for orders greater than 4 the
c    formulae are approximate rather than exact so incurring less roundoff.
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,i
      ds=q*dabs(dx)
      do 10 j=1,maxo1
      if(ds.gt.step(j)) go to 10
      i=j
      go to 15
   10 continue
      i=maxo
   15 c(1)=0.d0
      go to (1,2,3,4,5,6,7,8),i
    1 b(1)=dx
      return
    2 c(2)=dx
      b(1)=dx
      b(2)=.5d0*dx
      b(3)=1.d0
      return
    3 c(2)=.5d0*dx
      c(3)=dx
      b(1)=c(2)
      b(2)=-dx
      b(3)=-2.d0
      b(4)=.16666666666667d0*dx
      b(5)=4.d0
      b(6)=1.d0
      return
    4 c(2)=.01d0*dx
      c(3)=.6d0*dx
      c(4)=dx
      b(1)=c(2)
      b( 2)=-.17461224489790d+02*dx
      b( 3)=-.10343618513324d+01
      b( 4)= .59691275167780d+02*dx
      b( 5)=-.10140620414448d+01
      b( 6)= .30814908546230d-01
      b( 7)=-.25555555555556d+01*dx
      b( 8)=-.11165449632656d+01
      b( 9)=-.22568165070006d+00
      b(10)=-.49077733860351d-01
      return
    5 c( 2)= 1.1111111111111d-04*dx
      c( 3)= 3.0d-01*dx
      c( 4)= 7.5d-01*dx
      c( 5)= dx
      b( 1)=c(2)
      b( 2)=-.40470000000000d+03*dx
      b( 3)=-.10007412898443d+01
      b( 4)= .25301250000000d+04*dx
      b( 5)=-.10004446420631d+01
      b( 6)= .74107010523195d-03
      b( 7)=-.11494333333333d+05*dx
      b( 8)=-.10004929965491d+01
      b( 9)= .52629261224803d-03
      b(10)=-.12029545422812d-03
      b(11)= .92592592592593d-01*dx
      b(12)= .00000000000000d+00
      b(13)= .47619047619048d+01
      b(14)= .42666666666667d+01
      b(15)= .77142857142857d+00
      return
    6 c(2)=3.3333333333333d-03*dx
      c(3)=.2d0*dx
      c(4)=.6d0*dx
      c(5)=9.3333333333333d-01*dx
      c(6)=dx
      b( 1)=c(2)
      b( 2)=-.58000000000000d+01*dx
      b( 3)=-.10344827586207d+01
      b( 4)= .64600000000000d+02*dx
      b( 5)=-.10216718266254d+01
      b( 6)= .30959752321982d-01
      b( 7)=-.62975802469136d+03*dx
      b( 8)=-.10226149961576d+01
      b( 9)= .24906685695466d-01
      b(10)=-.37737402568887d-02
      b(11)=-.54275714285714d+04*dx
      b(12)=-.10225567867765d+01
      b(13)= .25375487829097d-01
      b(14)=-.31321559234596d-02
      b(15)= .12921040478749d-03
      b(16)= .53571428571429d-01*dx
      b(17)= .00000000000000d+00
      b(18)= .61868686868687d+01
      b(19)= .77777777777778d+01
      b(20)= .40909090909091d+01
      b(21)=-.38888888888889d+00
      return
    7 c(2)=5.2083333333333d-03*dx
      c(3)=1.6666666666667d-01*dx
      c(4)=.5d0*dx
      c(5)=dx
      c(6)=8.3333333333333d-01*dx
      c(7)=dx
      b( 1)=c(2)
      b( 2)=-.25000000000000d+01*dx
      b( 3)=-.10666666666667d+01
      b( 4)= .26166666666667d+02*dx
      b( 5)=-.10421204027121d+01
      b( 6)= .61228682966918d-01
      b( 7)=-.64500000000000d+03*dx
      b( 8)=-.10450612653163d+01
      b( 9)= .51262815703925d-01
      b(10)=-.77519379844961d-02
      b(11)=-.93549382716049d+02*dx
      b(12)=-.10450293206756d+01
      b(13)= .48394546673620d-01
      b(14)=-.11877268228307d-01
      b(15)=-.39590894094358d-03
      b(16)= .35111904761905d+03*dx
      b(17)=-.10446476812124d+01
      b(18)= .52479782656724d-01
      b(19)=-.71200922221468d-02
      b(20)=-.61029361904114d-03
      b(21)= .27463212856852d-02
      b(22)= .46666666666667d-01*dx
      b(23)= .57857142857143d+01
      b(24)= .78571428571429d+01
      b(25)= .00000000000000d+00
      b(26)= b(23)
      b(27)= .10000000000000d+01
      return
    8 c(2)=.14814814814815d0*dx
      c(3)=.22222222222222d0*dx
      c(4)=.33333333333333d0*dx
      c(5)= .5d0*dx
      c(6)=.66666666666667d0*dx
      c(7)=.16666666666667d0*dx
      c(8)=dx
      c(9)=.83333333333333d0*dx
      c(10)=dx
      b( 1)=c(2)
      b( 2)= .55555555555556d-01*dx
      b( 3)= .30000000000000d+01
      b( 4)= .83333333333333d-01*dx
      b( 5)= .00000000000000d+00
      b( 6)= .30000000000000d+01
      b( 7)= .12500000000000d+00*dx
      b( 8)= .00000000000000d+00
      b( 9)= .00000000000000d+00
      b(10)= .30000000000000d+01
      b(11)= .24074074074074d+00*dx
      b(12)= .00000000000000d+00
      b(13)=-.20769230769231d+01
      b(14)= .32307692307692d+01
      b(15)= .61538461538461d+00
      b(16)= .90046296296295d-01*dx
      b(17)= .00000000000000d+00
      b(18)=-.13881748071980d+00
      b(19)= .24832904884319d+01
      b(20)=-.21182519280206d+01
      b(21)= .62467866323908d+00
      b(22)=-.11550000000000d+02*dx
      b(23)=-.35064935064935d+00
      b(24)= .50389610389610d+01
      b(25)=-.28398268398268d+01
      b(26)= .52813852813853d+00
      b(27)=-.34632034632035d+01
      b(28)=-.44097222222222d+00*dx
      b(29)=-.14173228346457d+00
      b(30)= .53385826771654d+01
      b(31)=-.35905511811023d+01
      b(32)= .70866141732284d-01
      b(33)=-.45354330708661d+01
      b(34)=-.31496062992126d-01
      b(35)= .18060975609756d+01*dx
      b(36)=-.54692775151925d-01
      b(37)= .47967589466576d+01
      b(38)=-.22795408507765d+01
      b(39)= .48615800135044d-01
      b(40)=-.34031060094530d+01
      b(41)=-.40513166779204d-01
      b(42)= .48615800135044d+00
      b(43)= .48809523809524d-01*dx
      b(44)= .65853658536585d+00
      b(45)= .66341463414634d+01
      b(46)= .52682926829268d+01
      i=10
      return
      end

      subroutine steps(eps)
c*** computes 8 dimensionless step sizes for rks integration
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      ps=dlog(eps)
      fac=1.d0
      do 2 n=1,8
      fn=n+1
      fac=fac*fn
      x=(dlog(fac)+ps)/fn
      x=dexp (x)
      s=x
      do 1 i=1,n
    1 s=x*dexp(-s/fn)
    2 step(n)=s
      return
      end

      subroutine drspln(i1,i2,x,y,q,f)
      implicit real*8(a-h,o-z)
c   rspln computes cubic spline interpolation coefficients
c   for y(x) between grid points i1 and i2 saving them in q.  the
c   interpolation is continuous with continuous first and second
c   derivitives.  it agrees exactly with y at grid points and with the
c   three point first derivitives at both end points (i1 and i2).
c   x must be monotonic but if two successive values of x are equal
c   a discontinuity is assumed and seperate interpolation is done on
c   each strictly monotonic segment.  the arrays must be dimensioned at
c   least - x(i2), y(i2), q(3,i2), and f(3,i2).  f is working storage
c   for rspln.
c                                                     -rpb
      dimension x(223),y(223),q(3,223),f(3,223),yy(3)
      equivalence (yy(1),y0)
      data yy/3*0.d0/
      j1=i1+1
      y0=0.d0
c   bail out if there are less than two points total.
      if(i2-i1)13,17,8
 8    a0=x(j1-1)
c   search for discontinuities.
      do 3 i=j1,i2
      b0=a0
      a0=x(i)
      if(a0-b0)3,4,3
 3    continue
 17   j1=j1-1
      j2=i2-2
      go to 5
 4    j1=j1-1
      j2=i-3
c   see if there are enough points to interpolate (at least three).
 5    if(j2+1-j1)9,10,11
c   only two points.  use linear interpolation.
 10   j2=j2+2
      y0=(y(j2)-y(j1))/(x(j2)-x(j1))
      do 15 j=1,3
      q(j,j1)=yy(j)
 15   q(j,j2)=yy(j)
      go to 12
c   more than two points.  do spline interpolation.
 11   a0=0.d0
      h=x(j1+1)-x(j1)
      h2=x(j1+2)-x(j1)
      y0=h*h2*(h2-h)
      h=h*h
      h2=h2*h2
c   calculate derivitive at near end.
      b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
      b1=b0
c   explicitly reduce banded matrix to an upper banded matrix.
      do 1 i=j1,j2
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
 1    b0=f(3,i)
c   take care of last two rows.
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
c   calculate derivitive at far end.
      y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
      q(3,i)=(y0*h2a+h2b)/(h*h2*(h-2.d0*a0))
      q(2,i)=f(1,i)-q(1,i)*q(3,i)
c   solve upper banded matrix by reverse iteration.
      do 2 j=j1,j2
      k=i-1
      q(1,i)=f(3,k)-q(3,k)*q(2,i)
      q(3,k)=f(2,k)-q(2,k)*q(1,i)
      q(2,k)=f(1,k)-q(1,k)*q(3,k)
 2    i=k
      q(1,i)=b1
c   fill in the last point with a linear extrapolation.
 9    j2=j2+2
      do 14 j=1,3
 14   q(j,j2)=yy(j)
c   see if this discontinuity is the last.
 12   if(j2-i2)6,13,13
c   no.  go back for more.
 6    j1=j2+2
      if(j1-i2)8,8,7
c   there is only one point left after the latest discontinuity.
 7    do 16 j=1,3
 16   q(j,i2)=yy(j)
c   fini.
 13   return
      end

      subroutine rkdot(f,s,h,nvec,ni)
c*** performs dot product with rks coefficients ***
      implicit real*8(a-h,o-z)
      common/shanks/b(46),c(10),dx,step(8),stepf,maxo,in
      dimension s(8),f(8),h(nvec,10)
      goto (1,2,3,4,5,6,7,8,9,10),ni
    1 do 21 j=1,nvec
   21 f(j)=s(j)+b(1)*h(j,1)
      return
    2 do 22 j=1,nvec
   22 f(j)=s(j)+b(2)*(h(j,1)+b(3)*h(j,2))
      return
    3 do 23 j=1,nvec
   23 f(j)=s(j)+b(4)*(h(j,1)+b(5)*h(j,2)+b(6)*h(j,3))
      return
    4 do 24 j=1,nvec
   24 f(j)=s(j)+b(7)*(h(j,1)+b(8)*h(j,2)+b(9)*h(j,3)+b(10)*h(j,4))
      return
    5 do 25 j=1,nvec
   25 f(j)=s(j)+b(11)*(h(j,1)+b(12)*h(j,2)+b(13)*h(j,3)+b(14)*h(j,4)+
     +b(15)*h(j,5))
      return
    6 do 26 j=1,nvec
   26 f(j)=s(j)+b(16)*(h(j,1)+b(17)*h(j,2)+b(18)*h(j,3)+b(19)*h(j,4)+
     +b(20)*h(j,5)+b(21)*h(j,6))
      return
    7 do 27 j=1,nvec
   27 f(j)=s(j)+b(22)*(h(j,1)+b(23)*h(j,3)+b(24)*h(j,4)+b(25)*h(j,5)+
     +b(26)*h(j,6)+b(27)*h(j,7))
      return
    8 do 28 j=1,nvec
   28 f(j)=s(j)+b(28)*(h(j,1)+b(29)*h(j,3)+b(30)*h(j,4)+b(31)*h(j,5)+
     +b(32)*h(j,6)+b(33)*h(j,7)+b(34)*h(j,8))
      return
    9 do 29 j=1,nvec
   29 f(j)=s(j)+b(35)*(h(j,1)+b(36)*h(j,3)+b(37)*h(j,4)+b(38)*h(j,5)+
     +b(39)*h(j,6)+b(40)*h(j,7)+b(41)*h(j,8)+b(42)*h(j,9))
      return
   10 do 30 j=1,nvec
   30 f(j)=s(j)+b(43)*(h(j,1)+h(j,10)+b(44)*(h(j,4)+h(j,6))+
     +b(45)*h(j,5)+b(46)*(h(j,7)+h(j,9)))
      return
      end

      subroutine entry(w,imax,kei)
      implicit real*8(a-h,o-z)
      common/mtab/we(1200),de(1200),ke(1200),wtry(600),bm(600),um(600)
      call detqn(w,kei,dei,0)
      indx=min0(max0(2*(kei-ke(1)),1),imax)
      if(indx.eq.1.and.we(1).lt.w) goto 10
      if(indx.eq.imax.and.we(imax).gt.w) goto 10
      if(kei.ne.ke(indx)) goto 5
      if(we(indx).gt.w) goto 10
      indx=indx+1
      if(we(indx).lt.w) goto 10
      return
    5 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      indx=indx+1
   10 we(indx)=w
      ke(indx)=kei
      de(indx)=dei
      return
      end

      subroutine bfs(l,xsq,eps,fp)
c  this routine calculates spherical bessel function of the ist kind.
c  fp is equivalent to (r*dj/dr)/j
c  where r is radius and j is the sbf of order l and argument x=k*r
c  the technique employs the continued fraction approach
c  described in w. lentz's article in applied qptics, vol.15, #3, 1976
      implicit real*8(a-h,o-z)
c  28-8-92: l and lp1 added to real*8 statement (GN)
      real*8 numer,nu,l,lp1
      lp1=l+1
      if(xsq.le.0.d0) goto 10
      x=dsqrt(xsq)
      rx=2.0d0/x
      nu=lp1-0.5d0
      rj=nu*rx
      rx=-rx 
      nu=nu+1
      denom=nu*rx
      numer=denom+1.0d0/rj
      rj=rj*numer/denom
    5   nu=nu+1
        rx=-rx
        a3=nu*rx
        denom=a3+1.d0/denom
        numer=a3+1.d0/numer
        ratio=numer/denom
        rj=rj*ratio
        if(dabs(dabs(ratio)-1.d0).gt.eps) goto 5
      fp=rj*x-lp1
      return
c  series solution
   10 f=1.d0
      fp=l
      a=1.d0
      b=l+lp1
      c=2.d0
      d=l+2.d0
   15   a=-a*xsq/(c*(b+c))
        f=f+a
        fp=fp+a*d
        if(dabs(a*d).lt.eps) goto 20
        c=c+2.d0
        d=d+2.d0
        goto 15
   20 fp=fp/f
      return
      end


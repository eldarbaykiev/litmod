c routines only called by vertex/meemum, could be combined with nlib.f

      subroutine lpopt0 (idead)
c-----------------------------------------------------------------------
c lpopt0 - calls lp minimization after a call to initlp. lpopt0
c does the minimization, writes error messages if necessary.

c this is an utterly stupid formulation of the lp problem because i modified
c the lp code to impose the implicit constraint that phase amounts were between
c 0 and 1, but did not impose the constraint the sum of the amounts must be
c be 1 (which would be unwise in any case). this requires that compositions and
c g's must be normalized to the total number of moles (and therefore lots of 
c extra bookkeeping). 
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw,lw,k,idead,inc

      parameter (liw=2*k1+3,lw=2*(k5+1)**2+7*k1+5*k5)  

      double precision ax(k5),x(k1),clamda(k1+k5),w(lw),oldt,oldp

      integer is(k1+k5),iw(liw)

      logical quit
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)  

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision g
      common/ cst2 /g(k1)

      integer jphct,istart
      common/ cst111 /jphct,istart

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      save ax, x, clamda, w, is, iw
c-----------------------------------------------------------------------
      if (.not.usv) then 

         inc = istct - 1

         oldt = t
         oldp = p
c                                logarithmic_p option
         if (lopt(14)) p = 1d1**p
c                                t_stop option
         if (t.lt.nopt(12)) t = nopt(12)

         call gall

         do k = 1, jphct
            c(k) = g(k+inc)/ctot(k+inc)
         end do

      end if 
c                                 idead = -1 tells lpnag to save parameters
c                                 for subsequent warm starts
      idead = -1
c                                 optimize by nag
      call lpnag (jphct,hcp,a,k5,b,c,is,x,ax,
     *            clamda,iw,liw,w,lw,idead,l6,istart)

      if (idead.gt.0) then
c                                 look for severe errors                                            
         call lpwarn (idead,'LPOPT ')
c                                 on severe error do a cold start.
c                                 necessary?
         istart = 0

      else if (hcp.eq.1.or.iopt(10).eq.0.or.isoct.eq.0.or.usv) then 
c                                 no refinement, find the answer
         call yclos0 (x,is,jphct) 
c                                 final processing, .true. indicates static
         call rebulk (.true.)

      else
c                                 find discretization points
c                                 for refinement
         call yclos1 (clamda,x,is,jphct,quit)
c                                 returns quit if nothing to refine
         if (quit) then 
c                                 final processing, .true. indicates static
            call rebulk (.true.)

         else 
c                                 reoptimize with refinement
            call reopt (idead)
c                                 final processing, .false. indicates dynamic
            if (idead.eq.0) call rebulk (.false.) 

         end if 
         
      end if 
     
      if (.not.usv) then 
         t = oldt
         p = oldp
      end if 

      end 

      subroutine reopt (idead)
c-----------------------------------------------------------------------
c reopt - given the results of an initial optimization for lpopt, reopt
c iteratively refines the solution by generating pseudocompounds in the
c neighborhood of the initial optimization.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer liw, lw, iter, iref, i, j, id, idead, ids, jstart, inc, 
     *        opt

      logical first

      parameter (liw=2*k21+3,lw=2*(k5+1)**2+7*k21+5*k5)  

      double precision  ax(k5), x(k21), clamda(k21+k5), w(lw)

      integer is(k21+k5), iw(liw)
c                                 -------------------------------------
c                                 global variables
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision xa,b,xc
      common/ cst313 /xa(k5,k1),b(k5),xc(k1)

      integer ikp
      common/ cst61 /ikp(k1)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      double precision g
      common/ cst2  /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c DEBUG
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus
c-----------------------------------------------------------------------
c                                 the pseudocompounds to be refined
c                                 are identified in jdv(1..npt)
      iter = 1
      jphct = 0
      iref = 0 
      jcoct = 1
      inc = istct - 1
      opt = npt
c                                 --------------------------------------
c                                 first iteration
      first = .true.

      do i = 1, npt

         id = jdv(i) + inc


c         if (ikp(id).eq.0) then #refine endmember change

         if (id.le.ipoint) then
c                                 the point is a true compound
            jphct = jphct + 1
c                                 jkp indicates which phase a point is associated with
            jkp(jphct) = -id
c                                 hkp indicates which refinement point
            hkp(jphct) = i

            g2(jphct) = g(id)/ctot(id)

            do j = 1, icp
               cp2(j,jphct) = cp(j,id)/ctot(id)
            end do 

         else 
c                                 the point is a pseudocompound, refine it
            call resub (i,id,ikp(id),iref,iter,first)

         end if
c                                 reset jdv in case of exit
         jdv(i) = jphct

      end do 

      if (iref.eq.0) return
 
      do 
c                                 iter is incremented before the operations,
c                                 i.e., on the nth iteration, iter is n+1
         iter = iter + 1
c                                 cold start
         jstart = 0 
c                                 set idead = 0 to prevent lpnag from
c                                 overwriting warm start parameters
         idead = 0 
c                                 do the optimization
         call lpnag (jphct,icp,cp2,k5,b,g2,is,x,ax,
     *               clamda,iw,liw,w,lw,idead,l6,jstart)
c                                 warn if severe error
         if (idead.gt.0) then

            call lpwarn (idead,'REOPT ')
            exit

         end if 
c                                 analyze solution, get refinement points
         call yclos2 (clamda,x,is,iter,opt)
c                                 save the id and compositions
c                                 of the refinement points, this
c                                 is necessary because resub rewrites
c                                 the xcoor array.
         call saver 

         if (iter.gt.iopt(10)) exit 

         jphct = 0 
         iref = 0 
         jcoct = 1

         first = .true.

c DEBUG DEBUG 
c         switch(1) = .true.
c         jcount(3) = 0
c         jcount(4) = 0
c         write (*,*) 'iter =',iter,' gmax =', gmax

c                                 generate new pseudocompounds
         do i = 1, npt

            ids = lkp(i)

            if (ids.lt.0) then 
c                                 the point is a true compound
               jphct = jphct + 1
               jkp(jphct) = ids
               hkp(jphct) = mkp(i)
               ids = -ids
               g2(jphct) = g(ids)/ctot(ids)

               do j = 1, icp
                  cp2(j,jphct) = cp(j,ids)/ctot(ids)
               end do 

            else
c                                 the point is a pseudocompound 
               call resub (mkp(i),i,ids,iref,iter,first)

            end if
c                                 reset jdv in case of exit
            jdv(i) = i 

         end do 

c        if (switch(1)) then  
c            jcount(5) = jcount(5) + jcount(3)
c            jcount(6) = jcount(6) + jcount(4)
c        end if 

      end do 

      end

      subroutine resub (jd,id,ids,iref,iter,first)
c----------------------------------------------------------------------
c subroutine to generate new pseudocompound compositions around the
c pseudocompound id of solution ids in iteration iter. ifst is the 
c pointer to the first pseudocompound of the solution ids and ilst
c points to the last. 
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      logical bad, first, keep

      double precision xxnc, ysum, res0, gmin

      integer i, j, k, l, m, ids, id, jd, iter, kcoct, iref, last 
c                                 -------------------------------------
c                                 functions
      double precision gsol1, ydinc
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive g and compositions
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 temporary subdivision limits:
      double precision xmn,xmx,xnc
      common/ cxt108 /xmn(mst,msp),xmx(mst,msp),xnc(mst,msp)
c                                 coordinates output by subdiv
      integer ntot,npairs
      common/ cst86 /ntot,npairs

      double precision simp,prism
      common/ cxt86 /simp(k13),prism(k24)

      integer ncoor,mcoor,ndim
      common/ cxt24 /ncoor(h9),mcoor(h9),ndim(mst,h9)

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)
c                                 option values
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      integer ineg
      common/ cst91 /ineg(h9,m15)

      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      save last 
c----------------------------------------------------------------------

      if (iter.eq.1) then
c                                first iteration id array points to 
c                                original compound arrays:
         call getolx (ids,id)

      else
c                                on subsequent iterations get the y's
c                                stored in the ycoor array by routine 
c                                saver, these are reindexed copies of the
c                                coordinates originally saved in zcoor
c                                below. 
         call getxy0 (ids,id)

      end if
c DEBUG should be eq 3, to use this need a different iopt value, i.e., 
c iopt(31) is the number of metastable iteration points for yclos2
      if (.not.first.and.ids.eq.last.and.iopt(31).eq.99) then

         keep = .false.
c                                check if the point lies within the 
c                                limits of the previous refinement point, 
c                                if it does, then skip the point. 
         do i = 1, istg(ids)

            do j = 1, ndim(i,ids)

               if (x(i,j).gt.xmx(i,j).or.x(i,j).lt.xmn(i,j)) then 
                  keep = .true.
                  exit
               end if

            end do
           
            if (keep) exit 
 
         end do

         if (.not.keep) return 

      end if 
c                                load the subdivision limits into
c                                temporary limit arrays:
      last = ids 

      res0 = nopt(24)/nopt(21)**iter
      
      do i = 1, istg(ids)

         do j = 1, ndim(i,ids)

            xnc(i,j) = xncg(ids,i,j)*res0
            xxnc = xnc(i,j)*reachg(ids)

            if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
               xmn(i,j) = x(i,j) - xxnc
               xmx(i,j) = x(i,j) + xxnc

            else
c                                 conformal
               xmn(i,j) = ydinc (x(i,j),-xxnc,imdg(j,i,ids),j,i,ids)
               xmx(i,j) = ydinc (x(i,j),xxnc,imdg(j,i,ids),j,i,ids)

            end if 
c                                 changed feb 6, 2012 from xmng/xmxg
c                                 to allow hardlimits. JADC
            if (xmn(i,j).lt.xmno(ids,i,j)) xmn(i,j) = xmno(ids,i,j)
            if (xmx(i,j).gt.xmxo(ids,i,j)) xmx(i,j) = xmxo(ids,i,j)

         end do 
      end do

c DEBUG DEBUG 
c      jcount(1) = 0d0
c      jcount(2) = 0d0
      gmin = 1d99 

      call subdiv (ids,.true.)

      do 10 i = 1, ntot

         jphct = jphct + 1

         if (jphct.gt.k21) call error (58,x(1,1),k21,'resub')
c                                 convert to compositional corrdinates 
c                                 required by routine gsol, y coordinates
c                                 are placed in first array of cxt7,
c                                 store a copy of x coordinates in 
c                                 1-d array zcoor
         jkp(jphct) = ids
         hkp(jphct) = jd
         jcoor(jphct) = jcoct - 1
         kcoct = jcoct + ncoor(ids)
c                                 counter for number of non 0 or 1 compositions

         if (kcoct.gt.k20) call error (59,x(1,1),k20,'resub')

         l = (i-1)*mcoor(ids)
         m = 0

         do j = 1, istg(ids)

            ysum = 0d0

            do k = 1, ndim(j,ids)

               m = m + 1

               x(j,k) = prism(l+m)
               ysum = ysum + x(j,k)
               zcoor(jcoct) = x(j,k)

               if (x(j,k).lt.xmno(ids,j,k).and.
     *             x(j,k).gt.xmxo(ids,j,k)) then 
c                                 the composition is out of range
                  jphct = jphct - 1
                  jcoct = kcoct - ncoor(ids)
                  goto 10

               end if

               jcoct = jcoct + 1

            end do

            x(j,ispg(ids,j)) = 1d0 - ysum
            zcoor(jcoct) = x(j,ispg(ids,j))
            jcoct = jcoct + 1

         end do 

         call xtoy (ids)

         if (ksmod(ids).eq.5) then
c                                 this is an el cheapo filter for redundant
c                                 compositions, a better method would be to
c                                 do the subdivision properly.
            bad = .false.

            do j = 1, ndep(ids)

               if (y(knsp(lstot(ids)+j,ids)).gt.0d0.and.
     *             y(knsp(lstot(ids)+j,ids)).le.y(ineg(ids,j))) then
c                                 reject composition 
                  jphct = jphct - 1
                  jcoct = kcoct - ncoor(ids)
                  bad = .true. 
                  exit 

               end if 

            end do 
      
            if (bad) cycle 

         end if


c                                 call gsol to get g of the solution, gsol also
c                                 computes the p compositional coordinates
         g2(jphct) = gsol1(ids)
c                                 this is to check for invalid site fractions
c                                 arising from h&p's equipartition models.
c                                 this is redundant with the z calculations 
c                                 made in gsol1, but to avoid modifying the code
c                                 the test is repeated here, the test will be 
c                                 obsolete once equipartition is phased out.
         if (lopt(5).and.lrecip(ids)) then 
c                                 check for invalid site fractions, this is only necessary
c                                 for H&P models that assume equipartition (which is not 
c                                 implemented). 
            call zchk (pa,ids,bad)

            if (bad) then
               jphct = jphct - 1
               jcoct = jcoct - ncoor(ids)
               cycle
            end if

         end if
c                                 use the coordinates to compute the composition 
c                                 of the solution
         call csol (ids,iter,bad,gmin)

         if (bad) then 
               jphct = jphct - 1
               jcoct = jcoct - ncoor(ids)
               cycle
         end if

         iref = iref + 1

10    continue

      first = .false.

      end 

      subroutine csol (id,iter,bad,gmin)
c-----------------------------------------------------------------------
c csol computes chemical composition of solution id from the macroscopic
c endmember fraction array y or p0a (cxt7), these arrays are prepared by a prior
c call to function gsol. the composition is loaded into the array cp2 at
c position jphct.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      external deltag

      integer i,j,id,iter
      logical bad 

      double precision ctot2, deltag, dg, gmin
c                                 -------------------------------------
c                                 global variables:
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision ctot
      common/ cst3  /ctot(k1)
c                                 adaptive coordinates
      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c                                 bookkeeping variables
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer jend
      common/ cxt23 /jend(h9,m4)

      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      ctot2 = 0d0

      do i = 1, icp
         cp2(i,jphct) = 0d0
      end do  

      if (lrecip(id).or.lorder(id)) then 
c                                 solutions with dependent endmembers, p0a 
c                                 contains the p's. for ksmod=8 these are a 
c                                 reformulation of the p's to eliminate the ordered 
c                                 endmembers. p0a is constructed in function gsol.
         do i = 1, lstot(id) 
            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + p0a(i) * cp(j,jend(id,2+i))
            end do 
            ctot2 = ctot2 + p0a(i)*ctot(jend(id,2+i))
         end do 

      else 
c                                 general case (y coordinates)
         do i = 1, mstot(id)

            do j = 1, icp 
               cp2(j,jphct) = cp2(j,jphct) + y(i) * cp(j,jend(id,2+i))
            end do

            ctot2 = ctot2 + y(i)*ctot(jend(id,2+i)) 

         end do 
         
      end if
c DEBUG DEBUG 

      bad = .false. 

      if (mus) then 

         dg = deltag(jphct,iter)/ctot2

         if (dg.lt.gmin.or.dg.lt.nopt(32)) then 
            jcount(4) = jcount(4) + 1 
            if (dg.lt.gmin) gmin = dg
         else 
            bad = .true.
            jcount(3) = jcount(3) + 1
         end if 

      end if 
c                                  normalize the composition and free energy
      g2(jphct) = g2(jphct)/ctot2

      do j = 1, icp 
         cp2(j,jphct) = cp2(j,jphct)/ctot2
      end do  

      end 

      subroutine sortin 
c-----------------------------------------------------------------------
c sort the first npt values of jdv
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, imin

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      do j = 1, npt-1

         imin = jdv(j)

         do i = j+1, npt

            if (jdv(i).lt.imin) then 
               imin = jdv(i)
               jdv(i) = jdv(j)
               jdv(j) = imin
            end if
 
         end do 

      end do 

      end 

      subroutine lpwarn (idead,char)
c----------------------------------------------------------------------
c write warning messages from lpnag as called by routine 'char',
c set flags ier and idead, the optimization is a total fail if
c idead set to 1.
c----------------------------------------------------------------------
      implicit none

      integer idead, iwarn91, iwarn42, iwarn90

      character*6 char     

      double precision c

      save iwarn91, iwarn42, iwarn90

      data iwarn91, iwarn42, iwarn90/0,0,0/
c----------------------------------------------------------------------
c                                             look for errors                                            
      if (idead.eq.2.or.idead.gt.4.and.iwarn91.lt.6) then 
c                                             unbounded solution, or
c                                             other programming error.
         call warn (91,c,idead,char) 
         iwarn91 = iwarn91 + 1
         if (iwarn91.eq.5) call warn (49,c,91,'LPWARN')

      else if (idead.eq.3.and.iwarn42.lt.6) then 
c                                             no feasible solution
         call warn (42,c,idead,char)
         iwarn42 = iwarn42 + 1
         if (iwarn42.eq.6) call warn (49,c,42,'LPWARN')

      else if (idead.eq.4.and.iwarn90.lt.6) then 
c                                             iteration count exceeded,
c                                             probable cause no feasible
c                                             solution.
         call warn (90,c,idead,char) 
         iwarn90 = iwarn90 + 1
         if (iwarn90.eq.5) call warn (49,c,90,'LPWARN')

      end if

      end 

      subroutine saver 
c----------------------------------------------------------------------
c subroutine to save a copy of adaptive pseudocompound x(i,j) compositions
c in the temporary array ycoor (also lcoor) used by resub to generate
c the new zcoor array for the subsequent iteration.  
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables
      integer i, j, k, kcoct, id, ids, itic
c                                 -------------------------------------
c                                 global variables:
c                                 adaptive z coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      kcoct = 0

      do i = 1, npt

         id = jdv(i)
         ids = jkp(id)
         lkp(i) = ids
c                                 cycle on a compound
         if (ids.lt.0) cycle
c                                 it's a solution:
         lcoor(i) = kcoct
         itic = 0

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
               itic = itic + 1
               if (kcoct+itic.gt.k22) then 
                  call error (60,ycoor(1),k22,'saver')
               end if 
               ycoor(lcoor(i)+itic) = zcoor(jcoor(id)+itic)
            end do 
         end do 

         kcoct = kcoct + itic

      end do 

      end 

      subroutine getxz (jd,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the zcoor array loaded in resub.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, jd, id, ids, icoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                  xcoordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c----------------------------------------------------------------------
      icoor = jcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            icoor = icoor + 1
            x(i,j) = zcoor(icoor)
            x3(jd,i,j) = x(i,j)
         end do 
      end do 

      end 

      subroutine getxy0 (ids,id)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the ycoor array loaded in saver
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 interim storage array
      integer lcoor,lkp
      double precision ycoor
      common/ cxt14 /ycoor(k22),lcoor(k19),lkp(k19)
c----------------------------------------------------------------------
      jcoor = lcoor(id)

      do i = 1, istg(ids)
         do j = 1, ispg(ids,i)
            jcoor = jcoor + 1
            x(i,j) = ycoor(jcoor)
         end do 
      end do 

      end 

      logical function solvs1 (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of solution
c ids, called only for final solution vales by avrger.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c-----------------------------------------------------------------------
      solvs1 = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle 

         if (dabs(cp3(i,id1) - cp3(i,id2))/dcp(i,ids).gt.soltol) then 
            solvs1 = .true.
            exit 
         end if 
         
      end do 

      end 

      logical function solvs2 (id1,id2)
c-----------------------------------------------------------------------
c function to test if a solvus separates two pseudocompounds of solution
c ids, intermediate solution values. 
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer id1,id2,i

      integer icomp,iphct,icp,istct
      common/ cst6 /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct
c-----------------------------------------------------------------------
      solvs2 = .false.

      do i = 1, icp

         if (dabs(cp2(i,id1) - cp2(i,id2)).gt.soltol) then
            solvs2 = .true.
            exit
         end if  

      end do 

      end 

      subroutine avrger 
c----------------------------------------------------------------------
c avrger combines discretization points into a single solution
c composition. on output

c     np  - is the number of solutions, 
c     ncpd - is the number of true compounds
c     ntot - np+ncpd

c this routine is unecessarily complicated, because it assumes
c pseudocompounds are not ordered by solution (but they are
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      logical solvs1, check
c                                 -------------------------------------
c                                 local variables
      integer idsol(k5),kdsol(k5,k5),ids,isite,xidsol,xkdsol,irep,
     *        i,j,jdsol(k5,k5),jd,k,l,nkp(k5),xjdsol(k5)

      double precision bsol(k5,k5),cpnew(k5,k5),xx,xb(k5), 
     *                 bnew(k5),xnew(k21,mst,msp)
c                                 -------------------------------------
c                                 global variables:
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ikp
      common/ cst61 /ikp(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv
c-----------------------------------------------------------------------
c                                first check if solution endmembers are
c                                among the stable compounds:
      do i = 1, ntot
         if (kkp(i).lt.0) then 
            if (ikp(-kkp(i)).ne.0) then 
c                                we have an endmember
               nkp(i) = ikp(-kkp(i))
            else
               nkp(i) = kkp(i)
            end if 
         else 
            nkp(i) = kkp(i)
         end if 
      end do 
c                                check if any solutions
      do i = 1, ntot
         if (nkp(i).gt.0) goto 10
      end do 

      if (usv) then 
         do i = 1, ntot
            do j = i+1, ntot
               if (nkp(i).eq.nkp(j)) goto 10 
            end do 
         end do
      end if 

      np = 0
      ncpd = ntot

      goto 99
c                                figure out how many solutions
c                                are present:
10    np = 0
      ncpd = 0
c                                solvus_tolerance
      soltol = nopt(8)

      do 30 i = 1, ntot
          
         if (nkp(i).lt.0) then
c                                 the pseudocompound is a true compound
            ncpd = ncpd + 1 
            idsol(ntot) = ncpd
            bsol(ntot,ncpd) = amt(i)
            kdsol(ntot,ncpd) = nkp(i)       
            jdsol(ntot,ncpd) = i   

         else 

            do j = 1, np
c                                 compare the compound to the np solutions 
c                                 identfied so far:        
               if (kdsol(j,1).eq.nkp(i)) then 
c                                 if match check for a solvus
                  if (.not.solvs1(i,jdsol(j,idsol(j)),nkp(i))) then
c                                 the pseudocompound matches a solution
c                                 found earlier.
                     idsol(j) = idsol(j) + 1
                     bsol(j,idsol(j)) = amt(i)
                     jdsol(j,idsol(j)) = i  
                     
                     goto 30 
                     
                  end if
 
            end if 

            end do
c                                 the pseudocompound is a new solution 
c                                 phase.
            np = np + 1
            idsol(np) = 1
            kdsol(np,1) = nkp(i)
            jdsol(np,1) = i 
            bsol(np,1) = amt(i)

         end if    

30    continue  
c                                 check if a solution occurs more than once
c                                 but the occurences are not sequential (this
c                                 can only occur if an endmember is immiscible 
c                                 with a general composition
      if (np.gt.2) then
 
         do i = 1, np

            check = .false.
            irep = 0

            do j = i+1, np
               if (kdsol(j,1).ne.kdsol(i,1)) then
                  check = .true.
               else 
                  irep = irep + 1
               end if 
            end do 

            if (check.and.irep.gt.0) then 

               l = i + 1

               if (kdsol(l,1).ne.kdsol(i,1)) then 
c                                 not in sequence, find the next occurence
                  do j = i+2, np 
                     if (kdsol(i,1).eq.kdsol(j,1)) exit
                  end do 
c                                 swap phase at i+1 with the one at j
                  xidsol = idsol(l)
                  xkdsol = kdsol(l,1)
                  do k = 1, xidsol
                     xb(k) = bsol(l,k)
                     xjdsol(k) = jdsol(l,k)
                  end do 

                  idsol(l) = idsol(j)
                  kdsol(l,1) = kdsol(j,1)
                  do k = 1, idsol(j)
                     bsol(l,k) = bsol(j,k)
                     jdsol(l,k) = jdsol(j,k)
                  end do 

                  idsol(j) = xidsol
                  kdsol(j,1) = xkdsol
                  do k = 1, xidsol
                     bsol(j,k) = xb(k)
                     jdsol(j,k) = xjdsol(k)
                  end do 

               end if 
            end if 
         end do 
      end if 
c                                 if a solution is represented by
c                                 more than one pseudocompound get
c                                 the everage composition
      do i = 1, np 
c                                 initialize
         bnew(i) = 0d0

         do j = 1, icomp
            cpnew(j,i) = 0d0
         end do 

         ids = kdsol(i,1)
         isite = istg(ids)

         do j = 1, isite
            do k = 1, ispg(ids,j)
               xnew(i,j,k) = 0d0
            end do 
         end do 

         do j = 1, idsol(i)
            bnew(i) = bnew(i) + amt(jdsol(i,j))
         end do 

         do j = 1, idsol(i)

            jd = jdsol(i,j)
c                                conditional in case zero mode
c                                is off:
            if (bnew(i).gt.0d0) then 

               xx =  amt(jd)/bnew(i)
c                                save the new compositions
               do k = 1, icomp
                  cpnew(k,i) = cpnew(k,i) + xx*cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = xnew(i,k,l) + xx*x3(jd,k,l)
                  end do 
               end do 
            
            else 
c                               
               do k = 1, icomp
                  cpnew(k,i) = cp3(k,jd)
               end do 

               do k = 1, isite
                  do l = 1, ispg(ids,k)
                     xnew(i,k,l) = x3(jd,k,l)
                  end do 
               end do 

            end if 

         end do 

      end do
c                                now reform the arrays kdv and b
      do i = 1, np

         amt(i) = bnew(i)
         kkp(i) = kdsol(i,1)
         ids = kkp(i)

         do j = 1, icomp
            cp3(j,i) = cpnew(j,i)
         end do

         do j = 1, istg(ids)
            do k = 1, ispg(ids,j)
c                                 set x's for sollim
               x(j,k) = xnew(i,j,k)
c                                 set x's for global storage
               x3(i,j,k) = x(j,k) 

            end do 
         end do 
c                                 check composition against solution model ranges
c                                 if auto_refine is on:
         call sollim (ids)

      end do

      if (.not.usv) then 

         do i = 1, ncpd
            k = np + i
            l = kdsol(ntot,i)
            amt(k) = bsol(ntot,i)
            kkp(k) = l
c                                for the sake of completeness load
c                                compound composition into cp3 array
            do j = 1, icomp
               cp3(j,k) = cp(j,-l)
            end do 
         end do
 
      else 

         irep = 1
         idsol(1) = kdsol(ntot,1)

         do i = 2, ncpd

            check = .false.

            do j = 1, irep

               if (kdsol(ntot,i).eq.idsol(j)) then 
                  check = .true.
                  exit 
               end if 
               
               if (check) exit

            end do 

            if (.not.check) then 
               irep = irep + 1
               idsol(irep) = kdsol(ntot,i)
            end if 
             
         end do 
         
         do i = 1, irep

            k = np + i
            l = idsol(i)

            do j = 1, ncpd
               if (kdsol(ntot,j).eq.l) amt(k) = amt(k) + bsol(ntot,j)
            end do 

            amt(k) = bsol(ntot,i)
            kkp(k) = l
c                                for the sake of completeness load
c                                compound composition into cp3 array
            do j = 1, icomp
               cp3(j,k) = cp(j,-l)
            end do 
         end do 

         ncpd = irep

      end if 

      ntot = np + ncpd

99    end 

      subroutine sollim (ids)
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'

      integer ids,i,j

      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /y(m4),z(m4),pa(m4),p0a(m4),x(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)

      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 interval limits conformal transformation
      integer intv
      double precision yint, yfrc
      common/ cst47 /yint(5,ms1,mst,h9),yfrc(4,ms1,mst,h9),intv(4)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------
c                                 set stable flag
      stable(ids) = .true.
c                                 check x-ranges
      do i = 1, istg(ids)
         do j = 1, ispg(ids,i) - 1
c                                 low limit:
            if (x(i,j).lt.xlo(j,i,ids)) then

               xlo(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).gt.xmno(ids,i,j).and.
     *            (x(i,j).le.xmng(ids,i,j).and..not.lopt(3))) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmng(ids,i,j) = xmng(ids,i,j) - nopt(10)
                     if (xmng(ids,i,j).lt.0d0) xmng(ids,i,j) = 0d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then
c                                 assymmetric stretching towards xmin
                     yint(1,j,i,ids) = yint(1,j,i,ids) - nopt(10)
                     if (yint(1,j,i,ids).lt.0d0) yint(1,j,i,ids) = 0d0
                     xmng(ids,i,j) =  yint(1,j,i,ids)

                  else 
c                                 symmetric modes, don't reset, but
c                                 set xmn to prevent future warnings
                     xmng(ids,i,j) = 0d0 
                     relax(ids) = .false.

                  end if 

                  limit(ids) = .true.

               end if 
            end if 
c                                 high limit:
            if (x(i,j).gt.xhi(j,i,ids)) then
               xhi(j,i,ids) = x(i,j)
c                                 check if solution is at an unnatural limit
               if (x(i,j).lt.xmxo(ids,i,j).and.
     *            (x(i,j).ge.xmxg(ids,i,j).and..not.lopt(3))) then
c                                 relax limits according to subdivsion model
                  if (imdg(j,i,ids).eq.0) then 
c                                 cartesian
                     xmxg(ids,i,j) = xmxg(ids,i,j) + nopt(10)
                     if (xmxg(ids,i,j).gt.1d0) xmxg(ids,i,j) = 1d0

                  else if (imdg(j,i,ids).eq.1.or.imdg(j,i,ids).eq.4)then
c                                 assymmetric stretching
                     yint(2,j,i,ids) = yint(2,j,i,ids) + nopt(10)
                     if (yint(2,j,i,ids).gt.1d0) yint(2,j,i,ids) = 1d0
                     xmxg(ids,i,j) = yint(2,j,i,ids)

                  else 
c                                 symmetric modes, don't reset
c                                 set xmx to prevent future warnings
                     xmxg(ids,i,j) = 1d0 
                     relax(ids) = .false.

                  end if 

                  limit(ids) = .true.

               end if 
            end if 
         end do 
      end do  

      end 

      subroutine outlim 
c----------------------------------------------------------------------
c subroutine to extract compositional range of endmembers in stable phases
c for auto_refine option.
c----------------------------------------------------------------------
      implicit none 

      include 'perplex_parameters.h'
c                                 -------------------------------------
c                                 local variables:
      integer i, j, k, ibad1, ibad2, ibad3, igood

      logical bad1, bad2, good, reach

      double precision num
c                                 -------------------------------------
      double precision goodc, badc
      common/ cst20 /goodc(3),badc(3)
c                                 global variables:
c                                 working arrays
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 solution limits and stability
      logical stable,limit,relax
      double precision xlo,xhi
      common/ cxt11 /xlo(m4,mst,h9),xhi(m4,mst,h9),stable(h9),limit(h9),
     *               relax(h9)
c                                 solution model counter
      integer isoct
      common/ cst79 /isoct
c                                 solution model names
      character fname*10, aname*6, lname*22
      common/ csta7 /fname(h9),aname(h9),lname(h9)
c                                 endmember pointers
      integer jend
      common/ cxt23 /jend(h9,m4)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      logical refine
      common/ cxt26 /refine

      integer isec,icopt,ifull,imsg,io3p
      common/ cst103 /isec,icopt,ifull,imsg,io3p
c DEBUG 
      integer jcount
      logical switch
      common/ debug1 /jcount(10),switch(10)

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
      
      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      double precision xmng, xmxg, xncg, xmno, xmxo, reachg
      common/ cxt6r /xmng(h9,mst,msp),xmxg(h9,mst,msp),xncg(h9,mst,msp),
     *               xmno(h9,mst,msp),xmxo(h9,mst,msp),reachg(h9)
c----------------------------------------------------------------------
      ibad1 = 0 
      ibad2 = 0 
      ibad3 = 0 
      igood = 0 
      rewind (n10)
      if (lopt(11)) rewind (n11)

      if (isoct.eq.0) goto 99

      bad1 = .false.
      bad2 = .false.
      good = .false.

      do i = 1, isoct

         if (.not.stable(i)) then
            bad1 = .true.
            ibad1 = ibad1 + 1
         else
            good = .true.
            igood = igood + 1
         end if
 
         if (limit(i)) then 
            bad2 = .true.
            if (relax(i)) then 
               ibad2 = ibad2 + 1
            else 
               ibad3 = ibad3 + 1
            end if 
         end if 

      end do 

      if (.not.refine) write (n10,*) ibad1,0,igood
c                                 write solutions present that are 
c                                 not stable
      if (bad1) then 

         write (*,1000)
         if (lopt(11)) write (n11,1000)
     
         do i = 1, isoct
            if (.not.stable(i)) then 
               write (*,'(5x,a)') fname(i) 
               if (.not.refine) write (n10,'(a)') fname(i)
               if (lopt(11)) write (n11,'(5x,a)') fname(i) 
            end if 
         end do
      end if 

      if (.not.good) goto 99
c                                 write solutions that are on an internal
c                                 limit
      if (bad2.and.icopt.gt.3) then 
c                                 adaptive minimization
         if (ibad2.gt.0) then 
c                                 solutions whose limits could be relaxed
            write (*,1080) 
            if (lopt(11)) write (n11,1080) 
            do i = 1, isoct
               if (limit(i).and.relax(i)) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

         if (ibad3.gt.0) then 
c                                 solutions whose limits could NOT be relaxed
            write (*,1090) 
            if (lopt(11)) write (n11,1090) 
            do i = 1, isoct
               if (limit(i).and.(.not.relax(i))) then
                  write (*,'(5x,a)') fname(i) 
                  if (lopt(11)) write (n11,'(5x,a)') fname(i)
               end if  
            end do
         end if

      else if (bad2) then 
c                                 non-adaptive minimization,
c                                 solutions on internal limits 
         write (*,1010) 
         if (lopt(11)) write (n11,1010) 
         do i = 1, isoct
            if (limit(i)) then
               write (*,'(5x,a)') fname(i) 
               if (lopt(11)) write (n11,'(5x,a)') fname(i)
            end if  
         end do

      end if 

      reach = .false.

      do i = 1, isoct

         if (int(reachg(i)*2d0/nopt(21)-1d0).gt.0) reach = .true.

         if (.not.stable(i)) cycle

         if (.not.refine) then

            write (n10,'(a)') fname(i)

            do j = 1, istg(i)
               do k = 1, ispg(i,j)-1
                  write (n10,*) xlo(k,j,i),xhi(k,j,i)
               end do
            end do 

         end if 

         if (istg(i).eq.1) then 
c                                 single site solution
            write (*,1020) fname(i)
            if (lopt(11)) write (n11,1020) fname(i)

            do j = 1, ispg(i,1) - 1
            
               if (ksmod(i).eq.5) then
               
                  write (*,1070) j,xlo(j,1,i),xhi(j,1,i)
                                  
               else
                
                  write (*,1030) names(jend(i,2+j)),
     *                           xlo(j,1,i),xhi(j,1,i)
               end if
                
               if (lopt(11)) write (n11,1030) 
     *                       names(jend(i,2+j)),xlo(j,1,i),xhi(j,1,i)
     
            end do 

         else
c                                 reciprocal solution
            write (*,1040) fname(i)
            if (lopt(11)) write (n11,1040) fname(i)

            do j = 1, istg(i)
               
               write (*,1050) j
               if (lopt(11)) write (n11,1050) j

               if (ispg(i,j).eq.1) then 

                  write (*,1060)
                  if (lopt(11)) write (n11,1060) 

               else

                  do k = 1, ispg(i,j) - 1
                     write (*,1070) k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                     if (lopt(11)) write (n11,1070) 
     *                              k,xlo(k,j,i),xhi(k,j,i)
c    *                     ,names(jend(i,2+indx(i,j,k)))
                  end do 

               end if 

            end do

         end if 

      end do 

      if (reach) then 

         write (*,1100)

         do i = 1, isoct
            if (int(reachg(i)*2d0/nopt(21)-1d0).eq.0) cycle
            write (*,1110) fname(i), int(reachg(i)*2d0/nopt(21)-1d0)
         end do 

      end if 

99    if (goodc(1)+badc(1).gt.0d0) then
         num = badc(1)/(badc(1)+goodc(1))*1d2
         write (*,1120) num, badc(1) + goodc(1)
         if (num.gt.1d-1) call warn (53,num,i,'OUTLIM')
         write (*,1140) dfloat(jcount(1))/(badc(1)+goodc(1))
         goodc(1) = 0d0
         badc(1) = 0d0 
      end if 

      if (lopt(32).and.jcount(4)+jcount(3).gt.0) then 
         write (*,1130) 1d2*jcount(3)/(jcount(4)+jcount(3))
      else
         write (*,'(/,a,/)') 'Metastability filtering is off.'
      end if 

      jcount(3) = 0
      jcount(4) = 0
      jcount(1) = 0

      close (n10)
      if (lopt(11)) close (n11)

1000  format (/,'WARNING: The following solutions were input, but are',
     *          ' not stable:',/)
1010  format (/,'WARNING: The following solutions have compositions on',
     *          ' an internal limit (i.e., 0<x<1)',/,'(see ranges ',
     *          'below to determine which limits should be relaxed or',
     *        /,'if executing in auto_refine mode inrease auto_refine',
     *          '_slop in perplex_option.dat):',/)
1020  format (/,'Endmember compositional ranges for model: ',a,//,5x,
     *        'Endmember   Minimum   Maximum')
1030  format (5x,a8,4x,f7.5,3x,f7.5)
1040  format (/,'Site fraction ranges for multisite model: ',a)
1050  format (/,'  Site ',i1,/,5x,'Species   Minimum   Maximum   ')
c     *          'Endmember with this species')
1060  format (8x,'Dummy site generated by model reformulation',/)
1070  format (8x,i1,6x,f7.5,3x,f7.5,3x,12(a8,1x))
1080  format (/,'WARNING: The compositions of the following solutions ',
     *        'reached internal limits',/,'that were automatically ',
     *        'relaxed.',/)
1090  format (/,'WARNING: The compositions of the following solutions ',
     *        'reached internal limits',/,'that could not be ',
     *        'automatically relaxed. To avoid this problem change',/,
     *        'the subdivision mode for ',
     *        'the solutions or increase auto_refine_slop.',/)
1100  format (/,'The following solution models have non-zero reach_',
     *          'increment:',//,t30,'reach_increment')
1110  format (4x,a,t35,i2)
1120  format (/,'The failure rate during speciation (order-disorder) ',
     *        'calculations is ',f7.3,'%',/,'out of a total of ',f12.0,
     *        ' calculations.',/)
1130  format (/,'Metastability filtering eliminated ',f7.3,'% of the ',
     *        'trial compositions prior',/,'to optimization.',/)
1140  format (/,'Average number of iterations per speciation ',
     *          'calculation:',f5.1,/)
      end 

      subroutine sorter (kdbulk,ico,jco,output)
c----------------------------------------------------------------------
c sorter compares assemblages to those already defined and reorders 
c the phases if the assemblage has been identified earlier
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,l,m,kdbulk,ico,jco,ids,ioct,inct

      logical output, reord, match, nomtch, ok 

      double precision cpt(k5,k5),xt(k5,mst,msp),bt(k5)
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk

      integer idasls,iavar,iasct,ias
      common/ cst75  /idasls(k5,k3),iavar(3,k3),iasct,ias
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
c                                 look for a match with known assemblages
      match = .false.
 
      do i = 1, iasct

         if (np.ne.iavar(1,i).or.ncpd.ne.iavar(2,i)) cycle 

         nomtch = .false.

         do j = 1, ntot

            ok = .false.

            do k = 1, ntot 

               if (idasls(k,i).eq.kkp(j)) then 

                  ok = .true.
c                                 check that the phase occurs the same 
c                                 number of times in each assemblage:
                  inct = 0 
                  ioct = 0 

                  do l = 1, np
                     if (kkp(l).eq.kkp(j)) inct = inct + 1
                     if (idasls(l,i).eq.kkp(j)) ioct = ioct + 1
                  end do 

                  if (ioct.ne.inct) then
                     nomtch = .true.  
                     exit 
                  end if 

               end if 

            end do

            if (.not.ok) nomtch = .true. 
            if (nomtch) exit 

         end do 

         if (nomtch) cycle

         match = .true.
c                                 check if reordering is necessary
         reord = .false.

         do j = 1, ntot
            if (kkp(j).eq.idasls(j,i)) cycle
            reord = .true.
            exit
         end do

         if (reord) then 
c                                 reorder the result arrays of the
c                                 current occurence to match initial 
c                                 occurence:
            do j = 1, ntot

               do k = 1, ntot

                  if (kkp(k).eq.idasls(j,i)) then
c                                 load temporary array
                  bt(j) = amt(k)

                     if (kkp(k).gt.0) then 

                        do l = 1, icomp
                           cpt(l,j) = cp3(l,k)
                        end do

                        do l = 1, istg(kkp(k))
                           do m = 1, ispg(kkp(k),l)
                              xt(j,l,m) = x3(k,l,m) 
                           end do 
                        end do 
                     end if 
c                                 this eliminates immiscible phases
                     kkp(k) = 0

                     exit 
 
                  end if 

               end do 

            end do
c                                 reload final arrays from temporary
            do j = 1, ntot

               amt(j) = bt(j)
               ids = idasls(j,i)
               kkp(j) = ids

               if (ids.gt.0) then 

                  do k = 1, icomp
                     cp3(k,j) = cpt(k,j)
                  end do

                  do k = 1, istg(ids)
                     do l = 1, ispg(ids,k)
                        x3(j,k,l) = xt(j,k,l) 
                     end do 
                  end do 
               end if 
            end do 

         end if 

         if (ibulk.gt.k2) call error (183,0d0,k2,'SORTER')
         ibulk = ibulk + 1
         iap(ibulk) = i
         kdbulk = ibulk

         exit  

      end do 

      if (.not.match) then 
c                                 the assemblage is new:
         iasct = iasct + 1
         if (iasct.gt.k3) call error (184,0d0,k3,'SORTER')

         do i = 1, ntot
            idasls(i,iasct) = kkp(i)
         end do

         ibulk = ibulk + 1
         if (ibulk.gt.k2) call error (183,0d0,k2,'BLKMAT')
         kdbulk = ibulk 
         iap(ibulk) = iasct 

         iavar(1,iasct) = np
         iavar(2,iasct) = ncpd
         iavar(3,iasct) = np + ncpd

      end if 
                                
      if (output) call outbl1 (ico,jco)
     
      end 

      subroutine outbl1 (ico,jco)
c----------------------------------------------------------------------
c output data for compositions and phases of assemblage ibulk
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ico,jco,i,j,k,ids
c                                 -------------------------------------
c                                 global variables
      integer iap,ibulk
      common/ cst74  /iap(k2),ibulk
c                                 composition and model flags
c                                 for final adaptive solution
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)
c                                 x coordinate description
      integer istg, ispg, imlt, imdg
      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 i/o
      integer io3,io4,io9
      common / cst41 /io3,io4,io9

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c----------------------------------------------------------------------
      if (io4.eq.1) return
c                                graphics output  
      write (n5,'(3(i8,1x))') ico,jco,iap(ibulk)
c                                phase molar amounts
      write (n5,1010) (amt(i),i=1,np+ncpd)
c                                solution phase compositions
      do i = 1, np
         ids = kkp(i)
         write (n5,1010) ((x3(i,j,k),k=1,ispg(ids,j)),j=1,istg(ids))
      end do 
c                                dependent potentials
      if (jpot.ne.1) write (n5,1010) (mu(i),i=1,jbulk)

1010  format (20(g16.8,1x))

      end 

      subroutine yclos1 (clamda,x,is,jphct,quit)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement. this routine is only called as preparation
c for iterative refinement.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer jphct, i, j, k, is(*), idsol(k5), kdv(h9), nsol, 
     *        mpt, iam, id, inc, jdsol(k5,k5), ldv1, ldv2,
     *        kdsol(k5), max, mcpd

      external ffirst

      logical solvus, quit

      double precision clamda(*), x(*),  slam(h9), clam

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision wmach(9)
      common /ax02za/wmach

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer isoct
      common/ cst79 /isoct

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer ikp
      common/ cst61 /ikp(k1)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus
c----------------------------------------------------------------------
      npt = 0 
      nsol = 0
      inc = istct - 1
      quit = .true.
c                                 solvus_tolerance_II, 0.25
      soltol = nopt(25)

      do i = 1, jphct

         if (is(i).ne.1) then 
c                                 make a list of found phases:
            id = i + inc
c                                 currently endmember compositions are not 
c                                 refined (this is probably a mistake, but 
c                                 seems to work fine), so use id > ipoint
c                                 to identify a refineable point

c                                 modified 3/4/2015 to refine endmember compositions. JADC
c           if (ikp(id).ne.0) then  #refine endmember change
            if (ikp(id).ne.0.and.id.gt.ipoint) then

               quit = .false.

               do j = 1, nsol
                  if (ikp(id).eq.idsol(j)) then 
                     kdsol(j) = kdsol(j) + 1
                     jdsol(kdsol(j),j) = id
                     goto 10
                  end if 
               end do
c                                 new phase, add to list
               nsol = nsol + 1
               idsol(nsol) = ikp(id)
               jdsol(1,nsol) = id
               kdsol(nsol) = 1

            end if 
c                                 new point, add to list
10          npt = npt + 1
            jdv(npt) = i

         end if 

      end do
c DEBUG DEBUG 
      if (lopt(32)) then

         if (npt.eq.icp) then 

            call getmus (1)

         else

            mus = .false. 

         end if

      end if 

      ldv1 = 0
      ldv2 = 0 
      clam = 1d99

      do i = 1, isoct
         slam(i) = 1d99
         kdv(i) = 0 
      end do 
c                                 perp 6.6.3, make a list of metastable
c                                 phases, this list includes two compounds
c                                 and the least metastable composition of
c                                 each solution.      
      do 20 i = 1, jphct
c DEBUG
         if (is(i).ne.1.or.clamda(i).lt.wmach(3)) cycle 
         if (is(i).ne.1) cycle 

         id = i + inc 
         iam = ikp(id)
c                                modified to allow endmembers, Mar 4, 2015. JADC
c        if (iam.ne.0.and.id.gt.0*ipoint) then   #refine endmember change
         if (iam.ne.0.and.id.gt.ipoint) then  

            if (clamda(i).lt.slam(iam)) then
c                                the composition is more stable
c                                than the previous composition 
c                                of the solution. check if it's 
c                                one of the stable solutions                               
               do j = 1, nsol
                  if (iam.eq.idsol(j)) then
c                                it's already stable, only accept
c                                it if its further than the solvus
c                                tolerance from any of the stable
c                                compositions.
                     do k = 1, kdsol(j)
                        if (.not.solvus(jdsol(k,j),id,iam)) goto 20
                     end do 
                  end if
               end do
c                                the composition is acceptable
               slam(iam) = clamda(i)
               kdv(iam) = i

            end if

         else 
c                                a compound, save lowest 2, changed from 1
c                                march 4, 2015, JADC
            if (clamda(i).lt.clam) then 
c                                put the old min into ldv2
               ldv2 = ldv1
c                                and save the current in ldv1
               ldv1 = i

               clam = clamda(i)

            end if 
 
         end if 

20    continue
c                                 load the metastable compounds directly
c                                 into jdv
      if (ldv2.ne.0) then 
         npt = npt + 2
         jdv(npt-1) = ldv2
         jdv(npt) = ldv1
         mcpd = 2
      else if (ldv1.ne.0) then 
         npt = npt + 1
         jdv(npt) = ldv1
         mcpd = 1
      else 
         mcpd = 0 
      end if 
c                                 load the metastable solutions into kdv
      mpt = 0 

      do i = 1, isoct

         if (kdv(i).eq.0) cycle
         mpt = mpt + 1
         kdv(mpt) = kdv(i)
         slam(mpt) = slam(i)

      end do 

      if (mpt.le.iopt(12)) then 
c                                 less metastable refinement points than
c                                 iopt(12)
         max = mpt

      else 
c                                 sort the metastable points to
c                                 find the most stable iopt(12) points
         max = iopt(12)

         call ffirst (slam,kdv,1,mpt,max,h9,ffirst)

      end if 
 
      do i = 1, max

         jdv(npt+i) = kdv(i)
c                                 a metastable solution to be refined
c        if (kdv(i)+inc.gt.ipoint) quit = .false.
         if (ikp(kdv(i)+inc).ne.0) quit = .false.
      end do

      if (quit) then 
c                                 zero mode filter and 
c                                 save amounts for final processing
         mpt = npt - mcpd 
         npt = 0 

         do i = 1, mpt
            if (x(jdv(i)).lt.nopt(9)) cycle 
            npt = npt + 1
            jdv(npt) = jdv(i)
            amt(npt) = x(jdv(i)) 
         end do 

      else 

         npt = npt + max
c                                 sort the phases, why? don't know, but it's 
c                                 necessary
         call sortin 

      end if 

      end 

      recursive subroutine ffirst (a, ind, left, right, k, n, dumsub)
c-----------------------------------------------------------------------
c find the k smallest values of array between indices left and right
c from http://en.wikipedia.org/wiki/Selection_algorithm
c-----------------------------------------------------------------------
      implicit none

      integer left, right, k, n, pivot, opivot, partit, ind(n)

      external dumsub

      double precision a(n)

      if (right.gt.left) then 

         opivot = left + (right-left)/2
         pivot = partit (a, ind, left, right, opivot, n)

         if (pivot.gt.k) then 
            call dumsub (a,ind,left,pivot-1,k,n,dumsub)
         else if (pivot.lt.k) then 
            call dumsub (a,ind,pivot+1,right,k-pivot,n,dumsub)
         end if 

      end if 

      end 

      integer function partit (a, ind, left, right, opivot, n)

      implicit none

      integer left, right, n, pivot, opivot, iold, ind(n), i

      double precision a(n), value, oldval

      value = a(opivot)
c                                 swap a(opivot) with a(right)
      iold = ind(opivot)
      a(opivot) = a(right)
      ind(opivot) = ind(right)
      a(right) = value
      ind(right) = iold

      pivot  = left

      do i = left, right-1

         if (a(i).le.value) then
 
            iold = ind(pivot)
            oldval = a(pivot)
            a(pivot) = a(i)
            ind(pivot) = ind(i)
            a(i) = oldval
            ind(i) = iold
            
            pivot = pivot + 1
         end if 
      end do 
c                                 swap a(right) with a(pivot)
      iold = ind(pivot)
      oldval = a(pivot)
      a(pivot) = a(right)
      ind(pivot) = ind(right)
      a(right) = oldval 
      ind(right) = iold
      
      partit = pivot

      end 

      subroutine subst (a,ipvt,n,b,ier)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'
 
      double precision a(k8,k8),b(k8),x(k8),sum

      integer ipvt(k8),ip,i,j,n,ii,ier
c----------------------------------------------------------------------
c                                 solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)
      do i = 2, n

         sum = 0d0

         do j = 1, i - 1
            sum = a(i,j)*x(j)+sum
         end do 

         ip = ipvt(i)
         x(i) = b(ip)-sum

      end do 
c                                 solve ux = y for x:
      if (a(n,n).eq.0d0) then
c                                 this check should be superfluous,
c                                 but reopt requires it. should check
c                                 what's with factor. 
         ier = 1
         goto 99
      end if 

      x(n) = x(n)/a(n,n)

      do ii = 1, n - 1

         i = n-ii

         sum = 0d0

         do j = i + 1, n
            sum = a(i,j)*x(j)+sum
         end do 

         if (a(i,i).eq.0d0) then
c                                 as above.
            ier = 1
            goto 99
         end if 

         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)

      end do 
      b(n) = x(n)
 
99    end

      subroutine subst1 (n)
c-----------------------------------------------------------------------
c subst uses the lu decomposition of the matrix 'a' contained
c in the array 'a' to solve ax = b for x. subst is modified from the
c the subroutine of the same name listed by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
 
c input     a- an n by n array containing the non-zero elements of
c              the u and l decompositions of a, as output by factor.
c           n- the dimension of the matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the coefficient a(n,k).
c           b- the vector b.
c output    b- the solution vector x.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      double precision x(k5), sum

      integer n, i, j, im1, ip1, nm1, ii, ip

      integer ipvt
      double precision a,b
      common/ cst301 /a(k5,k5),b(k5),ipvt(k5) 

c                            solve ly = b for y:
      ip = ipvt(1)
      x(1) = b(ip)

      do i = 2, n
         sum = 0d0
         im1 = i-1
         do j = 1, im1
            sum = a(i,j)*x(j)+sum
         end do 
         ip = ipvt(i)
         x(i) = b(ip)-sum
      end do 
c                            solve ux = y for x:
      x(n) = x(n)/a(n,n)
      nm1 = n-1

      do ii = 1, nm1
         i = n-ii
         ip1 = i+1
         sum = 0d0
         do j = ip1, n
            sum = a(i,j)*x(j)+sum
         end do
         x(i) = (x(i)-sum)/a(i,i)
         b(i) = x(i)
      end do 

      b(n) = x(n)
 
      end

      subroutine factr1 (n,ier)
c-----------------------------------------------------------------------
c factr1 is a subroutine which calculates the triangular
c decompositions of the matrix 'a'. factor is modified from
c the subroutine of the same name given by conte and de boor
c in 'elementary numerical analysis', mcgraw-hill, 1980.
c factor uses scaled partial pivoting.
c
c input     a- an n by n array containing the elements of matrix a.
c           n- the dimension of the matrix a.
c output    a- an n by n array containing the upper, u, and lower, l,
c              triangular decompositions of input matrix a.
c        ipvt- a vector indicating that row ipvt(k) was used to
c              eliminate the a(n,k).
c         ier- a flag, zero if a is of rank = n, and 1 if a is of
c              lower rank.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,n,ip1,istr,ier

      double precision temp,ratio,tmax,rmax

      integer ipvt
      double precision a,d,x(k5)
      common/ cst301 /a(k5,k5),d(k5),ipvt(k5)
c-----------------------------------------------------------------------
      ier = 0
c                            initialize ipvt,d
      do i = 1, n
         ipvt(i) = i
         rmax = 0d0
         do j = 1,n
            rmax = dmax1(rmax,dabs(a(i,j)))
         end do 
c                            ax = b is singular if rmax = 0
         if (dabs(rmax).lt.1d-5) goto 9000
         x(i) = rmax
      end do 
c                            begin decomposition:
      do i = 1, n - 1
c                            determine pivot row (istr).
         rmax = dabs(a(i,i))/x(i)
         istr = i
         ip1 = i + 1

         do j = ip1, n
            tmax = dabs(a(j,i))/x(j)
            if (tmax.le.rmax) cycle
            rmax = tmax
            istr = j
         end do 

         if (dabs(rmax).lt.1d-5) goto 9000
c                            if istr gt i, make i the pivot row
c                            by interchanging it with row istr.
         if (istr.gt.i) then 
            j = ipvt(istr)
            ipvt(istr) = ipvt(i)
            ipvt(i) = j
            temp = x(istr)
            x(istr) = x(i)
            x(i) = temp
            do j = 1, n
               temp = a(istr,j)
               a(istr,j) = a(i,j)
               a(i,j) = temp
            end do 
         end if 
c                            eliminate x(k) from rows k+1,...,n.
         do j = ip1,n
            a(j,i) = a(j,i)/a(i,i)
            ratio = a(j,i)
            do k = ip1, n
               a(j,k) = a(j,k)-ratio*a(i,k)
            end do 
         end do 
 
      end do 
     
      if (dabs(a(n,n)).lt.1d-5) ier = 1

      return
c                           algoritmic singularity.
9000  ier = 1
 
      end

      subroutine ufluid (fo2)
c----------------------------------------------------------------------
c subroutine ufluid computes the potential of the components
c of a saturated fluid phase. if the mole fraction of a component is les
c less than 1.d-38 the chemical potential is set to -9.9d09.
c ufluid may call one of three molecular fluid equations of state, or
c alternatively users may supply their own routines, however,
c the routines currently in use return the log of a components fugacity
c which is then added to the reference state potential computed by the
c function gphase.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i
 
      double precision xf(2), fo2, fs2, gcpd, gzero 

      external gcpd, gzero

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      double precision f
      common/ cst11 /f(3)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug
c-----------------------------------------------------------------------
c                           compute the chemical potentials of
c                           fluid components in fluid saturated
c                           systems.
      call cfluid (fo2,fs2)

      if (idfl.ne.0) then
   
         uf(idfl) = gcpd (idfl,.false.) + r * t * f(idfl)

      else

         xf(1) = 1d0 - xco2
         xf(2) = xco2
 
         do i = 1, 2

            if (iff(i).ne.0) then 

               if (xf(i).lt.1d-38) then 

                  uf(i) = -1d10

               else 

                  uf(i) = gzero (i) + r * t * f(i)

               end if

            end if 

         end do

      end if 
 
      end

      subroutine uproj
c----------------------------------------------------------------------
c subroutine uproj computes the potentials of saturated phase components
c and saturated components.
  
c the energies of saturated components are projected through
c saturated volatile components.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
 
      integer i,j,k,l,ict,ll,i1,id

      double precision uss(h6), fo2, gph, u, gphase

      external gphase

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c-----------------------------------------------------------------------
      fo2 = 0d0
c                                 compute the chemical potentials of
c                                 saturated phase components.
      if (ifct.gt.0) call ufluid (fo2)

      do i = 1, isat
c                                 determine stable saturated composants
c                                 and the corresponding chemical potentials
         ict = isct(i)

         ll = icp+i

         do j = 1, ict

            k = ids(i,j)
            gph = gphase (k)
            
            if (ifct.gt.0) then 
               do l = 1, 2
c                                 legendre transform for saturated phase
c                                 component potentials
                  if (iff(l).ne.0) gph = gph - cp(iff(l),k)*uf(l)
               end do 
            end if 

            uss(j) = gph 

            if (i.gt.1) then 
c                                 if multiple component saturation constraints
c                                 apply saturation hierarchy legendre transform:
               i1 = i-1
               do l = 1, i1
                  uss(j) = uss(j)-cp(icp+l,k)*us(l)
               end do
            end if 

            g(k) = uss(j)
            uss(j) = uss(j)/cp(ll,k)
         end do 
c                                 if O2, check if fo2 has been 
c                                 determined by a fluid phase routine,
c                                 if so, add the transform:
         if (io2.eq.i) then 
            do j = 1, ict 
               uss(j) = uss(j) + r*t*fo2
            end do 
         end if 
c                           now find stable "composant":

         u = uss(1)

         id = 1

         if (ict.ne.1) then 
            do j = 2, ict
               if (uss(j).gt.u) cycle  
               id = j
               u = uss(j)
            end do 
         end if 
c                               save the id of the stable composant.
         idss(i) = ids(i,id)
c                               and its chemical potential.
         us(i) = u
c                               in case a phase in the component
c                               saturation space is an endmember of
c                               a solution transform the endmember G's:
         do j = 1, ict
            k = ids(i,j)
            g(k) = g(k) - cp(icp+i,k)*u
         end do 

      end do 

      end

      subroutine initlp 
c--------------------------------------------------------------------
c initialize arrays and constants for lp minimizarion
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,j,k,l,inc,id,jk

      logical bad

      double precision u,s,vol,tot 

      double precision ctot
      common/ cst3  /ctot(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp  

      integer jphct,istart
      common/ cst111 /jphct,istart

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk

      integer ldt,ldq
      common /be04nb/ldt,ldq

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer ipoint,kphct,imyn
      common/ cst60  /ipoint,kphct,imyn
c-----------------------------------------------------------------------
      inc = istct - 1

      tloop = 40
      ploop = 40
      dv(1) = (vmax(1)-vmin(1))/(ploop-1)
      dv(2) = (vmax(2)-vmin(2))/(tloop-1)
c                                 load arrays for lp solution
      jphct = iphct - inc

      if (.not.usv) then
c                                 pressure and temperature are allowed 
c                                 EoS variables
         ctotal = 0d0

         do i = 1, icp
            ctotal = ctotal + cblk(i)
         end do 
c                                 composition constraint
         do i = 1, icp
            b(i) = cblk(i)/ctotal
         end do 

         do i = 1, jphct
            id = i + inc
            iam(i) = id
            do j = 1, icp
               a(j,i) = cp(j,id)/ctot(id)
            end do
         end do

         ldt = icp + 1
         ldq = icp + 1

      else
c                                 using usv formulation
         kphct = 0 

         do i = 1, jphct

            id = i + inc
            jk = 0
c                                 for each static chemical cpd
            do j = 1, tloop
c                                 for each temperature
               do k = 1, ploop
c                                 for each pressure
c                 
                  jk = jk + 1  

                  call pt4sv (jk)  

                  call getusv (u,s,vol,id,bad)
  
          
                  kphct = kphct + 1
                  
                  tot = ctot(id) + s + vol
                  
                  do l = 1, jbulk
                    a(l,kphct) = cp(l,id)/tot
                  end do

                  a(icp+1,kphct) = s/tot
                  a(icp+2,kphct) = vol/tot
                  c(kphct) = u/tot
                  iam(kphct) = id
                  jam(kphct) = jk

               end do 

            end do 

         end do 

         jphct = kphct

         ldt = jbulk + 1
         ldq = jbulk + 1

      end if 
c                                 cold start istart = 0
      istart = 0

      end 

      subroutine gall 
c-----------------------------------------------------------------------
c subroutine gall computes molar free energies of all static compounds.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'
c DEBUG
      integer i, j, k, id, isp, ins(nsp), ibug
      double precision gsol0
      external gsol0

      double precision dg1, gval, dg, gzero, g0(k5), gex, gfesi,
     *                 gfesic, gfecr1, gerk, x1(5), gproj, ghybrid

      external gerk, gzero, gex, gfesi, gfesic, gproj, ghybrid

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision g
      common/ cst2 /g(k1)

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision thermo,uf,us 
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      integer ifct,idfl
      common/ cst208 /ifct,idfl

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug

      integer isoct
      common/ cst79 /isoct

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      integer ixp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1)

      integer jend
      common/ cxt23 /jend(h9,m4)

      integer jspec
      common/ cxt8 /jspec(h9,m4)
c                                 working arrays
      double precision z, pa, p0a, x, w, y, wl
      common/ cxt7 /x(m4),y(m4),pa(m4),p0a(m4),z(mst,msp),w(m1),
     *              wl(m17,m18)
c                                 model type
      logical lorder, lexces, llaar, lrecip
      common/ cxt27 /lorder(h9),lexces(h9),llaar(h9),lrecip(h9)

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      integer ksmod, ksite, kmsol, knsp
      common/ cxt0  /ksmod(h9),ksite(h9),kmsol(h9,m4,mst),knsp(m4,h9)

      integer ideps,icase,nrct
      common/ cxt3i /ideps(j4,j3,h9),icase(h9),nrct(j3,h9)
c                                 endmember names
      character names*8
      common/ cst8  /names(k1)

      save ibug
      data ibug/1/
c-----------------------------------------------------------------------
c                                 compute the chemical potential
c                                 of the projected components.
      call uproj
c                                 first do the endmembers:
      do id = kphct+1, ipoint

         g(id) = gproj (id) 

      end do 
c                                 now do solutions:
      do i = 1, isoct
c                                 check if normal solution:
         if (.not.llaar(i).and.(ksmod(i).eq.7.or.ksmod(i).eq.5.or.
     *       ksmod(i).eq.2.or.ksmod(i).eq.24.or.ksmod(i).eq.25.or.
     *       ksmod(i).eq.28)) then 
c                                 it's normal margules or ideal:
            do j = 1, jend(i,2)
c                                 initialize with excess energy, dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               do k = 1, lstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 

               id = id + 1

            end do 

         else if (ksmod(i).eq.0) then
c                                 it's a fluid compound, the way 
c                                 things are now it must have two
c                                 components.
            do j = 1, mstot(i)
               g0(j) = gzero(jend(i,2+j))
            end do 

            do j = 1, jend(i,2)

               call fexces (id,gval)

               g(id) = g0(1) * sxs(ixp(id)+1) + g0(2) * sxs(ixp(id)+2) 
     *                                        + gval
               id = id + 1

            end do 

         else if (ksmod(i).eq.30.or.ksmod(i).eq.31) then 
c                                 call nastia's BCC/FCC model
c                                 after Lacaze and Sundman
            do j = 1, jend(i,2)

               g(id) = gfesic (sxs(ixp(id)+1),
     *                         sxs(ixp(id)+3),sxs(ixp(id)+4),
     *                         g(jend(i,3)),g(jend(i,4)),
     *                         g(jend(i,5)),g(jend(i,6)),ksmod(i))
  
               id = id + 1

            end do

         else if (lrecip(i).and.lorder(i)) then
c                                 prismatic solution with ordering 
c            sxs(1) = .6407726359
c            sxs(2) = .4912112378e-1
c            sxs(3) = .1781709887e-1
c            sxs(4) = .2657701269
c            sxs(5) = .2972736393e-1
c            sxs(6) = -.320834956e-2
c eq
c            sxs(1) = .7012131600d0
c            sxs(2) = .513d-1
c            sxs(3) = .1882040400d0
c            sxs(4) = .2473868400d0
c            sxs(5) = .5216d-4
c            sxs(6) = -.1881562000d0

c                                 compute margules coefficients
            call setw (i)
c                                 compute enthalpy of ordering
            call oenth (i)
c                                 now for each compound:
            do j = 1, jend(i,2)


           ibug = 99
              if (ibug.eq.0) then 
c DELETE ME DELETE ME DEBUG 1050 K 4000
         p0a(1) =                   -0.0065200000
         p0a(2) =                    0.2235000000
         p0a(3) =                    0.4769502876
         p0a(4) =                    -.21934
         p0a(5) =                    0.0928182680
         p0a(6) =                    0.5919781530
         p0a(7) =                      0.24302
         p0a(8) =                    0.0354697124
         p0a(9) =                      0.07056
         p0a(10) =                   -0.3901582680
         p0a(11) =                   -0.1182781530
               do k = 1, nstot(i) 
                  pa(k) = p0a(k)
               end do 
         dg = gsol0(i)

             else if (ibug.eq.2) then 

         p0a(1) =  .8124e-1              
         p0a(2) =   .79822           
         p0a(3) =   .4341293766            
         p0a(4) =   -.6628e-1             
         p0a(5) =  .9605e-1                 
         p0a(6) =    .464433782e-1              
         p0a(8) =  -.3898027548         
         p0a(7) =   0.                

               do k = 1, nstot(i) 
                  pa(k) = p0a(k)
               end do 
         dg = gsol0(i)

            end if 
c                                 assign x's
               do k = 1, nstot(i) 
c DELETE ME DELETE ME DEBUG 
                  p0a(k) = sxs(ixp(id)+k)
                  pa(k) = p0a(k)
               end do 
c                                 get the speciation energy effect
               call specis (dg,i)
c                                 and internal endmember dqf
               call gexces (id,dg1)

               g(id) = dg + dg1 
c                                 add in g from real endmembers, this
c                                 must include the g for the disordered equivalent
c                                 of the ordered species
               do k = 1, lstot(i)

                  g(id) = g(id) + g(jend(i,2+k)) * p0a(k)

               end do 

               id = id + 1

            end do 

         else if (lorder(i)) then
c                                 compute margules coefficients
            call setw (i)
c                                 compute enthalpy of ordering
            call oenth (i)
c                                 now for each compound:
            do j = 1, jend(i,2)
c                                 for speciation models gexces
c                                 evaluates only endmember sconf 
c                                 and internal dqf's
               call gexces (id,g(id))
c                                 gmech
               do k = 1, lstot(i)
                  p0a(k) = sxs(ixp(id)+k)
                  pa(k) = p0a(k)
                  g(id) = g(id) + g(jend(i,2+k)) * p0a(k)
               end do

               call specis (dg,i)

               g(id) = g(id) + dg

               id = id + 1

            end do 

         else if (llaar(i)) then 
c                                 compute margules coefficients
            call setw (i)
c                                 because the hp van laar may have p-t
c                                 dependent volumes, the full expression
c                                 must be evaluated here:                   
            do j = 1, jend(i,2)
c                                 initialize with dqf,
c                                 and configurational entropy terms
               call gexces (id,g(id))

               do k = 1, lstot(i) 
                  x(k) =  sxs(ixp(id)+k)
                  g(id) = g(id) + g(jend(i,2+k)) * x(k)
               end do 
c                                 add the real excess energy
               g(id) = g(id) + gex(i,x)

               id = id + 1

            end do 

         else if (ksmod(i).eq.23) then
c                                 toop melt:
            do j = 1, nstot(i)
               g0(j) = gzero (jend(i,2+j))
            end do 

            do j = 1, jend(i,2)
               call toop (id,g(id))
               do k = 1, nstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 
               id = id + 1
            end do 

         else if (ksmod(i).eq.26) then
c                                 H2O-CO2-Salt:
            do j = 1, jend(i,2)
        
               call hcneos (g(id),sxs(ixp(id)+1),
     *                      sxs(ixp(id)+2),sxs(ixp(id)+3))

               do k = 1, nstot(i) 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 

               id = id + 1

            end do

         else if (ksmod(i).eq.39) then
c                                 generic hybrid EoS
c                                 initialize pointer array
            isp = mstot(i)

            do j = 1, isp
               ins(j) = jspec(i,j)
            end do 

            do j = 1, jend(i,2)

               g(id) = 0d0

               do k = 1, isp
c                                 load composition array and pointers 
                  x(k) = sxs(ixp(id)+k)
c                                 sum pure species g's
                  g(id) = g(id) + g(jend(i,2+k)) * x(k)

               end do
c                                 compute and add in activities
               g(id) = g(id) + ghybrid (x,ins,isp)

               id = id + 1

            end do 

         else if (ksmod(i).eq.29) then 
 
            do j = 1, jend(i,2)
c                                 BCC Fe-Si Lacaze and Sundman
               g(id) = gfesi(sxs(ixp(id)+1),g(jend(i,3)),g(jend(i,4)))
  
               id = id + 1

            end do 

         else if (ksmod(i).eq.32) then 
 
            do j = 1, jend(i,2)
c                                 BCC Fe-Cr Andersson and Sundman
               g(id) = gfecr1 (sxs(ixp(id)+1),g(jend(i,3)),g(jend(i,4)))

               id = id + 1

            end do  


         else if (ksmod(i).eq.41) then 

            do j = 1, jend(i,2)
c                                 hybrid MRK ternary COH fluid
               call rkcoh6 (sxs(ixp(id)+2),sxs(ixp(id)+1),g(id)) 

               do k = 1, 3 
                  g(id) = g(id) + g(jend(i,2+k)) * sxs(ixp(id)+k)
               end do 

               id = id + 1

            end do 

         else if (ksmod(i).eq.40) then 

            do j = 1, jend(i,2)
c                                 MRK silicate vapor
               g(id) = 0d0

               do k = 1, nstot(i) 
                  x1(k) = sxs(ixp(id)+k)
                  g(id) = g(id) + gzero(jend(i,2+k)) * x1(k)
               end do 

               g(id) = g(id) + gerk(x1)   

               id = id + 1

            end do  

         end if 

      end do 

      end

      logical function solvus (id1,id2,ids)
c-----------------------------------------------------------------------
c function to test if a solvus separates two static pseudocompounds of
c solution ids.
c-----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i, id1, id2, ids

      integer icomp,istct,iphct,icp
      common/ cst6 /icomp,istct,iphct,icp 

      double precision cp
      common/ cst12 /cp(k5,k1)

      double precision dcp,soltol
      common/ cst57 /dcp(k5,k19),soltol

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c-----------------------------------------------------------------------
      solvus = .false.

      do i = 1, icp

         if (dcp(i,ids).eq.0d0) cycle 

         if (dabs(cp(i,id1)-cp(i,id2))/dcp(i,ids).gt.soltol) then 
            solvus = .true.
            exit
         end if 

      end do 

      end 

      subroutine grxn (gval) 
c-----------------------------------------------------------------------
c grxn computes the free energy of univariant equilibria
c defined by the data in commonn block cst21 which is initialized
c in the subprogram balanc.  grxn is partially redundant with
c the function gphase but because of the frequency that these
c these routines are used a significant increase in efficiency is
c gained by maintaining separate functions.
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer j

      double precision gval, gproj

      external gproj

      integer iffr,isr
      double precision vuf,vus
      common/ cst201 /vuf(2),vus(h5),iffr,isr

      integer idr,ivct
      double precision vnu
      common/ cst25 /vnu(k7),idr(k7),ivct
c-----------------------------------------------------------------------
c                                 compute potentials of saturated phases
c                                 and components, note that in this
c                                 version of vertex the stoichiometry of
c                                 such components may vary.

c                                 no saturated phase components and no
c                                 saturated components:
      if (iffr.eq.1.and.isr.eq.1) goto 10
c                                 note that this call to uproj makes a
c                                 subsequent call in gall redundant if
c                                 sfol1 is used to trace a univariant
c                                 curve.
      call uproj
c                                 compute free energy change of the rxn
10    gval = 0d0

      do j = 1, ivct
         gval = gval + vnu(j) * gproj (idr(j))
      end do 

      end

      subroutine yclos0 (x,is,jphct)
c----------------------------------------------------------------------
c subroutine to save optimization results for non-iterative refinement
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i,jphct,is(*)

      double precision x(*)
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)
c----------------------------------------------------------------------

      npt = 0 

      do i = 1, jphct

         if (is(i).eq.1.or.x(i).lt.nopt(9)) cycle  
c                                 acceptable cases 0 active, between bounds
c                                                  2 active, upper bound 
            npt = npt + 1
            jdv(npt) = i 
            amt(npt) = x(i)
 
      end do
    
      end 

      subroutine setx3 (ind,id,ids)
c----------------------------------------------------------------------
c subroutine to recover geometric reciprocal solution compositions (x(i,j))
c from the xcoor (reciprocal) or sxs (single site) arrays loaded in soload
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ids, jcoor, ind
c                                 x coordinate description
      integer istg, ispg, imlt, imdg

      common/ cxt6i /istg(h9),ispg(h9,mst),imlt(h9,mst),imdg(ms1,mst,h9)
c                                 stored x coordinate
      double precision xcoor
      integer icoor
      common/ cxt10 /xcoor(k18),icoor(k1)
c                                 single site solution coordinates:
      integer ixp
      double precision sxs,exces
      common/ cst304 /sxs(k13),exces(m3,k1),ixp(k1)
c                                  x-coordinates for the final solution
      double precision x3
      common/ cxt16 /x3(k21,mst,msp)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      integer lstot,mstot,nstot,ndep,nord
      common/ cxt25 /lstot(h9),mstot(h9),nstot(h9),ndep(h9),nord(h9)
c----------------------------------------------------------------------
      if (id.gt.ipoint) then 
c                                 a normal solution
         if (istg(ids).eq.1) then 
c                                 one site solution
            do j = 1, nstot(ids)
               x3(ind,1,j) = sxs(ixp(id)+j) 
            end do 

         else if (ispg(ids,1).gt.1) then 
c                                 multi-site solution
            jcoor = icoor(id)

            do i = 1, istg(ids)
               do j = 1, ispg(ids,i)
                  jcoor = jcoor + 1
                  x3(ind,i,j) = xcoor(jcoor)
               end do 
            end do 

         else 
c                                 a dummy site
            x3(ind,1,1) = 1d0
          
            do j = 1, nstot(ids)
               x3(ind,2,j) = sxs(ixp(id)+j) 
            end do 

         end if 

      else 
c                                 an endmember 
         call endcp (ind,id,ids)

      end if

      end 

      double precision function ginc0 (dt,dp,id)
c-----------------------------------------------------------------------
c id indicates a static pseudocompound, use ginc for dynamic compound
c-----------------------------------------------------------------------
      implicit none

      double precision dt, dp, gphase

      external gphase

      integer id

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps
c----------------------------------------------------------------------
      p = p + dp 
      t = t + dt 

      ginc0 = gphase (id)

      p = p - dp 
      t = t - dt

      end 

      subroutine getusv (u,s,v,id,bad)
c-----------------------------------------------------------------------
c getusv computes u,s,v by centered finite differences from the Gibbs energy

c the difference increments are

c dt0, dp0 for 1st order derivatives (entropy,volume and enthalpy)
c-----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer id

      logical bad

      double precision dt0,dp0,dp1,dp2,g0,v,ginc0,s,u

      double precision p,t,xco2,u1,u2,tr,pr,r,ps
      common/ cst5 /p,t,xco2,u1,u2,tr,pr,r,ps

      save dt0
      data dt0/0.5d0/
c----------------------------------------------------------------------

      dp0 = 0.5d-3 * p 
      dp1 = 0.5d-2 * p
      dp2 = 0.5d-1 * p
            
      g0 = ginc0(0d0,0d0,id)
c                                 straight derivatives:
c                                 first order
      if (p-dp0.le.0d0) then 

         v = (ginc0(0d0,dp0,id) - g0)/dp0
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment if invalid v
     *   v = (ginc0(0d0,dp1,id) - g0)/dp1
         if (v.lt.0d0.or.dabs(v).gt.1d9)  
c                                 expand increment more if invalid v
     *   v = (ginc0(0d0,dp2,id) - g0)/dp2

      else 

         v = (ginc0(0d0,dp0,id) - ginc0(0d0,-dp0,id))/dp0/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp1.gt.0d0)  
c                                 expand increment if invalid v
     *   v = (ginc0(0d0,dp1,id) - ginc0(0d0,-dp1,id))/dp1/2d0
         if ((v.lt.0d0.or.dabs(v).gt.1d9).and.p-dp2.gt.0d0)  
c                                 expand increment more if invalid v
     *   v = (ginc0(0d0,dp2,id) - ginc0(0d0,-dp2,id))/dp2/2d0

      end if 
c                                 in case the evaluating routine fails
c                                 on both calls to ginc 
      if (v.eq.0d0) v = 1d0

      s = (ginc0(-dt0,0d0,id) - ginc0(dt0,0d0,id))/dt0/2d0

      u = g0 + t*s - p*v
 
      if (s.lt.0d0.or.v.lt.0d0) then 
         bad = .true.
      else
         bad = .false.
      end if 

      end 

      subroutine pt4sv (ind)
c--------------------------------------------------------------------
c subroutine to recover the p-t condition of discreitization in USV
c calculations.
c---------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer ind,tind,pind

      double precision v,tr,pr,r,ps
      common / cst5 /v(l2),tr,pr,r,ps

      double precision vmax,vmin,dv
      common/ cst9 /vmax(l2),vmin(l2),dv(l2)

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c----------------------------------------------------------------------

      tind = (ind-1)/tloop
      pind = ind - tind*tloop - 1
      v(2) = vmin(2) + tind*dv(2)
      v(1) = vmin(1) + pind*dv(1)

      end 

      subroutine yclos2 (clamda,x,is,iter,opt)
c----------------------------------------------------------------------
c subroutine to identify pseudocompounds close to the solution for 
c subsequent refinement, for iteration > 1. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      external ffirst

      integer i, is(*), id, jmin(k19), kmin(k19), opt, kpt, mpt, iter, 
     *        tic

      double precision clamda(*), clam(k19), x(*)

      logical stable(k19)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      double precision wmach(9)
      common /ax02za/wmach

      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer hkp,mkp
      common/ cst72 /hkp(k21),mkp(k19)

      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot

      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct

      save tic
      data tic/0/
c----------------------------------------------------------------------
c                                 npt is the number of points refined
c                                 from the previous cycle, opt is the
c                                 number of points in the original 
c                                 solution.
      do i = 1, k19
         jmin(i) = 0 
         clam(i) = 1d99
         stable(i) = .false.
      end do 

      npt = 0
      mpt = 0 
      kpt = 0

      do i = 1, jphct
c                                 id indicates the original refinement
c                                 point.
         id = hkp(i)
c                                 check the stability of all points 
         if (is(i).ne.1.and.x(i).gt.0d0) then 
c                                 a stable point, add to list
            npt = npt + 1
            jdv(npt) = i
            stable(id) = .true.

         else if (clamda(i).lt.clam(id)) then
c DEBUG
            if (clamda(i).lt.wmach(3)) cycle 

            if (jkp(id).gt.0) then
c                                 it's a solution, keep the  
c                                 least metastable point
               jmin(id) = i
               clam(id) = clamda(i)
               mpt = mpt + 1
  
            else 
c                                 it's a compound, keep all 
c                                 as they hardly cost anything. 
               kpt = kpt + 1
               kmin(kpt) = i

            end if 

         end if 

      end do

c DEBUG DEBUG

      if (lopt(32)) then 

         if (npt.eq.hcp) then 

            call getmus (iter)

         else

           mus = .false. 

         end if

      end if 

      if (iter.le.iopt(10)) then
c                                 if not done iterating, add the metastable
c                                 phases, first the compounds:
         do i = 1, kpt
            npt = npt + 1 
            jdv(npt) = kmin(i)
         end do 
c                                 make a list of the solutions
         kpt = 0 

         do i = 1, opt

            if (jmin(i).eq.0) then
               cycle
            else if (stable(hkp(jmin(i)))) then 
               cycle 
            end if 

            kpt = kpt + 1
            jmin(kpt) = jmin(i)
            clam(kpt) = clam(i)

         end do
c DEBUG DEBUG 
c should be gt 0
         if (kpt.gt.iopt(31)) then 
c                                 sort if too many

            kpt = iopt(31)

            call ffirst (clam,jmin,1,kpt,iopt(31),k19,ffirst)

         end if 

         do i = 1, kpt

           npt = npt + 1
           jdv(npt) = jmin(i)

         end do 
c                                 sort the phases, this is only necessary if
c                                 metastable phases have been added
         call sortin
c                                 make a pointer to the original refinement 
c                                 point
         do i = 1, npt
            mkp(i) = hkp(jdv(i))
         end do 

      else  
  
         mpt = npt 
         npt = 0  
c                                 check zero modes the amounts
         do i = 1, mpt

            if (x(jdv(i)).ge.nopt(9)) then 
               npt = npt + 1
               amt(npt) = x(jdv(i))
               jdv(npt) = jdv(i)
            else if (lopt(13).and.x(jdv(i)).lt.-nopt(9)
     *                             .and.tic.lt.5) then 

               call warn (2,x(jdv(i)),i,'YCLOS2')
               tic = tic + 1

               if (tic.eq.5) call warn (49,x(1),2,'YCLOS2')

            end if 

         end do 

         if (npt.ne.icp) then
 
            fulrnk = .false.
         
         else 

            fulrnk = .true.

         end if 

      end if 

      end

      subroutine getmus (iter) 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, j, id, ier, ipvt(k8), iter

      double precision comp(k8,k8)

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------

       mus = .false.

       if (iter.gt.1) then 

         do i = 1, hcp

            id = jdv(i) 

            do j = 1, hcp 
               comp(i,j) = cp2(j,id)
            end do 

            mu(i) = g2(id)

         end do

      else

         do i = 1, hcp

            id = jdv(i) 

            do j = 1, hcp 
               comp(i,j) = a(j,id)
            end do 

            mu(i) = c(id)

         end do

      end if 

      call factor (comp,hcp,ipvt,ier)

      if (ier.ne.0) return

      call subst (comp,ipvt,hcp,mu,ier)

      if (ier.eq.0) mus = .true. 

      end 

      double precision function deltag (id,iter)
c----------------------------------------------------------------------
c return the difference g phase id - g current solution, i.e., < 0 phase
c id is metastable. 
c----------------------------------------------------------------------
      implicit none

      include 'perplex_parameters.h'

      integer i, id, iter

      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      double precision a,b,c
      common/ cst313 /a(k5,k1),b(k5),c(k1)

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk
c----------------------------------------------------------------------
      if (iter.gt.0) then 

         deltag = g2(id)

         do i = 1, hcp

            deltag = deltag - cp2(i,id)*mu(i)

         end do

      end if
 
      end 

      subroutine rebulk (static)
c----------------------------------------------------------------------
c upon successful completion of an optimization with either static or
c dynamic pseudocompounds rebulk:
c     1) loads the generic arrays cp3, cptot, ctot3, kkp and x3
c     2) computes the amounts of saturated component phases.
c     3) computes the dependent potentials.
c     4) checks for solvi and homogenizes miscible phases.
c----------------------------------------------------------------------
      implicit none
 
      include 'perplex_parameters.h'

      integer i,j,k,id,ier,tictoc,ipvt(k8)

      logical static, bad

      double precision c(k5),u,comp(k8,k8)

      integer iff,idss,ifug
      common/ cst10  /iff(2),idss(h5),ifug
 
      integer icomp,istct,iphct,icp
      common/ cst6  /icomp,istct,iphct,icp

      double precision cp
      common/ cst12 /cp(k5,k1)

      integer jbulk
      double precision cblk
      common/ cst300 /cblk(k5),jbulk
c                                 adaptive x(i,j) coordinates
      integer jcoct, jcoor, jkp
      double precision zcoor
      common/ cxt13 /zcoor(k20),jcoor(k21),jkp(k21),jcoct
c                                 compositions of stable adaptive
c                                 coordinates (and solution ids).
      integer kkp,np,ncpd,ntot
      double precision cp3,amt
      common/ cxt15 /cp3(k0,k19),amt(k19),kkp(k19),np,ncpd,ntot
c                                 options from perplex_option.dat
      integer iopt
      logical lopt
      double precision nopt
      common/ opts /nopt(i10),iopt(i10),lopt(i10)

      integer ikp
      common/ cst61 /ikp(k1)

      integer npt,jdv
      logical fulrnk
      double precision cptot,ctotal
      common/ cst78 /cptot(k19),ctotal,jdv(k19),npt,fulrnk

      integer iam,jam,tloop,ploop
      common/ cst55 /iam(k1),jam(k1),tloop,ploop
c                                 hcp is different from icp only if usv
      integer hcp,idv
      common/ cst52  /hcp,idv(k7) 

      logical mus
      double precision mu, gmax
      common/ cst330 /mu(k8),gmax,mus

      integer jtest,jpot
      common/ debug /jtest,jpot

      integer ids,isct,icp1,isat,io2
      common/ cst40 /ids(h5,h6),isct(h5),icp1,isat,io2

      logical usv
      integer pindex,tindex
      common/ cst54 /pindex,tindex,usv

      double precision ctot
      common/ cst3  /ctot(k1)

      integer ipoint,kphct,imyn
      common/ cst60 /ipoint,kphct,imyn

      double precision thermo,uf,us
      common/ cst1 /thermo(k4,k10),uf(2),us(h5)

      double precision g
      common/ cst2 /g(k1)

      integer jphct
      double precision g2, cp2
      common/ cxt12 /g2(k21),cp2(k5,k21),jphct

      save tictoc
      data tictoc/0/
c----------------------------------------------------------------------
      do i = 1, npt

         if (static) then 

            id = iam(jdv(i))
c                                 set identifier flag
            if (id.le.ipoint) then
               kkp(i) = -id
            else  
               kkp(i) = ikp(id)
            end if 

            cptot(i) = ctot(id)
c                                 save g for potential calculation
            mu(i) = g(id)
           
            do j = 1, icomp
               if (j.gt.icp.and.usv) exit
               cp3(j,i) = cp(j,id)
            end do
c                                 set the x3 array
            if (ikp(id).ne.0) call setx3 (i,id,ikp(id))

            if (usv) then 
c                                 usv calculations are currently
c                                 not set up for solutions and use
c                                 only static optimization
               call pt4sv (jam(jdv(i)))

               call getusv (u,cp3(icp1,i),cp3(jbulk,i),id,bad)

            end if 

         else 
c                                 getcmp assigns cp3, cptot, x3, and kkp
            call getcmp (i,jdv(i),jkp(jdv(i)))
c                                 convert normalized g's to molar g's,
c                                 before 3/5/2015 this was outside the
c                                 if construction. JADC.
            mu(i) = g2(jdv(i))*cptot(i)

         end if          
c                                 convert normalized amounts to molar 
c                                 amounts
         amt(i) = ctotal*amt(i)/cptot(i)


      end do 

      if (jbulk.gt.icp.and.(.not.usv)) then  
c                                 get the amounts of the saturated phases:
         do i = icp+1, jbulk
c                                 k is the saturated component pointer
            k = i - icp
c                                 initialize bulk                                 
            c(k) = cblk(i)
c                                 save chemical potentials from gproj
            mu(i) = us(k)
c                                 subtract the amount of the component in the 
c                                 phases in the thermodynamic c-space
            do j = 1, npt 
               c(k) = c(k)- amt(j)*cp3(i,j)
            end do 

         end do 

         do i = jbulk, icp+1, -1
c                                  cycle through the saturated phases
            npt = npt + 1
            id = idss(i-icp)
c                                  set case for solution in saturated component
c                                  space, the endmember composition is not set,
c                                  this is gonna cause problems, at least for 
c                                  meemum
            if (ikp(id).eq.0) then 
               kkp(npt) = -id
            else 
               write (*,1000) 
               stop 
            end if 
c                                  amount of the staurated phase
            amt(npt) = c(i-icp)/cp(i,id)
c                                  warn on undersaturation
            if (amt(npt).lt.nopt(9)) then 
               if (amt(npt).lt.-nopt(9).and.tictoc.lt.5) 
     *                             call warn (2,amt(npt),i,'REBULK')
               npt = npt - 1
               tictoc = tictoc + 1
               exit
            end if 
c                                  remove the saturated phase from 
c                                  the bulk composition.
            do j = icp+1, i - 1
               c(j) = c(j) - amt(npt)*cp(j,id)
            end do 
c                                  load the saturated phase composition 
            do j = 1, icomp
               cp3(j,i) = cp(j,id)
            end do           

         end do

      end if 

      ntot = npt

      if (usv.or.jpot.eq.0) then
c                                 compute chemical potentials
         if (npt-(jbulk-icp).ne.hcp) then 
c                                 not full rank
            do i = 1, hcp
               mu(i) = nopt(7)
            end do
          
         else 

            do i = 1, hcp
               do j = 1, hcp 
                  comp(i,j) = cp3(j,i)
               end do 
            end do 

            call factor (comp,hcp,ipvt,ier)

            if (ier.eq.1) then 

               do i = 1, hcp
                  mu(i) = nopt(7)
               end do
     
            else 
 
               call subst (comp,ipvt,hcp,mu,ier)

               if (ier.eq.1) then 
                  do i = 1, hcp
                     mu(i) = nopt(7)
                  end do
               end if 
             
            end if 

         end if 

      end if 
c                                 test for solvi and average
c                                 homogeneous phases.
      call avrger
      
1000  format (/,'**error ver901** solutions not allowed in saturated ',
     *   'component composition space',/,'in adaptive optimization ',
     *   'calculations, this limitation can be removed upon request.')

      end

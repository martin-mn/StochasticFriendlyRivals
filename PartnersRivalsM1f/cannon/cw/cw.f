* M1 strategies, CONTINUOUS space; WF pairwise comparison (synchronous
* Fermi imitation)
*
* [cw hdr] cw = DiskM1WF/Cannon/mw/sf1.f (byte copy: ref_sf1.f, sha256 in
* README.md) with the strategy space widened from the 16 binary memory-one
* strategies to the whole cube [eps,1-eps]^4 of EFFECTIVE memory-one
* strategies.  Everything not tagged [cw n] below is sf1.f verbatim;
* `diff ref_sf1.f cw.f` is the audit and CHANGES.md tabulates every hunk.
*   [cw 1] strategies live in LINEAGE slots (m1m2.F's class table)
*   [cw 2] payoffs from an ordered-pair cache keyed by lineage serials;
*          a miss solves the 4-state pair chain by GTH (pay4 below)
*   [cw 3] mutation / initial population: the exact-alpha biased kernel
*          (m1m2.F umut4, alpha = 1/3, ported to sf1's rb/irb buffer)
*   [cw 4] sampling: lineage abundance iab(ic) += icnt(ic) per sampled
*          generation (int64); the per-individual ef loop is sf1's
*   [cw 5] classification of each lineage's abundance at death and at the
*          end of a replicate: atom (E,N,R), stratum, corner (m1cls, C)
*   [cw 6] outputs: <isl>, w<isl>, h<isl> as described in README.md
*   [cw 7] -DBINMODE: sf1's 16-corner process THROUGH the lineage/cache
*          machinery; bit-identical to sf1.x (the audit of 1, 2, 4)
*   [cw 8] -DITEND=<n>: itend override for tests (default sf1's 1e7)
*   [cw 9] -DCLSREP=<k>: timing build only, m1cls called k times per
*          classification (measures the classifier's share of the time)
*
* A memory-one strategy answers the LAST ROUND'S OUTCOME, own action first:
* CC, CD, DC, DD.  Bit k of a corner index 0..15 is 1 for C (sf1's
* convention): ALLD = 0, Grim = 1, TFT = 5, WSLS = 9, ALLC = 15.  In the
* simulation a strategy is the effective vector s in [eps,1-eps]^4; the 16
* corners {eps,1-eps}^4 are exactly sf1's 16 binary strategies.
*
* The process, the games, the task decode, the RNG stream discipline, the
* Fermi update and the pay/ef accounting are sf1's; see ref_sf1.f's header
* for their history.  Continuum seeds 40000000 + (isl-1)*nper + irep;
* BINMODE keeps sf1's 10000000 so that it can be compared with sf1.x.
*
* usage:  ./cw.x <task index 1..ngam*neps>      (games.dat in the cwd)
* output (continuum):
*   unit 3 -> <task>  : nrep lines  gid ie u v Emax pay ef meanlive
*                       solves_per_gen  share(atom 0..7)
*   unit 4 -> w<task> : per-replicate diagnostics (see its header)
*   unit 7 -> h<task> : cell header; at the end '# atoms', '# strata',
*                       '# corners' sections of 'index count' lines and
*                       the closing '# total' line
* output (BINMODE): sf1.x's three files exactly; diagnostics on stdout.
C     [cw 8] itend override for short tests
#ifndef ITEND
#define ITEND 10000000
#endif
      program sb

      implicit real*8 (a-h,o-z)
      character*12 slurm,wfile,hfile
      character*200 lbuf,lb

      parameter (o0=0.d+00,oh=0.5d+00,o1=1.d+00,o2=2.d+00)
      parameter (n=1000,m=1)
      parameter (rn=dfloat(n),rrn=o1/rn,rrn1=o1/dfloat(n-1))
      parameter (rrm=oh/dfloat(m))
C     mgam is the game-list capacity ONLY -- it dimensions gcr/gcs/gct/gu/gv/
C     ggid/ggic/gtag and appears nowhere in the dynamics, the PRNG stream or
C     the payoff table, so changing it cannot move a number.  games.dat holds
C     the 512 sunflower games; 1024 is headroom.
      parameter (mgam=1024,neps=3,nu=1,nbet=1,nrep=10,nper=10)
C     [cw 7] seed base: sf1's in BINMODE, a disjoint one in the continuum
#ifdef BINMODE
      parameter (iseedb=10000000)
#else
      parameter (iseedb=40000000)
#endif
      dimension epsv(neps),uval(nu),bval(nbet)
      integer itv(nu)
C     the game list, read from games.dat at start-up
      dimension gcr(mgam),gcs(mgam),gct(mgam),gu(mgam),gv(mgam)
      integer ggid(mgam),ggic(mgam)
      character*24 gtag(mgam)
      parameter (ntpc=nrep/nper)
      parameter (nwin=10)
      parameter (imut=1)
      parameter (nst=16)
C     the 16 binary memory-one strategies, embedded in memory two (see the
C     header): mid(k) = 15 * sum_j bit_j(k) 16^j.  mid is filled at start-up
C     and is the ONLY place the memory-two world is mentioned.
      parameter (nrb=8192)
      parameter (thr=25.d+00)
C     [cw 1] lineage-table capacity.  nhi <= n always (a slot above the old
C     nhi is opened only when every free slot is used), so 2n+2 is headroom
      parameter (maxc=2*n+2)
C     [cw 2] pair-cache key = gid(ia)*ktbig + gid(ib), unique while gid < 2^31
      integer*8 ktbig
      parameter (ktbig=2147483648_8)

      dimension f(n),ps(4)
      data epsv /1.d-02,1.d-03,1.d-04/
      data uval /1.d-02/
      data bval /100.d+00/
C     run length per u; raise the 1e-4 entry to keep u*itend fixed
C     [cw 8] the literal 10000000 is now the macro ITEND (default 10000000)
      data itv  /ITEND/
C     [cw 1] s(n),sn(n) -> icls(n),iclsn(n): an individual holds a SLOT id
      integer iwin(nwin),mid(0:nst-1)
      integer*8 hist(0:nst-1),nsetl(4),ntot,ndist
      integer*8 htot(0:nst-1),ntotc
      integer*8 hbest
      dimension vlt(0:nst-1,0:nst-1),vht(0:nst-1,0:nst-1)
      dimension rb(nrb)
      character*4 bits

C     [cw 1] the lineage (class) table, as m1m2.F: slot ic holds the
C     strategy strat(:,ic), a task-monotone serial gid(ic) and the count
C     icnt(ic); ifree lists dead slots; mutants are staged in istg/sstr
C     and placed after the recount.  kid/kstg: BINMODE's corner index.
      integer icls(n),iclsn(n)
      integer icnt(maxc),ifree(maxc),istg(n)
      integer*8 gid(maxc),gcount,key
      dimension strat(4,maxc),sstr(4,n)
      integer kid(maxc),kstg(n)
C     [cw 2] ordered-pair payoff cache (m1m2.F's pkey/pc1/pc2)
      integer*8 pkey(maxc,maxc)
      dimension pc1(maxc,maxc),pc2(maxc,maxc)
C     [cw 4] lineage abundance over the sampled half (the spec's ab(ic))
      integer*8 iab(maxc)
C     [cw 5] per-replicate and per-task class counts, the corner cache and
C     the classification counters ncnt = (lineages classified, m1cls
C     calls, near ties, near-tie abundance, exact ties, exact-tie
C     abundance); near tie = m1cls neartie 1, exact tie = neartie 2
      integer*8 atomcnt(0:7),strcnt(0:80),corncnt(0:15),ncnt(6)
      integer*8 atot(0:7),strtot(0:80),corntot(0:15),nexp,nsum
      integer*8 nsolve,nlook
      integer icatom(0:15),icnear(0:15)
      dimension sv(4),share(0:7)

      common /epsval/ pmin,pmax
      common /phist/ hist
      common /thist/ htot
C     [cw 1,2] the big tables live in COMMON, as in m1m2.F
      common /bigcache/ pkey,pc1,pc2
      common /bigtab/ strat,sstr,gid,iab,icls,iclsn,icnt,ifree,istg,
     &                kid,kstg

      call get_command_argument(1,slurm)
      slurm=trim(adjustl(slurm))
      read(slurm,*) isl
      wfile='w'//slurm
      hfile='h'//slurm

C     ---- the game list ----------------------------------------------
C     one game per line:  gid  ic  cR  cS  cT  u  v  tag
C     cR/cS/cT are 17 significant digits and round-trip exactly; the u/v
C     columns are display only, so the dynamics never uses them.
C     ic is metadata (the donation cost index, 0 if not a donation game).
C     NOTE: DiskM2IN's games.dat is NOT interchangeable -- there ic>0 meant
C     'recompute the payoffs internally' and the cR/cS/cT columns of those
C     rows are zero placeholders.  Use Opt/mkgames.py from THIS project.
      open (unit=2,file='games.dat',status='old')
      ngam=0
   11 read (2,'(a)',end=12) lbuf
      lb=adjustl(lbuf)
      if (lb(1:1).eq.'#') goto 11
      if (len_trim(lb).eq.0) goto 11
      ngam=ngam+1
      if (ngam.gt.mgam) stop 'games.dat: too many games, raise mgam'
      read (lb,*) ggid(ngam),ggic(ngam),gcr(ngam),gcs(ngam),
     &            gct(ngam),gu(ngam),gv(ngam),gtag(ngam)
      goto 11
   12 close (2)
      if (ngam.le.0) stop 'games.dat: no games'

      open (unit=3,file=slurm)
      open (unit=4,file=wfile)
      open (unit=7,file=hfile)

C     which (beta, eps, u, c) cell this task carries, and which block of
C     replicates.  ntpc = nrep/nper tasks per cell; with nper = nrep the
C     task IS the cell and the abundance sum below is the cell's total,
C     otherwise the h-files of a cell's ntpc tasks are summed downstream
C     (Figures/wfabund.py does this).  nper = 1 is the right layout when a
C     single replicate is long enough to want its own wall clock.
C
C     task order, fastest first: replicate block, eps, u, beta, GAME.  The
C     game is outermost so that appending games to games.dat appends tasks
C     at the end and never moves an existing game's index.

      iblk =mod(isl-1,ntpc)+1
      ie   =mod((isl-1)/ntpc,neps)+1
      iu   =mod((isl-1)/(ntpc*neps),nu)+1
      ib   =mod((isl-1)/(ntpc*neps*nu),nbet)+1
      ig   =(isl-1)/(ntpc*neps*nu*nbet)+1
      if (ig.gt.ngam) stop 'task index beyond the game list'
      ic=ggic(ig)
      eps=epsv(ie)
      u=uval(iu)
      beta=bval(ib)
      itend=itv(iu)
      init=itend/2
      be=beta*rrm

C     pmin/pmax live in COMMON /epsval/ (shared with payf3), so they are
C     set here rather than declared as parameters
      pmin=eps
      pmax=o1-eps
      call payf2ini

C     THIS game's payoffs, (R,S,T,P) = (cR,cS,cT,0).  rr/ss/tt keep their
C     wf3 names so that the cache-miss line in the generation loop,
C     vl = rr*w0+ss*w1+tt*w2, is untouched -- with a donation row they are
C     (1-c,-c,1) and the arithmetic is bit-identical, not merely equal.

      rr=gcr(ig)
      ss=gcs(ig)
      tt=gct(ig)

C     the attainable optimum in self play: mutual R below the switch line
C     u+v = 1, alternation (S+T)/2 above it.  Efficiency is pay/emax.

      emax=dmax1(rr,oh*(ss+tt))

C     [cw 5] the game's (u,v) for the classifier, from the EXACT columns:
C     (R,S,T) = (1,u,1+v).  games.dat is on the 1/1024 grid, so these
C     equal the display columns; checked, not assumed.
      gamu=ss
      gamv=tt-o1
      if (rr.ne.o1) stop 'games.dat: cR is not 1'
      if (gamu.ne.gu(ig).or.gamv.ne.gv(ig))
     &   write (*,*) 'warning: games.dat u/v columns differ from cS/cT'

#ifdef BINMODE
      write (4,'(a)')
     & '# run: isl ie eps gid u v irun pay ef ndist ntot'
      write (4,'(a)')
     & '# win: isl eps gid irun rank code abund selfpay'//
     & ' selfef cc cd dc dd bits'
#else
C     [cw 6] the continuum w-file is per-replicate diagnostics
      write (4,'(a)')
     & '# diag: isl irun ncls ncall nneartie neartie_abund nhw'//
     & ' cache_hit_rate meanlive solves_per_gen gcount nexacttie'//
     & ' exacttie_abund'
      write (4,'(a)')
     & '#  ncls = lineages classified, ncall = m1cls calls (non-corner)'
     & //', nhw = class-table high water (maxc 2002); neartie = m1cls'
     & //' flag 1 (0<|margin|<=1e-9 scale), exacttie = flag 2 (margin'
     & //' exactly 0, a structural tie)'
#endif
      write (7,'(a,i7,a,i6,1x,a,a,e11.3,a,e11.3,a,f9.1)')
     & '# cell: isl ',isl,'  gid ',ggid(ig),trim(gtag(ig)),
     & '  eps ',eps,'  umut ',u,'  beta ',beta
      write (7,'(a,3f24.16,a,2f16.8,a,f24.16)')
     & '# game: cR cS cT ',rr,ss,tt,'   u v ',gu(ig),gv(ig),
     & '   Emax ',emax
      write (7,'(a,i5,a,i12,a,i3)')
     & '# run:  nper ',nper,'  itend ',itend,'  ie ',ie
      write (7,'(a,i4)') '# ic:   ',ic
C     [cw 6] sf1's histogram header line: BINMODE only
#ifdef BINMODE
      write (7,'(a)') '# code  count  (summed over nrep; code 0..15)'
#endif
      do is=0,nst-1
         htot(is)=0
      enddo
      ntotc=0

C     the game's payoffs are set above, before the headers that name them

C     the 16 x 16 payoff table, built once per task with dw1's own payf3 on
C     the embedded ids and dw1's own two payoff lines.  vlt(i,j) is i's
C     per-round payoff against j and vht(i,j) is j's against i, so the pair
C     (i,j) needs one call and the generation loop needs no min/max swap.
C     With 16 strategies the table is 136 distinct calls against dw1's hash
C     cache of 2^20 slots; the loop reads it in place of the cache and draws
C     the same random numbers in the same order.
C     [cw 7] mid is kept (the w-file's payf3 call); the table is built in
C     BINMODE only, where it fills the pair cache on a miss

      do k=0,nst-1
         mid(k)=0
         do j=0,3
            if (btest(k,j)) mid(k)=mid(k)+15*(16**j)
         enddo
      enddo
#ifdef BINMODE
      do ia=0,nst-1
         do ib=ia,nst-1
            call payf3(mid(ia),mid(ib),w0,w1,w2,w3)
            vl=rr*w0+ss*w1+tt*w2
            vh=rr*w0+tt*w1+ss*w2
            vlt(ia,ib)=vl
            vht(ia,ib)=vh
            vlt(ib,ia)=vh
            vht(ib,ia)=vl
         enddo
      enddo
#endif

C     [cw 5] the corner cache: the 16 corners classified once per task.
C     The stratum the classifier returns is checked against the index.
      do k=0,nst-1
         jst=0
         do j=1,4
            if (btest(k,j-1)) then
               sv(j)=pmax
               jst=jst+3**(j-1)
            else
               sv(j)=pmin
            endif
         enddo
         call cwcls1(sv,eps,gamu,gamv,icatom(k),kst,icnear(k))
         if (kst.ne.jst) stop 'm1cls: corner stratum mismatch'
      enddo
      do k=0,7
         atot(k)=0
      enddo
      do k=0,80
         strtot(k)=0
      enddo
      do k=0,15
         corntot(k)=0
      enddo
#ifndef BINMODE
      write (7,'(a,16i2)') '# corner atoms (code 0..15): ',
     &      (icatom(k),k=0,15)
#endif

C     [cw 2] empty the pair cache once per task.  gid is monotone over the
C     whole task (never reset between replicates), so no stale key can
C     ever match and the cache needs no per-replicate clearing.
      do jc=1,maxc
         do kc=1,maxc
            pkey(kc,jc)=0
         enddo
      enddo
      gcount=0

C     ================= replicate loop =================

      do irep=1,nper

      irun=(iblk-1)*nper+irep
      iseed=iseedb+(isl-1)*nper+irep
      call init_genrand(iseed)
      irb=nrb

      do is=0,nst-1
         hist(is)=0
      enddo
      z=o0
      ef=o0
C     [cw 1,4,5] per-replicate lineage and class state
      do is=1,maxc
         icnt(is)=0
         iab(is)=0
      enddo
      do k=0,7
         atomcnt(k)=0
      enddo
      do k=0,80
         strcnt(k)=0
      enddo
      do k=0,15
         corncnt(k)=0
      enddo
      do k=1,6
         ncnt(k)=0
      enddo
      sdiv=o0
      nsolve=0
      nlook=0

C     initial population: 4 fair coins per player -- the memory-one analogue
C     of dw1's 16, i.e. uniform on the 16 strategies
C     [cw 3] continuum: n independent draws of the exact-alpha kernel
C     [cw 1] one lineage slot per individual, as m1m2.F

      do i=1,n
#ifdef BINMODE
C        [cw 7] sf1's 4 fair coins; the slot records the corner index kid
C        and its effective vector (pmin/pmax)
         is=0
         do k=1,4
            if (irb.ge.nrb) then
               call mtfill(rb,nrb)
               irb=0
            endif
            irb=irb+1
            if (rb(irb).ge.oh) is=ibset(is,k-1)
         enddo
         kid(i)=is
         do k=1,4
            strat(k,i)=pmin
            if (btest(is,k-1)) strat(k,i)=pmax
         enddo
#else
         call umut4(rb,irb,strat(1,i))
#endif
         gcount=gcount+1
         gid(i)=gcount
         icls(i)=i
         icnt(i)=1
      enddo
      nhi=n
      nlive=n
      nhw=n

C     dynamics

      do it=1,itend

C        fitnesses: play against m random others

         do i=1,n
            f(i)=o0
         enddo

         do i=1,n
C           [cw 1] si=s(i) -> the slot ia
            ia=icls(i)
            do l=1,m
   43          if (irb.ge.nrb) then
                  call mtfill(rb,nrb)
                  irb=0
               endif
               irb=irb+1
               j=1+int(rb(irb)*rn)
               if (j.eq.i) goto 43
C              [cw 2] sj=s(j) -> the slot ib; payoffs via the pair cache
               ib=icls(j)
               key=gid(ia)*ktbig+gid(ib)
               if (pkey(ia,ib).ne.key) then
#ifdef BINMODE
C                 [cw 7] a miss reads sf1's own table entry
                  f1=vlt(kid(ia),kid(ib))
                  f2=vht(kid(ia),kid(ib))
#else
C                 [cw 2] a miss solves the 4-state pair chain (GTH) and
C                 forms sf1's own two payoff lines
                  call pay4(strat(1,ia),strat(1,ib),w0,w1,w2,w3)
                  f1=rr*w0+ss*w1+tt*w2
                  f2=rr*w0+tt*w1+ss*w2
#endif
                  pc1(ia,ib)=f1
                  pc2(ia,ib)=f2
                  pkey(ia,ib)=key
                  nsolve=nsolve+1
               endif

C              [cw 2] vlt(si,sj)/vht(si,sj) -> the cached pc1/pc2
               f(i)=f(i)+pc1(ia,ib)
               f(j)=f(j)+pc2(ia,ib)
            enddo
         enddo

C        sample the second half

         if (it.gt.init) then
            z=z+o1
            efg=o0
            do i=1,n
               efg=efg+f(i)
            enddo
            ef=ef+efg
C           [cw 4] hist(s(i))+=1 per individual -> iab(ic)+=icnt(ic) per
C           lineage; hist is rebuilt from iab when the lineage is
C           classified [cw 5]
            do is=1,nhi
               if (icnt(is).gt.0) iab(is)=iab(is)+icnt(is)
            enddo
            sdiv=sdiv+dfloat(nlive)
         endif

C        update: mutation or pairwise comparison
C        [cw 1] a mutant is STAGED (iclsn(i) = -its stage index) and gets
C        its slot after the recount below

         nstg=0

         do i=1,n

            if (irb.ge.nrb) then
               call mtfill(rb,nrb)
               irb=0
            endif
            irb=irb+1

            if (rb(irb).lt.u) then

               nstg=nstg+1
               istg(nstg)=i
               iclsn(i)=-nstg
#ifdef BINMODE
C              [cw 7] sf1's imut = 1 path: one integer draw supplies all
C              4 loci.  (sf1's dead imut /= 1 branch is not carried.)
               if (irb.ge.nrb) then
                  call mtfill(rb,nrb)
                  irb=0
               endif
               irb=irb+1
               kstg(nstg)=int(rb(irb)*16.d+00)
#else
C              [cw 3] one draw of the exact-alpha kernel, 4 uniforms
               call umut4(rb,irb,sstr(1,nstg))
#endif

            else

               if (irb.ge.nrb) then
                  call mtfill(rb,nrb)
                  irb=0
               endif
               irb=irb+1
               j=1+int(rb(irb)*rn)

C              logistic imitation probability, xf = be*df with
C              be = beta/2m and df a sum over 2m games, so the intensity
C              acting on a per-round payoff difference is beta.
C              |xf|>thr is decided
C              without dexp: genrand_real3 lies in
C              [2**-33, 1-2**-33], so for thr=25 the comparison
C              ran<prob has the same outcome as with prob=1 or 0.

               xf=be*(f(j)-f(i))
               if (xf.ge.thr) then
                  prob=o1
               else if (xf.le.-thr) then
                  prob=o0
               else if (xf.ge.o0) then
                  prob=o1/(o1+dexp(-xf))
               else
                  ex=dexp(xf)
                  prob=ex/(o1+ex)
               endif

               if (irb.ge.nrb) then
                  call mtfill(rb,nrb)
                  irb=0
               endif
               irb=irb+1
C              [cw 1] sn(i)=s(j) / s(i) -> the slot of j / of i
               if (rb(irb).lt.prob) then
                  iclsn(i)=icls(j)
               else
                  iclsn(i)=icls(i)
               endif

            endif

         enddo

C        [cw 1] new generation (replaces s(i)=sn(i)): reference-count,
C        recycle dead slots, place the staged mutants -- m1m2.F's step.
C        [cw 5] a slot that dies with iab > 0 is classified here, before
C        it can be reused.

         do is=1,nhi
            icnt(is)=0
         enddo
         do i=1,n
            is=iclsn(i)
            if (is.gt.0) icnt(is)=icnt(is)+1
         enddo

         nfree=0
         nh=0
         nkeep=0
         do is=1,nhi
            if (icnt(is).eq.0) then
               nfree=nfree+1
               ifree(nfree)=is
               if (iab(is).gt.0) then
#ifdef BINMODE
                  hist(kid(is))=hist(kid(is))+iab(is)
#endif
                  call clsab(strat(1,is),iab(is),eps,gamu,gamv,
     &                 icatom,icnear,atomcnt,strcnt,corncnt,ncnt)
                  iab(is)=0
               endif
            else
               nh=is
               nkeep=nkeep+1
            endif
         enddo

         kf=0
         isnext=nhi+1
         do is=1,nstg
            if (kf.lt.nfree) then
               kf=kf+1
               jc=ifree(kf)
            else
               jc=isnext
               isnext=isnext+1
            endif
            if (jc.gt.maxc) stop 'class table overflow'
            gcount=gcount+1
            gid(jc)=gcount
#ifdef BINMODE
            kid(jc)=kstg(is)
            do k=1,4
               strat(k,jc)=pmin
               if (btest(kstg(is),k-1)) strat(k,jc)=pmax
            enddo
#else
            do k=1,4
               strat(k,jc)=sstr(k,is)
            enddo
#endif
            iclsn(istg(is))=jc
            icnt(jc)=1
            if (iab(jc).ne.0) stop 'recycled slot with abundance'
            if (jc.gt.nh) nh=jc
         enddo
         nhi=nh
         nlive=nkeep+nstg
         if (nhi.gt.nhw) nhw=nhi

         do i=1,n
            icls(i)=iclsn(i)
         enddo

      enddo

C     [cw 5] classify every lineage still holding abundance (the live ones)
C     and assert the invariant sum(atomcnt) == (itend-init)*n exactly
      nlook=int(n,8)*int(m,8)*int(itend,8)
      do is=1,nhw
         if (iab(is).gt.0) then
#ifdef BINMODE
            hist(kid(is))=hist(kid(is))+iab(is)
#endif
            call clsab(strat(1,is),iab(is),eps,gamu,gamv,
     &           icatom,icnear,atomcnt,strcnt,corncnt,ncnt)
            iab(is)=0
         endif
      enddo
      if (gcount.ge.ktbig) stop 'gid overflow: pair-cache key unsafe'
      nexp=int(itend-init,8)*int(n,8)
      nsum=0
      do k=0,7
         nsum=nsum+atomcnt(k)
      enddo
      if (nsum.ne.nexp) then
         write (*,*) 'INVARIANT FAILED: isl irep sum expected ',
     &        isl,irep,nsum,nexp
         stop 'invariant sum(atomcnt) == (itend/2)*n violated'
      endif
      nsum=0
      do k=0,80
         nsum=nsum+strcnt(k)
      enddo
      if (nsum.ne.nexp) stop 'invariant on strcnt violated'

C     ---- output ------------------------------------------------------

C     pay = mean per-round payoff of a player, ef = pay/emax.  wf3 divided
C     by (1-c), which IS emax on the donation ray below the switch line.
      pay=ef*rrn*rrm/z
      ef=pay/emax

C     [cw 6,7] sf1's output block, verbatim, is BINMODE's; the continuum
C     writes its own summary line and diagnostics after the #else
#ifdef BINMODE
      ntot=0
      ndist=0
      do k=1,4
         nsetl(k)=0
      enddo
      do is=0,nst-1
         if (hist(is).ne.0) then
            ndist=ndist+1
            ntot=ntot+hist(is)
            htot(is)=htot(is)+hist(is)
            do k=1,4
               if (btest(is,k-1)) nsetl(k)=nsetl(k)+hist(is)
            enddo
         endif
      enddo
      do k=1,4
         ps(k)=(pmax*dble(nsetl(k))
     &         +pmin*dble(ntot-nsetl(k)))*rrn/z
      enddo

      ntotc=ntotc+ntot
      write (3,'(2i8,2f14.6,7f16.10)') ggid(ig),ie,gu(ig),gv(ig),
     &      emax,pay,ef,(ps(k),k=1,4)
      call flush(3)

C     the nwin most abundant strategies of the second half

      do iw=1,nwin
         iwin(iw)=-1
         hbest=-1
         do is=0,nst-1
            if (hist(is).gt.hbest) then
               idup=0
               do jw=1,iw-1
                  if (iwin(jw).eq.is) idup=1
               enddo
               if (idup.eq.0) then
                  hbest=hist(is)
                  iwin(iw)=is
               endif
            endif
         enddo
      enddo

      write (4,'(a,2i7,e11.3,i7,2f10.4,i7,2f12.7,i9,i16)') '# ',
     &      isl,ie,eps,ggid(ig),gu(ig),gv(ig),irun,pay,ef,
     &      ndist,ntot
      do iw=1,nwin
         is=iwin(iw)
         ab=dble(hist(is))/dble(ntot)
         call payf3(mid(is),mid(is),w0,w1,w2,w3)
         spay=rr*w0+ss*w1+tt*w2
         do k=1,4
            if (btest(is,k-1)) then
               bits(k:k)='C'
            else
               bits(k:k)='D'
            endif
         enddo
         write (4,'(i8,e11.3,i7,2i6,i7,f12.8,2f12.7,4f11.8,1x,a4)')
     &         isl,eps,ggid(ig),irun,iw,is,ab,spay,spay/emax,
     &         w0,w1,w2,w3,bits
      enddo
      call flush(4)
C     [cw 7] BINMODE diagnostics go to stdout, keeping sf1's files exact
      write (*,'(a,2i6,a,i12,a,6i12,a,f10.4,a,f12.8)')
     &   'BIN isl irep ',isl,irep,'  sum(atomcnt) ',nexp,
     &   '  ncls ncall nnear nearab ntie tieab ',(ncnt(k),k=1,6),
     &   '  meanlive ',sdiv/z,'  solves/gen ',
     &   dble(nsolve)/dble(itend)
#else
C     [cw 6] the continuum summary line and diagnostics
      ntot=nexp
      ntotc=ntotc+ntot
      do k=0,7
         share(k)=dble(atomcnt(k))/dble(ntot)
      enddo
      write (3,'(2i8,2f14.6,3f16.10,2f12.4,8f16.12)') ggid(ig),ie,
     &      gu(ig),gv(ig),emax,pay,ef,sdiv/z,
     &      dble(nsolve)/dble(itend),(share(k),k=0,7)
      call flush(3)
      write (4,'(a,i7,i4,2i12,i10,i16,i6,f12.8,f10.3,f10.4,i12,
     &          i12,i16)')
     &      '# ',isl,irun,ncnt(1),ncnt(2),ncnt(3),ncnt(4),nhw,
     &      o1-dble(nsolve)/dble(nlook),sdiv/z,
     &      dble(nsolve)/dble(itend),gcount,ncnt(5),ncnt(6)
      call flush(4)
#endif
C     [cw 5] task totals
      do k=0,7
         atot(k)=atot(k)+atomcnt(k)
      enddo
      do k=0,80
         strtot(k)=strtot(k)+strcnt(k)
      enddo
      do k=0,15
         corntot(k)=corntot(k)+corncnt(k)
      enddo

      enddo

C     ================= end replicate loop =============

C     the summed abundance of every strategy over this task's nper
C     replicates.  ntotc = nper*(itend/2)*n is the grand total, so
C     count/ntotc is the mean frequency of that strategy in the sampled
C     half of a run.  Sum the h-files of a cell's ntpc tasks for the cell.

C     [cw 6,7] sf1's histogram tail in BINMODE; the continuum tail below
#ifdef BINMODE
      ndist=0
      do is=0,nst-1
         if (htot(is).ne.0) ndist=ndist+1
      enddo
      write (7,'(a,i8,a,i16)') '# distinct ',ndist,'   total ',ntotc
      do is=0,nst-1
         if (htot(is).ne.0) write (7,'(i7,1x,i16)') is,htot(is)
      enddo
#else
C     [cw 6] the continuum h-file tail, in one burst: atoms (index =
C     4E+2N+R), strata (sum code_k 3^k, code 0: p_k=0, 1: p_k=1, 2:
C     interior), corners (sf1's bit convention), then the closing total.
      write (7,'(a)') '# atoms 8  (index 4E+2N+R, count)'
      do k=0,7
         write (7,'(i7,1x,i16)') k,atot(k)
      enddo
      write (7,'(a)') '# strata 81  (index sum code_k 3^k, count)'
      do k=0,80
         write (7,'(i7,1x,i16)') k,strtot(k)
      enddo
      write (7,'(a)') '# corners 16  (sf1 code 0..15, count)'
      do k=0,15
         write (7,'(i7,1x,i16)') k,corntot(k)
      enddo
      write (7,'(a,i16)') '# total ',ntotc
#endif
      call flush(7)

      stop
      end

**********************************************************************
*     [cw 3] the exact-alpha biased kernel: m1m2.F's umut4 (-DALPHA path)
*     ported to sf1's buffer (rb, 1-based irb, mtfill).  Each coordinate
*     is independently exactly pmin w.p. alpha, exactly pmax w.p. alpha,
*     and Uniform(pmin,pmax) otherwise, alpha = 1/3, from ONE uniform:
*       uw = (1-2 pmin)/(1-2 alpha),  ua = pmin - alpha*uw,
*       ran = ua + rb*uw  clipped to [pmin, pmax].
*     This deliberately does NOT reproduce the original s.f literal
*     ua = -1, uw = 3 (alpha = (1+eps)/3).  Returns the 4 memory-one
*     coordinates (CC, CD, DC, DD), not m1m2's 16-vector embedding.
**********************************************************************

      subroutine umut4(rb,irb,sv)
      implicit real*8 (a-h,o-z)
      parameter (nrb=8192)
      dimension rb(nrb),sv(4)
      common /epsval/ pmin,pmax
      parameter (ualph=1.d+00/3.d+00)

      uw=(1.d+00-2.d+00*pmin)/(1.d+00-2.d+00*ualph)
      ua=pmin-ualph*uw
      do k=1,4
         if (irb.ge.nrb) then
            call mtfill(rb,nrb)
            irb=0
         endif
         irb=irb+1
         ran=ua+rb(irb)*uw
         if (ran.lt.pmin) then
            sv(k)=pmin
         else if (ran.gt.pmax) then
            sv(k)=pmax
         else
            sv(k)=ran
         endif
      enddo
      return
      end

**********************************************************************
*     [cw 2] pay4: stationary distribution of the FOUR-state pair chain
*     of memory-one strategies sa (focal) and sb, by GTH state reduction
*     (subtraction-free; the same algorithm and storage as payf3).
*     States from sa's side: 1 CC, 2 CD, 3 DC, 4 DD.  In state i the
*     focal player cooperates w.p. sa(i); the co-player sees the SWAPPED
*     state (focal's CD is its DC) and cooperates w.p. sb(isw(i)).
*     w0..w3 = stationary P(CC), P(CD), P(DC), P(DD) from sa's side.
**********************************************************************

      subroutine pay4(sa,sb,w0,w1,w2,w3)
      implicit none
      double precision sa(4),sb(4),w0,w1,w2,w3
      double precision g(4,4),v(4),p,q,s,sinv,pin,t
      integer i,j,k,isw(4)
      data isw /1,3,2,4/

C---- g(j,i) = P(i -> j), destination first (payf3's layout)
      do i=1,4
         p=sa(i)
         q=sb(isw(i))
         g(1,i)=p*q
         g(2,i)=p*(1.d0-q)
         g(3,i)=(1.d0-p)*q
         g(4,i)=(1.d0-p)*(1.d0-q)
      enddo

C---- GTH state reduction
      do k=4,2,-1
         s=0.d0
         do j=1,k-1
            s=s+g(j,k)
         enddo
         sinv=1.d0/s
         do i=1,k-1
            pin=g(k,i)
            if (pin.ne.0.d0) then
               pin=pin*sinv
               g(k,i)=pin
               do j=1,k-1
                  g(j,i)=g(j,i)+pin*g(j,k)
               enddo
            endif
         enddo
      enddo

C---- back substitution
      v(1)=1.d0
      t=1.d0
      do k=2,4
         s=0.d0
         do i=1,k-1
            s=s+v(i)*g(k,i)
         enddo
         v(k)=s
         t=t+s
      enddo
      t=1.d0/t
      w0=v(1)*t
      w1=v(2)*t
      w2=v(3)*t
      w3=v(4)*t
      return
      end

**********************************************************************
*     [cw 5] clsab: classify one lineage and add its abundance nab.
*     The stratum is computed here with the classifier's exact rule
*     (s == pmin -> code 0, s == 1-pmin -> code 1, else 2); a corner is
*     looked up in the per-task corner cache, anything else calls m1cls
*     and the stratum it returns is checked.
*     ncnt = (lineages classified, m1cls calls, near ties, near-tie
*     abundance, exact ties, exact-tie abundance).  m1cls returns
*     neartie 1 for a rounding-sensitive Nash margin (0 < |m| <= 1e-9
*     scale) and 2 for a margin that is exactly 0.0 (a structural tie,
*     e.g. pi(ALLC,p) = pi(p,p) on the slab p1 = 1): counted apart.
**********************************************************************

      subroutine clsab(sv,nab,eps,gamu,gamv,icatom,icnear,
     &                 atomcnt,strcnt,corncnt,ncnt)
      implicit none
      double precision sv(4),eps,gamu,gamv,pmx
      integer*8 nab,atomcnt(0:7),strcnt(0:80),corncnt(0:15),ncnt(6)
      integer icatom(0:15),icnear(0:15)
      integer k,ist,icor,iat,inr,jst,mul,icode
      logical corner

      pmx=1.d0-eps
      ist=0
      icor=0
      mul=1
      corner=.true.
      do k=1,4
         if (sv(k).eq.eps) then
            icode=0
         else if (sv(k).eq.pmx) then
            icode=1
            icor=ibset(icor,k-1)
         else
            icode=2
            corner=.false.
         endif
         ist=ist+icode*mul
         mul=mul*3
      enddo
      if (corner) then
         iat=icatom(icor)
         inr=icnear(icor)
      else
         call cwcls1(sv,eps,gamu,gamv,iat,jst,inr)
         ncnt(2)=ncnt(2)+1
         if (jst.ne.ist) stop 'm1cls: stratum mismatch'
      endif
      if (iat.lt.0.or.iat.gt.7) stop 'm1cls: atom out of range'
      atomcnt(iat)=atomcnt(iat)+nab
      strcnt(ist)=strcnt(ist)+nab
      if (corner) corncnt(icor)=corncnt(icor)+nab
      ncnt(1)=ncnt(1)+1
      if (inr.eq.1) then
         ncnt(3)=ncnt(3)+1
         ncnt(4)=ncnt(4)+nab
      else if (inr.eq.2) then
         ncnt(5)=ncnt(5)+1
         ncnt(6)=ncnt(6)+nab
      else if (inr.ne.0) then
         stop 'm1cls: neartie flag out of range'
      endif
      return
      end

**********************************************************************
*     [cw 5] cwcls1: the ONE place the C classifier is called, through
*     an explicit bind(C) interface (no compiler-specific name mangling)
*       void m1cls(const double s[4], const double *eps, const double *u,
*                  const double *v, int *atom, int *stratum, int *neartie)
**********************************************************************

      subroutine cwcls1(sv,eps,gamu,gamv,iat,ist,inr)
      use iso_c_binding
      implicit none
      double precision sv(4),eps,gamu,gamv
      integer iat,ist,inr
      integer(c_int) ja,js,jn
#ifdef CLSREP
      integer irp
#endif
      interface
         subroutine m1cls(s,e,u,v,iatom,istr,inear)
     &              bind(C,name='m1cls')
         use iso_c_binding
         real(c_double) s(4),e,u,v
         integer(c_int) iatom,istr,inear
         end subroutine m1cls
      end interface

#ifdef CLSREP
C     [cw 9] timing build only: the call is repeated CLSREP-1 extra
C     times, so (t(CLSREP) - t(1))/(CLSREP-1) is the classifier's cost
      do irp=2,CLSREP
         call m1cls(sv,eps,gamu,gamv,ja,js,jn)
      enddo
#endif
      call m1cls(sv,eps,gamu,gamv,ja,js,jn)
      iat=ja
      ist=js
      inr=jn
      return
      end

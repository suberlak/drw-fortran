

	program main	

        parameter (NMAX=1600)
        implicit real*8 (a-h,o-z)
        character*100 fname
        common /object/fname

        real*8 jd(NMAX),mag(NMAX),emag(NMAX)
 
        real*8 tempmag(NMAX)

        real*8 jwant(10*NMAX),mwant(10*NMAX),ewant(10*NMAX)

        real*8 lpar(10),lmat(NMAX,10),lcovar(10,10)
        integer id(NMAX)

        real*8 bar0(10),bar1(10)
  
        real*8 tmag(10),temag(10),tjd(10)
        real*8 shift(10),scale(10)
        real*8 bshift(10),bscale(10)
 
        real*8 zin,lin,mbhin

        real*8 ln10

        real*8 om,ol,ok,ztab(100),dtab(100),ztabs(100),dtabs(100)
        integer ntab
        common /cosmo/om,ol,ok,ztab,dtab,ztabs,dtabs,ntab

        integer PR,FOR,NORMAL,NOISE,INFINITE
        common /fittype/PR,FOR,NORMAL,NOISE,INFINITE

        integer igrid
        common /flags/igrid

        common /limfits/pnoise,cnoise,snoise,pinf,cinf,sinf

c methods -- Press/Rybicki or forecasting
        PR = 0
        FOR= 1
c limits of covariance matrix in PR
c NOISE limit is tau->0, just diagonal like noise matrix
c INFINITE limit is tau->infinty, all elements same
        NORMAL   = 0
        NOISE    = 1
        INFINITE = 2

c set the cosmology
        om = 0.3
        ol = 0.7
        ok = 1.0-om-ol
        print*,'cosmology set to om= ',om
        print*,'                 ol= ',ol
        print*,'                 ok= ',ok
c inititalize the distances
        ntab = 100
        call distinit

        ln10 = log(10.0)

c it wants a file in the directory giving the filter number u=1 gri z=5
        open(unit=13,file='input.filter',form='formatted',status='old')
        read(13,*)ifilt
        close(unit=13)

c this gives a list of things to be run
        open(unit=9,file='process.dat',form='formatted',status='old')
c        open(unit=9,file='short.dat',form='formatted',status='old')
        read(9,*)nqso
        print*,'will try to do ',nqso,' qsos '
        do kk=1,nqso

        read(9,*)nline,rain,decin,zin,absmagi,fname

        call splint(ztab,dtab,dtabs,ntab,zin,dol)
        rh  = 3000.0
        dollum = rh*dol*(1.0+zin)
        dmod   = 5.0*log10(dollum*1.0D5)

	print*,' trying to do lightcurve file ',fname
        print*,' lines ',nline,' coord ',rain,decin,' z ',zin,' absmag '
     &      ,absmagi

        open(unit=13,file=fname,form='formatted',status='old')
c        read(13,*)nline
c        read(13,*)
c        read(13,*)
        ncurve = 1
c        if (nline*ncurve.gt.NMAX) then
        if (nline.gt.NMAX) then
          print*,' light curve too long for NMAX ',nline*ncurve,
     &        ' versus ',NMAX
          stop
          endif
        npt = 0
        do i=1,nline
          read(13,*)(tjd(k),tmag(k),temag(k),k=1,ncurve)
          if ((tmag(ifilt).gt.-30.0).and.(tmag(ifilt).le.30.0)) then
            igood = 1
            if (i.ne.1) then
              if (tjd(ifilt)-jd(npt).le.0.0001) then
                igood = 0
                print*,'bad point in ',fname,' at ',i,' jd ',tjd(ifilt)
                endif
              endif
            if (igood.eq.1) then
              npt = npt + 1
              jd(npt)   = tjd(ifilt)
              mag(npt)  = tmag(ifilt)
              emag(npt) = temag(ifilt)
              id(npt)   = 1
              endif
            endif
          enddo
        close(unit=13)
        print*,'got ',npt,' good points from entry ',kk,fname
c        do i=1,5
c          print*,jd(i),mag(i),emag(i),id(i)
c          enddo

c print out a likelihood plot
        igrid = 0
 
        print*,'resetting to one light curve after read '
        ncurve = 1

c we are using the fast method -- we cannot have duplicate
c epochs -- 
        ibad = 0
        do i=2,npt
          if (jd(i).le.jd(i-1)) then
            print*,'bad point ',jd(i-1),jd(i)
            estmag = mag(i-1)/emag(i-1)**2+mag(i)/emag(i)**2
            esterr =      1.0/emag(i-1)**2+   1.0/emag(i)**2
            estmag = estmag/esterr
            esterr = 1.0/sqrt(esterr)
            print*,'merge as ',jd(i),estmag,esterr
            ibad = 1
            endif
          enddo
        if (ibad.eq.1) then
          print*,'duplicated JD in ',fname,' at index ',i,' time ',jd(i)
          stop
          endif
c determine and subtract off the averages
        do j=1,ncurve
          bar0(j) = 0.0
          bar1(j) = 0.0
          enddo
        do i=1,npt
          j       = id(i) 
          bar0(j) = bar0(j) +    1.0/emag(i)**2
          bar1(j) = bar1(j) + mag(i)/emag(i)**2
          enddo
        do j=1,ncurve
          bar1(j) = bar1(j)/bar0(j)
          print*,'mean for curve ',j,' is ',bar1(j),' divided by ',
     &      bar0(j)
          enddo
        chicon = 0.0
        do i=1,npt
          mag(i) = mag(i)-bar1(id(i))
          chicon = chicon + (mag(i)/emag(i))**2
          enddo
c fit a line
        a1 = 0.0
        a2 = 0.0
        aa = 0.0
        ab = 0.0
        bb = 0.0
        do i=1,npt
          delt = jd(i)-jd(1)
          a1 = a1 + mag(i)/emag(i)**2
          a2 = a2 + delt*mag(i)/emag(i)**2
          aa = aa + 1.0/emag(i)**2
          ab = ab + delt/emag(i)**2
          bb = bb + delt*delt/emag(i)**2
          enddo
        det = aa*bb-bb*bb
        aai =  bb/det
        abi = -ab/det
        bbi =  aa/det
        b1  = a1*aai+abi*a2
        b2  = a1*abi+bbi*a2
        chilin = 0.0
        do i=1,npt
          tmodel = b1 + (jd(i)-jd(1))*b2
          chilin = chilin + ((mag(i)-tmodel)/emag(i))**2
          enddo

      

c warning -- seed is hard wired...
        iseed = -23424

        nlin    = 1
        nlindim = 10
        do i=1,npt
          lmat(i,1) = 1.0
          enddo
c do the PR fit
        call dofit(kk,jd,mag,emag,npt,bar1(1),lmat,lpar,lcovar,nlin,nlindim,PR ,40,
     1            plikea,chi2a,tbesta,tlowa,thiga,sbesta,slowa,shiga)
        avgmaga = bar1(1) + lpar(1)
        write(39,'(i5,1x,14(1x,g13.6),1x,i5,1x,a65)')kk,zin,dmod,absmagi,bar1(1),chicon,chilin,pnoise,pinf,plikea,cnoise,cinf,chi2a,snoise,sinf,npt,fname
c do the forecasting fit
        call dofit(0,jd,mag,emag,npt,bar1(1),lmat,lpar,lcovar,nlin,nlindim,FOR,41,
     1            plikeb,chi2b,tbestb,tlowb,thigb,sbestb,slowb,shigb)
c generate the Monte Carlo realization
        tau   = 10.0**tbesta
        sigma = 10.0**sbesta
        call gencurve(jd,tempmag,emag,npt,tau,sigma,avgmaga,iseed)
c rerun the two fits
        call dofit(0,jd,tempmag,emag,npt,bar1(1),lmat,lpar,lcovar,nlin,nlindim,PR ,50,
     1            plikec,chi2c,tbestc,tlowc,thigc,sbestc,slowc,shigc)
        call dofit(0,jd,tempmag,emag,npt,bar1(1),lmat,lpar,lcovar,nlin,nlindim,FOR,51,
     1            pliked,chi2d,tbestd,tlowd,thigd,sbestd,slowd,shigd)
c search for bad points by sequentially increasing errors on each point for
c fits where the chi^2 per dof exceeds some limit -- use the forecasting 
c method for this, since it should be faster - not as elegant as sequentially
c using PR to predict each point given the other fits, but ....
        print*,'fit has chi2a = ',chi2a
        if ((chi2a/float(npt-3).ge.2.0).and.(npt.ge.10)) then
          print*,'doing bad point search '
          do k=1,npt
            tempmag(k) = emag(k)
            enddo 
          dfitmax = 0.0
          do k=1,npt
            tempmag(k) = 1000.0
            call dofit(0,jd,mag,tempmag,npt,bar1(1),lmat,lpar,lcovar,nlin,nlindim,PR,91,
     1              pliket,chi2t,tbestt,tlowt,thigt,sbestt,slowt,shigt)
            tempmag(k) = emag(k)
            dfit       = chi2a - chi2t
c            write(81,'(10(1x,f13.6))')dfit,chi2t,chi2b,jd(k),mag(k),emag(k)
            if (dfit.ge.dfitmax) then
              cfitmax = chi2t
              dfitmax = dfit
              ifitmax = k
              endif 
            enddo
          write(80,'(i5,6(1x,f13.6),1x,a65)')npt,dfitmax,chi2a,cfitmax,jd(ifitmax),bar1(1)+mag(ifitmax),emag(ifitmax),fname
          endif
        enddo

      
        end


c build the inverse of an exponential covariance matrix in the format
c which will match the NR tridag solver, using the formulas from Rybicki
c check the inverse on the way out

c determines the log of the determinant while it is at it
c   determinant = Pi_i 1/(1-ri^2)
c so log(determinant) = - sum log(1-ri^2)

c a basic problem with this approach is that the matrices go singular as 
c the ri->1 (i.e. equal times, where times become equal as dt/tau->0)
c we can make the covariance matrix non-singular in one way by multiplying
c by the diagonal matrix with elements
c  sqrt(1-r1^2) ... sqrt(1-ri^2)sqrt(1-ri+1)^2 .... sqrt(1-rn-1^2)
c this has determinant Pi_i (1-ri^2) so the determinant of the
c   covariance matrix times this diagonal matrix is now unity
c note, however, that the elements of the matrix still diverge for any ri->1
c but more slowly -- we have to multiply the same diagonal matrix into the
c noise as well.

      subroutine buildcinv(atridag,btridag,ctridag,diag,detlog,jd,npt,sigma,tau)
        implicit real*8 (a-h,o-z)
        PARAMETER (NMAX=1600)
        real*8 atridag(*),btridag(*),ctridag(*),jd(*)
        real*8 diag(*)
        real*8 ri(NMAX)

c if isigma=1, put sigma into the diagonal matrix
        isigma = 1


        detlog = 0.0
      
        do i=1,npt-1
          if (jd(i+1).le.jd(i)) then
            print*,' this does not presently work if the data is not in time order '
            print*,' and will have a divergence for equal time data '
            print*,' offender is points ',i,i+1,jd(i),jd(i+1)
            stop
            endif
          arg  = abs(jd(i+1)-jd(i))/tau
          if (arg.ge.1.e-4) then
            ri(i)      = exp(-arg)
            ctridag(i) = -ri(i)/(1.0-ri(i)*ri(i))
            detlog     = detlog - log(1.0-ri(i)*ri(i))
          else
            arg2       = arg*arg
            ctridag(i) = -(0.5D0/arg)*(1.0D0-arg2/6.0D0 + (7.0D0/360.0D0)*arg2*arg2)
            detlog     = detlog + log(2.0*arg) - arg + arg2/6.0
            endif
c          print*,'buildit',i,jd(i),jd(i+1),tau,ri(i),ctridag(i)
          enddo
        if (isigma.eq.0) detlog = detlog - 2.0*float(npt)*log(sigma)

        do i=2,npt
          atridag(i) = ctridag(i-1)
          enddo
        btridag(1)   = 1.0-ri(1)    *ctridag(1)
        btridag(npt) = 1.0-ri(npt-1)*ctridag(npt-1)
        do i=2,npt-1
          btridag(i) = 1.0 - ri(i)*ctridag(i) - ri(i-1)*ctridag(i-1)
          enddo

c by making the diagonal element sigma2 it is taken care of
c when the noise matrix is multiplied by sigma2 -- essentially
c     (sigma2*C)^(-1) ( (sigma2*C)^(-1) + N^(-1))(-1) = C^(-1) ( C^(-1) + sigma2 N^(-1) )^(-1)
c this makes the 1/C+1/N matrix closer to unity to avoid the determinant becoming too small a number

        sigma2 = sigma*sigma
        do i=1,npt
          if (isigma.eq.0) then
            diag(i)    = 1.0
            atridag(i) = atridag(i)/sigma2
            btridag(i) = btridag(i)/sigma2
            ctridag(i) = ctridag(i)/sigma2
          else
            diag(i)    = sigma2
            endif
          enddo

c check the results since it costs little -- comment out to save some time
c        do i=1,npt
c          do j=1,npt
c            val = btridag(j)*sigma2*exp(-abs(jd(i)-jd(j))/tau)
c            if (j.ne.1  ) val = val + atridag(j)*sigma2*exp(-abs(jd(i)-jd(j-1))/tau)
c            if (j.ne.npt) val = val + ctridag(j)*sigma2*exp(-abs(jd(i)-jd(j+1))/tau)
c            if (i.eq.j) val = val-1.0
c            if (abs(val).ge.1.e-6) then
c              if (i.eq.j) val = val+1.0
c              print*,'inverse problem for buildcinv at ',i,j,val
c              endif
c            enddo
c          enddo

c        print*,'debug matrix old '
c        print*,'{{ ' ,btridag(1),',',ctridag(1),',0},'
c        print*,' {',atridag(2),',',btridag(2),',',ctridag(2),'},'
c        print*,' {0,',atridag(3),',',btridag(3),'}}'
c              
c        print*,'old={{ ' ,btridag(1),',',ctridag(1),',0,0},'
c        print*,' {',atridag(2),',',btridag(2),',',ctridag(2),',0},'
c        print*,' {0,',atridag(3),',',btridag(3),',',ctridag(3),'},'
c        print*,' {0,0,',atridag(4),',',btridag(4),'}}'
      
        return
        end

c this produces matrices which have the diagonal matrix multiplied into them
c and thus have unit deterinants for the covariance matrix (before doing the sigmas)
      subroutine buildcinvnew(atridag,btridag,ctridag,diag,detlog,jd,npt,sigma,tau,mode)
        implicit real*8 (a-h,o-z)
        PARAMETER (NMAX=1600)
        real*8 atridag(*),btridag(*),ctridag(*),jd(*)
        real*8 diag(*)
        real*8 ri(NMAX),rootc(NMAX)

        integer PR,FOR,NORMAL,NOISE,INFINITE
        common /fittype/PR,FOR,NORMAL,NOISE,INFINITE

c with this structure, including putting sigma2 into the diagonal matrix
c the determinant of C is unity and detlog = 0
        detlog = 0.0


        do i=1,npt-1
          if (jd(i+1).le.jd(i)) then
            print*,' this does not presently work if the data is not in time order '
            print*,' and will have a divergence for equal time data '
            print*,' offender is points ',i,i+1,jd(i),jd(i+1)
            stop
            endif
          if (mode.eq.NORMAL) then
            arg  = abs(jd(i+1)-jd(i))/tau
            ri(i)= exp(-arg)
          else
c just make it like the noise matrix
            if (mode.eq.NOISE) then
              arg   = 1.0
              ri(i) = 0.0
            else
c "limit" of infinite time scale
              tause= 1.e7
              arg  = abs(jd(i+1)-jd(i))/tause
              ri(i)= exp(-arg)
              endif
            endif
          if (arg.ge.1.e-4) then
            rootc(i) = sqrt(1.0-ri(i)*ri(i))
          else
            arg2     = arg*arg
            rootc(i) = sqrt(2.0*arg)*(1.0D0-0.5D0*arg+(5.0D0/24.0D0)*arg2-arg*arg2/16.0D0+(79.0D0/5760.0D0)*arg2*arg2)
            endif
          enddo

c the diagonal
        btridag(1)   = 1.0D0/rootc(1)
        btridag(npt) = 1.0D0/rootc(npt-1)
        do i=2,npt-1
          btridag(i) = (1.0-ri(i-1)*ri(i-1)*ri(i)*ri(i))/rootc(i)/rootc(i-1)
          if (btridag(i).lt.0.0) then
            print*,'something screwy - negative diagonal element ',i,btridag(i)
            print*,'  ri values ',ri(i-1),ri(i)
            print*,'  rootc values ',rootc(i-1),rootc(i)
            stop
            endif
          enddo 
c the lower elements
        atridag(2) = -ri(1)/rootc(1)
        do i=3,npt
          atridag(i) = -ri(i-1)*rootc(i-2)/rootc(i-1)
          enddo
c the upper elements
        ctridag(npt-1) = -ri(npt-1)/rootc(npt-1)
        do i=1,npt-2
          ctridag(i) = -ri(i)*rootc(i+1)/rootc(i)
          enddo
 
c        print*,'debug ',btridag(1),ctridag(1)
c        do i=2,npt-1
c          print*,'debug ',atridag(i),btridag(i),ctridag(i)
c          enddo
c        print*,'debug ',ctridag(npt-1),btridag(npt-1)
c        stop

c        print*,'{r1->',ri(1),',r2->',ri(2),',r3->',ri(3),'}'
c debug 
c        print*,'debug matrix '
c        print*,'new={{ ' ,btridag(1),',',ctridag(1),',0,0},'
c        print*,' {',atridag(2),',',btridag(2),',',ctridag(2),',0},'
c        print*,' {0,',atridag(3),',',btridag(3),',',ctridag(3),'},'
c        print*,' {0,0,',atridag(4),',',btridag(4),'}}'
c
c        print*,'diagonal '
c        print*,'diag={{ ' ,rootc(1),',0,0,0},'
c        print*,' {0,',rootc(1)*rootc(2),',0,0},'
c        print*,' {0,0,',rootc(2)*rootc(3),',0},'
c        print*,' {0,0,0,',rootc(3),'}}'
c        stop

c set up the diagonal matrix to be multiplied into the noise matrix
        sigma2    = sigma*sigma
        diag(1)   = sigma2*rootc(1)
        diag(npt) = sigma2*rootc(npt-1)
        do i=2,npt-1
          diag(i) = sigma2*rootc(i)*rootc(i-1)
          enddo


        return
        end

       
      SUBROUTINE tridag(a,b,c,r,u,n)
      implicit real*8 (a-h,o-z)
      INTEGER n,NMAX
      REAL*8 a(n),b(n),c(n),r(n),u(n)
      PARAMETER (NMAX=1600)
      INTEGER j
      REAL*8 bet,gam(NMAX)
      if(b(1).eq.0.) print*,'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.) print*,'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END

      SUBROUTINE ludcmp(a,n,np,indx,d)
      implicit real*8 (a-h,o-z)
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=1600,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) print*,'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END

      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit real*8 (a-h,o-z)
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END

      subroutine tridiagdet(cinvatrid,cinvbtrid,cinvctrid,npt,detlog)

      implicit real*8 (a-h,o-z)
      real*8 cinvatrid(*),cinvbtrid(*),cinvctrid(*)
      real*8 detlog

c this computes the log of the determinant of a tridiagonal matrix
c now one can compute the determinant using the recursion relation
c    D_n = s_nn D_n-1 - s_n,n-1 s_n-1,n D_n-2
c    D_0 = 1, D_-1 = 0
c however, this will diverge if not careful - so work out the determinant
c of     (tridiagonal)(diagonal)^(-1) = (   1    s12/s11   0..... )
c                                     = (s21/s22    1   s23/s22 0 ...)
c and the log of the determinant of the diagonal matrix, then combine
c supposed to be
c        dnow   = cinvbtrid(1)
c becomes after division (diagonal of the product is just unity)
        dnow   = 1.0
        dold1  = 1.0
        dold2  = 1.0
        detlog = log(cinvbtrid(1))
        do i=2,npt
c supposed to be
c          dnow   = cinvbtrid(i)*dold1 - cinvatrid(i)*cinvctrid(i-1)*dold2
c becomes after division
          dnow   = dold1 - (cinvatrid(i)/cinvbtrid(i))*(cinvctrid(i-1)/cinvbtrid(i-1))*dold2
          dold2  = dold1
          dold1  = dnow
          detlog = detlog + log(cinvbtrid(i))
c          print*,'debug ',cinvatrid(i),cinvbtrid(i),cinvctrid(i-1),dnow
          enddo

c        print*,'debug tridiagdet dnow ',dnow,' old1 ',dold1,' old2 ',dold2,' diagonal ',detlog
        detlog = detlog + log(dnow)

        return
        end


c this uses the same recurrance relation, but uses as its working variables
c the log of the sequential determinants -- also, with this version there
c is no point in taking out the diagonal elements
c     log D_n  = log D_n-1  + log s_nn + log( 1 - s_n,n-1 s_n-1,n D_n-2
c                                                       s_nn      D_n-1 )
c but we can do D_n-2/D_n-1 = exp( log D_n-2 - log D_n-1), which should be 
c much more stable
      subroutine tridiagdetlog(cinvatrid,cinvbtrid,cinvctrid,npt,detlog)

      implicit real*8 (a-h,o-z)
      real*8 cinvatrid(*),cinvbtrid(*),cinvctrid(*)
      real*8 detlog

      dold2  = 0.0
      dold1  = log(cinvbtrid(1))
      dnow   = dold1
      do i=2,npt
        factor = exp(dold2-dold1)
        arg    = 1.0-(cinvatrid(i)*cinvctrid(i-1)/cinvbtrid(i))*factor
        dnow   = dold1 + log(cinvbtrid(i)) + log(arg)
        dold2  = dold1
        dold1  = dnow
c        print*,'debug ',cinvatrid(i),cinvbtrid(i),cinvctrid(i-1),factor,arg,dnow
        enddo

      detlog = dnow

      return
      end


c this computes the likelihood the slow way, using LU
c decompositions, but can be modified for matrices with more
c complex structure
        subroutine slowlike(jd,mag,emag,npt,sigma,tau,plike,chi2,iwrite)
        parameter (NMAX=1600)
        implicit real*8 (a-h,o-z)

        real*8 jd(*),mag(*),emag(*)

        real*8 modmag(NMAX)
        real*8 cmatrix(NMAX,NMAX),cpnmatrix(NMAX,NMAX)
        integer index(NMAX)
        real*8 cplusninvdoty(NMAX)

c the 0.5 is here to make it match the definitions in Kelly 
        variance = sigma*sqrt(0.5*tau/365.0)

c this does it the slow way
c build the covariance plus noise matrix
        do i=1,npt
          do j=1,npt
            cmatrix(i,j)   = variance*variance*exp(-abs(jd(i)-jd(j))/tau)
            cpnmatrix(i,j) = cmatrix(i,j)
            enddo
          cpnmatrix(i,i) = cpnmatrix(i,i) + emag(i)*emag(i)
          enddo
c do the LU decomposition of the matrix
        np = NMAX
        call ludcmp(cpnmatrix,npt,np,index,d)
c now we want cpnmatrix^(-1)*mag = x, which is the same as
c    mag = cpnmatrix*x, so we solve this equation for x
        do i=1,npt
          cplusninvdoty(i) = mag(i)
          enddo
        call lubksb(cpnmatrix,npt,np,index,cplusninvdoty)
c now the chi^2 is just mag.cplusninvdoty
c the determinant of the matrix is just the product of its diagonal elements
c also work out the model light curve
        chi2   = 0.0
        detlog = 0.0
        do i=1,npt
          chi2   = chi2   + mag(i)*cplusninvdoty(i)
          detlog = detlog + log(abs(cpnmatrix(i,i)))
          enddo
c now plike = (-1/2) mag (C+N)^(-1) mag - (1/2) log |C+N|
c so we want
        plike = -0.5*chi2-0.5*detlog
        if (iwrite.eq.1) print*,'chi2 = ',chi2,npt,detlog,plike

c write out the model light curve
        if (iwrite.eq.1) then
          do i=1,npt
            modmag(i) = 0.0
            do j=1,npt
              modmag(i) = modmag(i) + cmatrix(i,j)*cplusninvdoty(j)
              enddo
            write(13,*)jd(i),mag(i),emag(i),modmag(i)
            enddo
          endif

        return 
        end

c this computes the likelihood the fast way, assuming the covariance
c matrix etc has the structure needed to use George Rybicki tridiagonal
c tricks

        subroutine fastlike(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,iwrite,mode)
        parameter (NMAX=1600)
        implicit real*8 (a-h,o-z)

        real*8 jd(*),mag(*),emag(*)

        real*8 cinvatrid(NMAX),cinvbtrid(NMAX),cinvctrid(NMAX)
        real*8 cpninvatrid(NMAX),cpninvbtrid(NMAX),cpninvctrid(NMAX)

        real*8 diag(NMAX)
        real*8 magnoise(NMAX)

        real*8 lpar(nlindim),lmat(NMAX,nlindim),lcovar(nlindim,nlindim)

        real*8 temp1(NMAX),temp2(NMAX)


        real*8 semp1(NMAX,nlindim),evec(NMAX),linnoise(NMAX)
        real*8 semp2(NMAX,nlindim)

        real*8 lcovarinv(nlindim,nlindim)
        real*8 linmod(NMAX),lincor(NMAX)

        integer PR,FOR,NORMAL,NOISE,INFINITE
        common /fittype/PR,FOR,NORMAL,NOISE,INFINITE

        integer index(nlindim)

        do i=1,npt
          magnoise(i) = mag(i)/emag(i)**2
          enddo

c        print*,'debug data {',mag(1),',',mag(2),',',mag(3),',',mag(4),',',mag(5),'}'
c        print*,'debug vec1 {',lmat(1,1),',',lmat(2,1),',',lmat(3,1),',',lmat(4,1),',',lmat(5,1),'}'
c        print*,'debug vec2 {',lmat(1,2),',',lmat(2,2),',',lmat(3,2),',',lmat(4,2),',',lmat(5,2),'}'

c the 0.5 is here to make it match the definitions in Kelly 
        if (mode.eq.NORMAL) then
          variance = sigma*sqrt(0.5*tau/365.0)
        else
          variance = sigma
          endif

c so we are going to evaluate mag^T (C+N)^(-1) mag
c rewrite as  mag^T C^(-1) (C^(-1)+N^(-1))^(-1) N^(-1)mag
c           = mag^T C^(-1) (C^(-1)+N^(-1))^(-1) magnoise
c build the inverse of the covariance matrix 1/C,  in tridiagonal form
c stored as will be used in the NR tridiagonal matrix equation solver
c diag is the diagonal matrix used to make the determinent of C unity -- needs
c to be multiplied into the noise
c        call buildcinv(cinvatrid,cinvbtrid,cinvctrid,diag,detlog1,jd,npt,variance,tau)
        call buildcinvnew(cinvatrid,cinvbtrid,cinvctrid,diag,detlog1,jd,npt,variance,tau,mode)
c        print*,'debug covar^(-1) '
c        print*,'  {{',cinvbtrid(1),',',cinvctrid(1),',0,0,0 },'
c        print*,'   {',cinvatrid(2),',',cinvbtrid(2),',',cinvctrid(2),',0,0},'
c        print*,'   {0,',cinvatrid(3),',',cinvbtrid(3),',',cinvctrid(3),',0},'
c        print*,'   {0,0,',cinvatrid(4),',',cinvbtrid(4),',',cinvctrid(4),'},'
c        print*,'   {0,0,0,',cinvatrid(5),',',cinvbtrid(5),'}}'
c add in the inverse noise to get M=(C^(-1)+N^(-1))
        do i=1,npt
          cpninvatrid(i) = cinvatrid(i)
          cpninvbtrid(i) = cinvbtrid(i) + diag(i)/emag(i)**2
          cpninvctrid(i) = cinvctrid(i)
          enddo
c        print*,'debug covar^(-1) + noise^(-1)'
c        print*,'  {{',cpninvbtrid(1),',',cpninvctrid(1),',0,0,0 },'
c        print*,'   {',cpninvatrid(2),',',cpninvbtrid(2),',',cpninvctrid(2),',0,0},'
c        print*,'   {0,',cpninvatrid(3),',',cpninvbtrid(3),',',cpninvctrid(3),',0},'
c        print*,'   {0,0,',cpninvatrid(4),',',cpninvbtrid(4),',',cpninvctrid(4),'},'
c        print*,'   {0,0,0,',cpninvatrid(5),',',cpninvbtrid(5),'}}'
c remember that magnoise is really N^(-1) mag
c so we now have M=1/C+1/N we need (1/M)(magnoise)=temp1 so rewrite
c as M*temp1 = magnoise, which we solve for r using the NR tridiagonal
c equation solver
        call tridag(cpninvatrid,cpninvbtrid,cpninvctrid,magnoise,temp1,npt)
c        print*,'temp1 ',(temp1(i),i=1,npt)
c where temp1 = (C^(-1)+N^(-1)) N^(-1) mag
c want  temp2 = C^(-1) temp1 = (C+N)^(-1) mag
c this leaves us with doing chi2 = mag^T C^(-1) temp1
        do i=1,npt
          temp2(i) = cinvbtrid(i)*temp1(i)
          if (i.ne.1)   temp2(i) = temp2(i) + cinvatrid(i)*temp1(i-1)
          if (i.ne.npt) temp2(i) = temp2(i) + cinvctrid(i)*temp1(i+1)
          enddo
c        print*,'temp2 ',(temp2(i),i=1,npt)
c now do the same for the linear parameters
        if (nlin.gt.0) then
c         semp1 = (1/C+1/N)^(-1) N^(-1) lmat
c               = (1/C+1/N)^(-1) linnoise
c verified that this way of filling the semp1 matrix workds
          do i=1,nlin
            do j=1,npt
              linnoise(j) = lmat(j,i)/emag(j)**2
              enddo
            call tridag(cpninvatrid,cpninvbtrid,cpninvctrid,linnoise,semp1(1,i),npt)
            enddo
c now need to multiply by 1/C
c         semp2 = (C+N)^(-1) lmat
          do j=1,nlin
            do i=1,npt
              semp2(i,j) = cinvbtrid(i)*semp1(i,j)
              if (i.ne.1)   semp2(i,j) = semp2(i,j) + cinvatrid(i)*semp1(i-1,j)
              if (i.ne.npt) semp2(i,j) = semp2(i,j) + cinvctrid(i)*semp1(i+1,j)
              enddo
            enddo 
          do j=1,nlin
            do k=1,nlin
              lcovarinv(j,k) = 0.0
              lcovar(j,k)    = 0.0
              do i=1,npt
                lcovarinv(j,k) = lcovarinv(j,k) + lmat(i,j)*semp2(i,k)
                enddo
              enddo
            lcovar(j,j) = 1.0
            enddo
c DEBUG
c          do i=1,nlin
c            print*,'debug lcovarinv ',(lcovarinv(i,j),j=1,nlin) 
c            enddo
c invert to get the covariance matrix lcovar = (L^T (C+N)^(-1) L)^(-1)
c also determine the determinant of the linear parameter covariance matrix
          call ludcmp(lcovarinv,nlin,nlindim,index,d)
          detlin = 0.0
          do i=1,nlin
            call lubksb(lcovarinv,nlin,nlindim,index,lcovar(1,i))
            detlin = detlin + log(lcovarinv(i,i))
c            print*,'debug lcovarinv after decomposition, detlin ',lcovarinv(i,i),detlin
            enddo
c work out the linear coefficients
          do i=1,nlin 
            lpar(i) = 0.0
            do j=1,nlin
              do k=1,npt
                lpar(i) = lpar(i) + lcovar(i,j)*lmat(k,j)*temp2(k)
                enddo
              enddo
c            print*,'got lpar ',i,lpar(i)
            enddo
c work out lmat*lpar
          do i=1,npt
            linmod(i) = 0.0
            lincor(i) = 0.0
            do j=1,nlin
              linmod(i) = linmod(i) +  lmat(i,j)*lpar(j)
              lincor(i) = lincor(i) + semp2(i,j)*lpar(j)
              enddo
            enddo
          chi2 = 0.0
          do i=1,npt
            chi2 = chi2 + temp2(i)*(mag(i)-linmod(i))
            enddo
c          print*,'got chi2 = ',chi2
        else
          detlin = 0.0
          chi2   = 0.0
          do i=1,npt
            chi2 = chi2 + temp2(i)*mag(i)
            enddo
          endif
c now compute the determinant -- need the determinants of
c this way of doing the determinants has problems for large tau -- use version in buildcinv
c        print*,'doing detlog1 '
c        call tridiagdet(cinvatrid  ,cinvbtrid  ,cinvctrid  ,npt,detlog1old)
c        print*,'done detlog1, got new ',detlog1,' old ',detlog1old
c        print*,'doing detlog2 '
c this version has accumulating errors that can make it become infinite
c        call tridiagdet(cpninvatrid,cpninvbtrid,cpninvctrid,npt,detlog2old)
        call tridiagdetlog(cpninvatrid,cpninvbtrid,cpninvctrid,npt,detlog2)
c        print*,'done detlog2, got ',detlog2,detlog2old
        detlog3 = 0.0
        do i=1,npt
          detlog3 = detlog3 - 2.0*log(emag(i))
          enddo
c detlin = log | L^T C^(-1) L |
c this is log |C+N| = log|1/C+1/N| - log|1/C| - log|1/N| 
        detlog = detlog2 - detlog1 - detlog3
        plike = -0.5*chi2 - 0.5*(detlog+detlin)
        if (iwrite.eq.1) print*,'sigma,tau',sigma,tau
        if (iwrite.eq.1) print*,' 1/noise ',detlog3,' 1/covariance ',detlog1,' 1/noise + 1/covar ',detlog2
        if (iwrite.eq.1) print*,'fastlike: like,chi2,detlog,detlin ',plike,chi2,detlog,detlin


        return
        end

c predicts the light curve
        subroutine predict(jd,mag,emag,id,npt,sigma,tau,shift,scale,jwant,mwant,ewant,nwant)
        parameter (NMAX=1600)
        implicit real*8 (a-h,o-z)

        real*8 jd(*),mag(*),emag(*)
        real*8 jwant(*),mwant(*),ewant(*)

        real*8 cmatrix(NMAX,NMAX),cpnmatrix(NMAX,NMAX)
        integer index(NMAX)
        real*8 cplusninvdoty(NMAX)

        real*8 covar(NMAX),temp(NMAX)

        integer id(*)

        real*8 shift(*),scale(*)

        print*,'will predict ',nwant,' points '

        variance  = sigma*sqrt(0.5*tau/365.0)
        variance2 = variance*variance

c this does it the slow way
c build the covariance plus noise matrix
        do i=1,npt
          do j=1,npt
            cmatrix(i,j)   = scale(id(i))*scale(id(j))*variance2*exp(-abs((jd(i)-shift(id(i)))-(jd(j)-shift(id(j))))/tau)
            cpnmatrix(i,j) = cmatrix(i,j)
            enddo
          cpnmatrix(i,i) = cpnmatrix(i,i) + emag(i)*emag(i)
          enddo
c        print*,'built matrix '
c do the LU decomposition of the matrix
        np = NMAX
        call ludcmp(cpnmatrix,npt,np,index,d)
c now we want cpnmatrix^(-1)*mag = x, which is the same as
c    mag = cpnmatrix*x, so we solve this equation for x
        do i=1,npt
          cplusninvdoty(i) = mag(i)
          enddo
        call lubksb(cpnmatrix,npt,np,index,cplusninvdoty)

        print*,'done rebuilding model '

c computes the model light curve C(t) (C+N)^(-1) mag
c                          error C(t) (C+N)^(-1) C(t)
c note -- it might be possible to do the error a faster way..
        do i=1,nwant
          do j=1,npt
            covar(j)   = variance2*scale(id(j))*exp(-abs((jd(j)-shift(id(j)))-jwant(i))/tau)
            temp(j)    = covar(j)
            enddo
          call lubksb(cpnmatrix,npt,np,index,temp)
          mwant(i) = 0.0
          ewant(i) = variance2
          do j=1,npt
            mwant(i) = mwant(i) + covar(j)*cplusninvdoty(j)
            ewant(i) = ewant(i) - covar(j)*temp(j)
            enddo
          ewant(i) = sqrt(ewant(i))
          enddo

        return 
        end

c this computes the likelihood with linear parameters the slow way, using LU
c decompositions, but can be modified for matrices with more
c complex structure
        subroutine slowlikelin(jd,mag,emag,id,npt,sigma,tau,shift,scale,lmat,lpar,lcovar,nlin,nlindim,plike,chi2)
        parameter (NMAX=1600)
        implicit real*8 (a-h,o-z)

        real*8 jd(*),mag(*),emag(*)

        real*8 modmag(NMAX)
        real*8 cmatrix(NMAX,NMAX),cpnmatrix(NMAX,NMAX)
        real*8 cpninvmatrix(NMAX,NMAX)
        integer index(NMAX)
        real*8 cplusninvdoty(NMAX)
        integer id(NMAX)

        real*8 lpar(nlindim),lmat(NMAX,nlindim),lcovar(nlindim,nlindim)
        real*8 tlmat1(NMAX,nlindim),tlmat3(NMAX,nlindim)
        real*8 tlmat2(nlindim,nlindim),tlmat2inv(nlindim,nlindim)

        real*8 cmat(nlindim,nlindim),cvec(nlindim)

        real*8 scale(*),shift(*)

c the 0.5 is here to make it match the definitions in Kelly 
        variance = sigma*sqrt(0.5*tau/365.0)

c this does it the slow way
c build the covariance plus noise matrix
        do i=1,npt
          do j=1,npt
            cmatrix(i,j)   = scale(id(i))*scale(id(j))*variance*variance*exp(-abs((jd(i)-shift(id(i)))-(jd(j)-shift(id(j))))/tau)
            cpnmatrix(i,j) = cmatrix(i,j)
            enddo
          cpnmatrix(i,i) = cpnmatrix(i,i) + emag(i)*emag(i)
          enddo

c        print*,'debug = built matrix nlin = ',nlin,' maximum ',nlindim 

c if linear parameters need inverse of cpnmatrix
        if (nlin.gt.0) then
          do i=1,npt
            do j=1,npt
              cpninvmatrix(i,j) = 0.0
              enddo
            cpninvmatrix(i,i) = 1.0
            enddo
          np = NMAX
          call ludcmp(cpnmatrix,npt,np,index,d)
          detlog = 0.0
          do j=1,npt
            detlog = detlog + log(abs(cpnmatrix(j,j)))
            call lubksb(cpnmatrix,npt,np,index,cpninvmatrix(1,j))
            enddo
c work out the optimal linear parameters
c   lpar = (lmat^T * cpninvmatrix * lmat)^(-1) lmat^T cpninvmatrix mag
c lets work out tlmat1 = lmat^T cpninvmatrix 
          do i=1,nlin
            do j=1,npt
              tlmat1(j,i) = 0.0
              do k=1,npt
                tlmat1(j,i) = tlmat1(j,i) + lmat(k,i)*cpninvmatrix(k,j)
                enddo
              enddo
            enddo
c          print*,'debug done 1 '
c now work out lpar = lmat^T cpninvmatrix mag = tlmat1 mag
          do i=1,nlin
            lpar(i) = 0.0 
            do j=1,npt
              lpar(i) = lpar(i) + tlmat1(j,i)*mag(j)
              enddo
c DEBUG
c            cvec(i) = lpar(i)
            enddo
c now workout tlmat2 = lmag^T cpninvmatrix lmat = tlmat1 lmat
          do i=1,nlin
            do j=1,nlin
              tlmat2(i,j) = 0.0
              do k=1,npt
                tlmat2(i,j) = tlmat2(i,j) + lmat(k,j)*tlmat1(k,i)
                enddo
              enddo
            enddo
c DEBUG -- save matrix to check inverse
c          do i=1,nlin
c            do j=1,nlin
c              cmat(i,j) = tlmat2(i,j)
c              enddo
c            enddo
c do the LU decomposition of the matrix
          np = NMAX
          call ludcmp(tlmat2,nlin,nlindim,index,d)
          detlin = 0.0
          do i=1,nlin
            detlin = detlin + log(abs(tlmat2(i,i)))
            enddo
c LU substitute to get the linear parameters
          call lubksb(tlmat2,nlin,nlindim,index,lpar)
c now if we just wanted to solve for the fit, we could use chi^2 = mag^T cpninvmatrix^(-1) (mag-lmat*lpar)
c but we also need the determinant, which means updating the matrix
c       cpninvmatrix -> cpninvmatrix - tlmat1 lcovar  tlmat1 
c so we actually need to invert tlmat2 (which will also be the covariance matrix of the linear parameters)
c note, however, that this could be done faster if we do not want the covariance matrix by simultaneously
c doing lcovar  tlmat1
          do i=1,nlin
            do j=1,nlin
              lcovar(i,j) = 0.0
              enddo
            lcovar(i,i) = 1.0
            enddo
          do j=1,nlin
            call lubksb(tlmat2,nlin,nlindim,index,lcovar(1,j))
            enddo
c DEBUG -- check that this is right
c          do i=1,nlin
c            do j=1,nlin
c              check = 0.0
c              do k=1,nlin
c                check = check + cmat(i,k)*lcovar(k,j)
c                enddo
c              print*,'inverse check ',i,j,check,lcovar(i,j)
c              enddo
c            enddo
c now compute lcovar *tlmat1
          do i=1,nlin
            do j=1,npt
              tlmat3(j,i) = 0.0
              do k=1,nlin
                tlmat3(j,i) = tlmat3(j,i) + lcovar(i,k)*tlmat1(j,k)
                enddo
              enddo
            enddo
c DEBUG check this is right -- if it acts on mag it should give the linpars
c          do i=1,nlin
c            check = 0.0
c            check2= 0.0
c            do j=1,npt
c              check = check + tlmat3(j,i)*mag(j)
c              enddo
c            do j=1,nlin
c              check2 = check2 + lcovar(i,j)*cvec(j)
c              enddo
c            print*,' linpar ',i,' new/old ',check,check2,lpar(i)
c            enddo  
c DEBUG - get the chi^2 you would have without removing the linear functions
c          chi20  = 0.0
c          do i=1,npt
c            do j=1,npt
c              chi20   = chi20   + mag(i)*mag(j)*cpninvmatrix(i,j)
c              enddo
c            enddo
c now update the matrix
          do i=1,npt
            do j=1,npt
              do k=1,nlin
                cpninvmatrix(i,j) = cpninvmatrix(i,j) - tlmat1(i,k)*tlmat3(j,k)
                enddo
              enddo
            enddo 
c determine the chi^2
          chi2  = 0.0
          do i=1,npt
            do j=1,npt
              chi2   = chi2   + mag(i)*mag(j)*cpninvmatrix(i,j)
              enddo
            enddo
c          print*,'debug: nolinear ',chi20,' with linear ',chi2
        else 
c NO LINEAR PARAMETERS
c do the LU decomposition of the matrix
          np = NMAX
          call ludcmp(cpnmatrix,npt,np,index,d)
c now we want cpnmatrix^(-1)*mag = x, which is the same as
c    mag = cpnmatrix*x, so we solve this equation for x
          do i=1,npt
            cplusninvdoty(i) = mag(i)
            enddo
          call lubksb(cpnmatrix,npt,np,index,cplusninvdoty)
c now the chi^2 is just mag.cplusninvdoty
c the determinant of the matrix is just the product of its diagonal elements
          chi2   = 0.0
          detlog = 0.0
          detlin = 0.0
          do i=1,npt
            chi2   = chi2   + mag(i)*cplusninvdoty(i)
            detlog = detlog + log(abs(cpnmatrix(i,i)))
            enddo
          endif

c now plike = (-1/2) mag (C+N)^(-1) mag - (1/2) log |C+N|
c so we want
        plike = -0.5*chi2-0.5*(detlog+detlin)
        print*,'slowlinlike: like, chi2, detlog, detlin ',plike,chi2,detlog,detlin 
c        if (iwrite.eq.1) print*,'chi2 = ',chi2,npt,detlog,plike


        return 
        end

      SUBROUTINE indexx(n,arr,indx)
      implicit real*8 (a-h,o-z)
      INTEGER n,indx(n),M,NSTACK
      REAL*8 arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK)  print*,'NSTACK too small in indexx'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END

        subroutine getlimit(tl,ptl,nt,tlow,thig,best)
          implicit real*8 (a-h,o-z)
          real*8 ptl(*),tl(*)
          real*8 ptlsum(1000)

          ptlsum(1) = ptl(1)
          do i=2,nt
            ptlsum(i) = ptlsum(i-1) + ptl(i)
            enddo
          do i=1,nt
            ptlsum(i) = ptlsum(i)/ptlsum(nt)
            enddo
          tlow = ptl(1)
          thig = ptl(nt)
          do i=1,nt-1
            if ((ptlsum(i).le.0.1585).and.(ptlsum(i+1).gt.0.1585)) then
              tlow = tl(i) + (tl(i+1)-tl(i))*(0.1585-ptlsum(i))/(ptlsum(i+1)-ptlsum(i))
              endif
            if ((ptlsum(i).le.0.500).and.(ptlsum(i+1).gt.0.500)) then
              best = tl(i) + (tl(i+1)-tl(i))*(0.5000-ptlsum(i))/(ptlsum(i+1)-ptlsum(i))
              endif
            enddo
          do i=nt,2,-1
            if ((ptlsum(i-1).le.0.8415).and.(ptlsum(i).gt.0.8415)) then
              thig = tl(i-1) + (tl(i)-tl(i-1))*(0.8415-ptlsum(i-1))/(ptlsum(i)-ptlsum(i-1))
              endif
            enddo

c          print*,'tlow,thig ',tlow,thig
c          do i=1,nt
c            write(41,*)tl(i),ptl(i),ptlsum(i)
c            enddo
c          stop

          return
          end

c use the forecasting approach -- assume you have the linear model from the
c full results
      subroutine forecast(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,iwrite)
           parameter (NMAX=1600)
           implicit real*8 (a-h,o-z)

           real*8 jd(*),mag(*),emag(*),lpar(*)
           real*8 lcovar(nlindim,nlindim)
           real*8 lmat(NMAX,nlindim)

           real*8 muse(NMAX)
           real*8 omega,omegaold

           omega0 = 0.5*sigma*sigma*(tau/365.0)

c remove the trends
           if (nlin.gt.0.0) then
             do i=1,npt
               muse(i) = mag(i)
               do j=1,nlin
                 muse(i) = muse(i) - lpar(j)*lmat(i,j)
                 enddo
               enddo
             endif

c do the forecasting 
           xstarold = muse(1)
c initialize using the 1 point result from Rybick/Press
c Kelly would use xhatold = 0 and omegaold = variance**2
           xhatold  = muse(1)*omega0/(omega0+emag(1)**2)
           omegaold = 1.0/(1.0/omega0+1.0/emag(1)**2)
           chi2     = 0.0
           plike    = 0.0
           do i=2,npt
             xstar = muse(i)
             afac  = exp(-(jd(i)-jd(i-1))/tau)
             denom = omegaold+emag(i-1)**2
             xhat  = afac*(xhatold + omegaold*(xstarold-xhatold)/denom)
             omega = omega0*(1.0-afac*afac) + afac*afac*omegaold*(1.0-omegaold/denom)
             chi2  = chi2  + (xhat-xstar)*(xhat-xstar)/denom
             plike = plike - 0.5*log(denom) - 0.5*(xhat-xstar)*(xhat-xstar)/denom
c save the lat values
             xstarold = xstar
             omegaold = omega
             xhatold  = xhat
             enddo 

           return
           end




c generates relationship between proper motion distances and redshift


	subroutine distinit

           implicit real*8 (a-h,o-z)

c om = density of normal matter, ol = cosmological constant
c ok = 1 - om - ok = curvarture density
c ztab(i) = list of redshifts, dtab(i) = proper motion distance to redshift
c ztabs, dtabs = spline coefficients

           real*8 om,ol,ok,ztab(100),dtab(100),ztabs(100),dtabs(100)
           integer ntab
           common /cosmo/om,ol,ok,ztab,dtab,ztabs,dtabs,ntab

           real*8 y(1),zmax,aokh,eps,hgus

           external dsfunc

           zmax = 5.0D0

           ok   = 1.0D0 - om - ol
           aokh = sqrt(abs(ok))

           nn   = 1
           eps  = 1.0D-6
           dz   = zmax/float(ntab-1)
           hgus = 0.1*dz 
           nmax = 1000
           do 100 i=1,ntab
             ztab(i) = dz*float(i-1)
             if (i.eq.1) then
               dtab(i) = 0.0D0
               y(i)    = 0.0D0
             else
               z1     = ztab(i-1)
               z2     = ztab(i)
               nbad   = 0
               nsteps = 0
               call difsub(dsfunc,nn,z1,z2,y,eps,hgus,nsteps,nbad,nmax)
               if (aokh.lt.1.0e-3) then
                 dtab(i)= y(1)*(1.0+ok*y(1)**2/6.0)
               else
                 if (ok.gt.0.0) then
                   dtab(i) = sinh(aokh*y(1))/aokh
                 else
                   dtab(i) = sin (aokh*y(1))/aokh
                   endif
                 endif
               endif
100          continue

c set up spline for converting z to d
           yp1 = 1.0e32
           yp2 = 1.0e32
           call spline(ztab,dtab,ntab,yp1,yp2,dtabs)
           call spline(dtab,ztab,ntab,yp1,yp2,ztabs)

           return
           end

        subroutine dsfunc(x,y,dy)
          implicit real*8 (a-h,o-z)
          real*8 x,y(*),dy(*)
          real*8 om,ol,ok,ztab(100),dtab(100),ztabs(100),dtabs(100)
          integer ntab
          common /cosmo/om,ol,ok,ztab,dtab,ztabs,dtabs,ntab

          dy(1) = 1.0/sqrt((1.0D0+om*x)*(1.0D0+x)**2-x*(2.0D0+x)*ol)

          return
          end

        function asinh(x)
          implicit real*8 (a-h,o-z)
          asinh = log(x+sqrt(1.0D0+x*x))
          return
          end

      SUBROUTINE DIFSUB(DF,NN,A,B,Y,EPS,HGUS,NSTEPS,NBAD,NMAX)
C
C     Note that after a call, the routine is ready to continue integrating
C     from the point it left off.
C     DF = SUBROUTINE DF(X,Y,DY). Computes rhs of DY = F(X,Y). Must be
C          declared external in calling program. Y and DY must be dimensioned.
C     NN = No. of eqns. If NN > 50, change dimension statement.
C     A  = Starting value of X. Updated to B on return.
C     B  = Final value of X. B < A is acceptable.
C     Y  = Solution vector. Must be initialized.
C     EPS= Accuracy parameter. It is a relative error tolerance per step for
C          Y's greater than 1, an absolute error tolerance for those less than
C          1. (see coding near statement 84).
C     HGUS=Initial guess for stepsize. If set to 0, uses default of .01*(B-A).
C          Returns last stepsize used.
C     NSTEPS= No. of steps taken. (Initialize to zero.)
C     NBAD= No. of steps that had to be retaken. (initialize to zero.)
C     NMAX= Max. no. of allowed steps (NSTEPS + NBAD).
C
      DOUBLE PRECISION A,E,H,HHH,ONE,T4,ABSH,HAIM,HSIGN,T1,TSAVE
      DOUBLE PRECISION B,EPS,HGUS,T2,TWO,DUM,ERRMAX,HH,T3,ZR
      DOUBLE PRECISION DY1,DYDUB,DYSAVE,Y,Y2,Y3,YDUB,YSAVE,YSING
      DIMENSION YSAVE(100),DYSAVE(100),
     1 YSING(100),YDUB(100),Y2(100),
     1 Y3(100),DY1(100),DYDUB(100)
      DIMENSION Y(NN)
      LOGICAL END
      DATA ZR/0.0E0/,ONE/1.0E0/,TWO/2.0E0/
      END=.FALSE.
      N=NN
      IF(HGUS.EQ.ZR)HGUS=(B-A)*0.01E0
      HSIGN=SIGN(ONE,B-A)
      HAIM=ABS(HGUS)
      DO 10 J=1,N
   10 YSAVE(J)=Y(J)
      TSAVE=A
   12 CALL DF(TSAVE,YSAVE,DYSAVE)
   11 H=SIGN(HAIM,HSIGN)
      IF(HAIM.GT.ABS(B-TSAVE))THEN
         H=B-TSAVE
         END=.TRUE.
      ENDIF
C HERE TAKE A STEP GIVING YSING,YDUB,DYDUB. BEGIN.
      HH=0.5E0*H
      HHH=0.25E0*H
      T1=TSAVE+HHH
      T2=T1+HHH
      T3=T2+HHH
      T4=TSAVE+H
      DO 61 J=1,N
   61 Y2(J)=YSAVE(J)+HH*DYSAVE(J)
      CALL DF(T2,Y2,DY1)
      DO 62 J=1,N
      Y3(J)=YSAVE(J)+HH*DY1(J)
   62 Y2(J)=Y2(J)+Y3(J)*2.E0
      CALL DF(T2,Y3,DY1)
      DO 63 J=1,N
      Y3(J)=YSAVE(J)+H*DY1(J)
   63 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T4,Y3,DY1)
      DO 64 J=1,N
      YSING(J)=(Y2(J)-YSAVE(J)+HH*DY1(J))/3.0E0
   64 Y2(J)=YSAVE(J)+HHH*DYSAVE(J)
      CALL DF(T1,Y2,DY1)
      DO 72 J=1,N
      Y3(J)=YSAVE(J)+HHH*DY1(J)
   72 Y2(J)=Y2(J)+Y3(J)*2.E0
      CALL DF(T1,Y3,DY1)
      DO 73 J=1,N
      Y3(J)=YSAVE(J)+HH*DY1(J)
   73 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T2,Y3,DY1)
      DO 74 J=1,N
   74 YDUB(J)=(Y2(J)-YSAVE(J)+HHH*DY1(J))/3.0E0
      CALL DF(T2,YDUB,DYDUB)
      DO 81 J=1,N
   81 Y2(J)=YDUB(J)+HHH*DYDUB(J)
      CALL DF(T3,Y2,DY1)
      DO 82 J=1,N
      Y3(J)=YDUB(J)+HHH*DY1(J)
   82 Y2(J)=Y2(J)+Y3(J)*2.E00
      CALL DF(T3,Y3,DY1)
      DO 83 J=1,N
      Y3(J)=YDUB(J)+HH*DY1(J)
   83 Y2(J)=Y2(J)+Y3(J)
      CALL DF(T4,Y3,DYDUB)
      ERRMAX=ZR
      ABSH=ABS(H)
      DO 84 J=1,N
      YDUB(J)=(Y2(J)-YDUB(J)+HHH*DYDUB(J))/3.0E0
C END OF TAKE A STEP BLOCK.
      DUM=DMAX1(ONE,ABS(YDUB(J)))
      E=ABS(YDUB(J)-YSING(J))/(15.E0*DUM)
      IF(E.GT.ERRMAX)ERRMAX=E
   84 CONTINUE
      ERRMAX=ERRMAX/EPS
      IF(ERRMAX.LE.ZR)THEN
         HAIM=ABSH*4.0
      ELSEIF(ERRMAX.LE.ONE)THEN
         HAIM=ABSH*ERRMAX**(-0.2)*.90
      ELSE
         HAIM=ABSH*ERRMAX**(-0.25)*.90
         NBAD=NBAD+1
         END=.FALSE.
         IF(NSTEPS+NBAD.GE.NMAX)GO TO 99
         GO TO 11
      ENDIF
C  STEP SUCCEEDED...
      NSTEPS=NSTEPS+1
      TSAVE=TSAVE+H
      IF(END)THEN
         A=TSAVE
         HGUS=HAIM
         DO 91 J=1,N
   91    Y(J)=(16.E0*YDUB(J)-YSING(J))*6.666666666666666E-2
         RETURN
      ELSE
         IF(NSTEPS+NBAD.GE.NMAX)GO TO 99
         DO 31 J=1,N
   31    YSAVE(J)=(16.E0*YDUB(J)-YSING(J))*6.666666666666666E-2
         GO TO 12
      ENDIF
c   99 PAUSE'NMAX EXCEEDED'
   99 print*,'NMAX EXCEEDED'
      stop
      END


      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)
      implicit real*8 (a-h,o-z)
      INTEGER N
      REAL*8 X,Y,XA(*),Y2A(*),YA(*)
      INTEGER K,KHI,KLO
      REAL*8 A,B,H
      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
         K=(KHI+KLO)/2
         IF(XA(K).GT.X)THEN
            KHI=K
         ELSE
            KLO=K
         ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.)  print*,'bad xa input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *         ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.
      RETURN
      END

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
      implicit real*8 (a-h,o-z)
      INTEGER N,NMAX
      REAL*8 YP1,YPN,X(*),Y(*),Y2(*)
      PARAMETER (NMAX=500)
      INTEGER I,K
      REAL*8 P,QN,SIG,UN,U(NMAX)
      IF (YP1.GT..99E30) THEN
         Y2(1)=0.
         U(1)=0.
      ELSE
         Y2(1)=-0.5
         U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
         SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
         P=SIG*Y2(I-1)+2.
         Y2(I)=(SIG-1.)/P
         U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *         /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT..99E30) THEN
         QN=0.
         UN=0.
      ELSE
         QN=0.5
         UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)
      DO 12 K=N-1,1,-1
         Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE
      RETURN
      END

c gets dls given dos,dol and omega
        function getdls(dol,dos)
          implicit real*8 (a-h,o-z)

          real*8 om,ol,ok,ztab(100),dtab(100),ztabs(100),dtabs(100)
          integer ntab
          common /cosmo/om,ol,ok,ztab,dtab,ztabs,dtabs,ntab

          getdls = sqrt(1.0+ok*dol)*dos-sqrt(1.0+ok*dos)*dol

          return
          end



        subroutine gencurve(jd,mag,emag,npt,tau,sigma,avgmag,iseed)
          implicit real*8 (a-h,o-z)
          real*8 jd(*),mag(*),emag(*)
          real gasdev
          real*8 omega0

          omega0 = 0.5*sigma*sigma*(tau/365.0)

          do i=1,npt
            if (i.eq.1) then
              signal = sqrt(omega0)*gasdev(iseed)
            else
              alpha  = exp(-(jd(i)-jd(i-1))/tau)
              signal = alpha*signal + sqrt(omega0*(1.0-alpha*alpha))*gasdev(iseed)
              endif
            mag(i) = avgmag + signal + emag(i)*gasdev(iseed) 
            enddo

        end


      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      REAL gasdev
CU    USES ran1
      INTEGER iset
      REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

        
       subroutine dofit(idq,jd,mag,emag,npt,avgin,lmat,lpar,lcovar,nlin,nlindim,
     1    ifit,iout,pbest,cbest,tbest,tlow,thig,sbest,slow,shig)
          parameter (NMAX=1600)

          implicit real*8 (a-h,o-z)
          real*8 jd(*),mag(*),emag(*)
          
          integer PR,FOR,NORMAL,NOISE,INFINITE
          common /fittype/PR,FOR,NORMAL,NOISE,INFINITE

          real*8 psave(251,251),tlsave(251),slsave(251)
          real*8 psave2(251,251)
          real*8 ptl(251),psl(251)
          real*8 ptlsum(251),pslsum(251)
          real*8 xm(251),pm(251)

          real*8 lpar(nlindim),lmat(NMAX,nlindim),lcovar(nlindim,nlindim)

          character*100 fname
          common /object/fname

          integer igrid
          common /flags/igrid

          common /limfits/pnoise,cnoise,snoise,pinf,cinf,sinf
      
          if (npt.lt.10) then
            iedge  =   0
            plike  = -10
            chi2   = -10
            cdof   = -10
            tbest  = -10
            sbest  = -10
            tlow   = -10
            slow   = -10
            thig   = -10
            shig   = -10
            avgmag = avgin
            smedian= -10
            tmedian= -10
            write(iout,'(i2,12(1x,f13.6),1x,a65)')iedge,avgmag,pbest,cbest,cdof,sbest,slow,shig,tbest,tlow,thig,tmedian,smedian,fname
            return
            endif

c          do i=1,5
c            write(*,'(10(1x,f13.6))')jd(i),mag(i),emag(i),lmat(i,1)
c            enddo

          if (ifit.eq.PR) then
            print*,'doing PR fit on ',npt,' points with ',nlin,' linear params '
          else
            print*,'doing FOR fit on ',npt,' points with ',nlin,' linear params '
            endif


c before doing the 2D fit, work out the best fit model using either just a noise matrix (limit tau->0)
c or infinite correlation (limit tau->infinity), but only for the PR case
          pnoise = -10
          cnoise = -10
          snoise = -10
          pinf   = -10
          cinf   = -10
          sinf   = -10
          if (ifit.eq.PR) then
            smin   = 0.001
            smax   = 1000.0
            sminl  = log10(smin)
            smaxl  = log10(smax)
            ns     = 241
            dsl    = (smaxl-sminl)/float(max(ns-1,1))
            tau    = 1.0
c best fitting added noise model
            print*,'doing noise model'
            pnoise = -1.e32
            do j=1,ns
              sl        = sminl + float(j-1)*dsl
              slsave(j) = sl
              sigma     = 10.0**sl 
              call fastlike(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,0,NOISE)
              psave(1,j) = plike
              if (plike.ge.pnoise) then
                pnoise = plike
                cnoise = chi2
                snoise = sl
                inoise = j
                endif
              enddo
c best fitting infinite correlation model
            print*,'doing infinite correlation '
            pinf = -1.e32
            do j=1,ns
              sl        = sminl + float(j-1)*dsl
              slsave(j) = sl
              sigma     = 10.0**sl
              call fastlike(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,0,INFINITE)
              psave(1,j) = plike
              if (plike.ge.pinf) then
                pinf = plike
                cinf = chi2
                sinf = sl
                iinf = j
                endif
              enddo
            endif

          tmin   = 0.1
          tmax   = 1.0e5
          tminl  = log10(tmin)
          tmaxl  = log10(tmax)

          smin   = 0.001
          smax   = 10.0
          sminl  = log10(smin)
          smaxl  = log10(smax)

          pbest  = -1.e32

          nshift = 0
c this will do the initial range at 25% resolution

100       continue
          nt     = 61
          dtl    = (tmaxl-tminl)/float(max(nt-1,1))
          ns     = 81
          dsl    = (smaxl-sminl)/float(max(ns-1,1))

          do i=1,nt
            tl        = tminl + float(i-1)*dtl
            tlsave(i) = tl
            tau       = 10.0**tl
            do j=1,ns
              sl        = sminl + float(j-1)*dsl
              slsave(j) = sl
              sigma     = 10.0**sl 
              if (ifit.eq.PR) then
                call fastlike(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,0,NORMAL)
              else
                call forecast(jd,mag,emag,npt,sigma,tau,lmat,lpar,lcovar,nlin,nlindim,plike,chi2,0)
                endif
              if ((igrid.eq.1).and.(nshift.eq.0).and.(idq.gt.0)) write(14,'(3(1x,i4),6(1x,g13.6),1x,i4)')idq,i-1,j-1,tl,sl,plike,plike,chi2,chi2,npt
c              if ((igrid.eq.1).and.(nshift.eq.1)) write(14,'(3(1x,i4),6(1x,g13.6),1x,i4,)')2,i-1,j-1,tl,sl,plike,plike,chi2,chi2,npt
c              if ((igrid.eq.1).and.(nshift.eq.2)) write(14,'(3(1x,i4),6(1x,g13.6),1x,i4,)')3,i-1,j-1,tl,sl,plike,plike,chi2,chi2,npt
              psave(i,j)  = plike
              if (plike.gt.pbest) then
                pbest = plike
                cbest = chi2
                sbest = sl
                tbest = tl
                abest = lpar(1)
                imax  = i
                jmax  = j
                endif
              enddo
            enddo
           if (igrid.eq.1) then
             print*,'pbest ',pbest,' chibest ',cbest,' for ',npt
             print*,'tbest ',tbest,' sbest ',sbest
             endif

c if on the edge of the grid, shift the grid by a decade and rerun
c do this up to two times
          if ((nshift.lt.2).and.(igrid.eq.0)) then
            iredo  = 0
            if (sbest.le.sminl+1) then
              iredo = 1
              sminl = sminl - 0.5
              smaxl = smaxl - 0.5
              endif
c only allaw one shift to higher amplitudes, there is a solution
c branch corresponding to just rescaling the errors with
c large sigma and large tau 
            if ((sbest.ge.smaxl-1).and.(nshift.eq.0)) then
              iredo = 1
              sminl = sminl + 0.5
              smaxl = smaxl + 0.5
              endif
            if (tbest.le.tminl+1) then
              iredo = 1
              tminl = tminl - 1
              tmaxl = tmaxl - 1
              endif
            if (tbest.ge.tmaxl-1) then
              iredo = 1
              tminl = tminl + 1
c do not shift to higher maximum times than 10^5 years
c              tmaxl = tmaxl + 1
c but improve sigma resolution
              sminl = sbest - 1.0
              smaxl = sbest + 1.0
              endif
            if (iredo.eq.1) then
              print*,'shifted grid to ',tminl,tmaxl,sminl,smaxl
              print*,'given max like ',pbest,' at ',tbest,sbest
              nshift = nshift + 1
              go to 100
              endif
            endif

          iedge = 0
          if ((imax.le.3).or.(imax.ge.nt-2)) then
            print*,' peak on edge of tau grid - no interpolation '
            iedge = 1
            endif
          if ((jmax.le.3).or.(jmax.ge.ns-2)) then
            print*,' peak on edge of sigma grid - no interpolation '
            iedge = 1
            endif
c          print*,'sbest,tbest ',sbest,tbest
c          print*,'sgrid ',sminl,smaxl,' tgrid ',tminl,tmaxl
c interpolation to get the max like
          if (iedge.eq.0) then
            pintmax = -1.e32
            do i=1,nt
              ptl(i) = -1.e32
              jintmax = 1
              pjmax   = psave(i,1)
c              print*,'doing tau ',tlsave(i)
              do j=2,ns-2
                p0 = psave(i,j)
                pl = psave(i,j-1)
                ph = psave(i,j+1)
                if (p0.ge.pjmax) then
                  pjmax   = p0
                  jintmax = j
                  endif
                if ((p0.ge.pl).and.(p0.ge.ph)) then
                  x0 = slsave(j)
                  xl = slsave(j-1)
                  xh = slsave(j+1)
                  dx = xh-x0
c fit as quadratic and maximize
                  xm(i) = x0 + 0.5*dx*(ph-pl)/(2.0*p0-ph-pl)
                  pm(i) = p0 + 0.125*(ph-pl)*(ph-pl)/(2.0*p0-ph-pl)
c                  print*,'sigmas ',xl,x0,xh,'like   ',pl,p0,ph,' pout ',xm(i),pm(i)
                  if (pm(i).gt.pintmax) then
                    pintmax = pm(i)
                    iintmax = i
                    endif
                  endif
                enddo
c              print*,'max was at ',jintmax,pjmax
              enddo
c now interpolate along the new quadratic (unless on edge)
            if ((iintmax.ne.1).or.(iintmax.ne.nt)) then
              p0 = pm(iintmax)
              pl = pm(iintmax-1)
              ph = pm(iintmax+1)
              y0 = tlsave(iintmax)
              yl = tlsave(iintmax-1)
              yh = tlsave(iintmax+1)
              dy = yh-y0
              pbestint = p0 + 0.125*(ph-pl)*(ph-pl)/(2.0*p0-ph-pl)
              ybestint = y0 + 0.5*dy*(ph-pl)/(2.0*p0-ph-pl)
c              print*,'taus ',yl,y0,yh,' ls ',pl,p0,ph,' out ',ybestint,pbestint
              dely     = (ybestint-y0)/dy
              x0       = xm(iintmax)
              xl       = xm(iintmax-1)
              xh       = xm(iintmax+1)
              xbestint = x0 + 0.5*dely*((xh-xl)+dely*(xh+xl-2.0*x0))
c              print*,' sigs ',xl,x0,xh,' out ',xbestint
c              print*,'input ',pbest,tbest,sbest
c              print*,'output',pbestint,ybestint,xbestint
              if (pbestint.lt.pbest) then
                print*,'interpolation problem -- output max like ',pbestint,' lower than input ',pbest
                print*,' this is usually a sign that the grid is poorly centered '
                iedge = 2
              else 
                print*,' interp input ',pbest,tbest,sbest
                pbest = pbestint
                tbest = ybestint
                sbest = xbestint  
                print*,' interp output ',pbest,tbest,sbest
                endif
              endif
            endif
          

c Bayesian errorbars and median
          avgmag = avgin + abest
          do i=1,nt
            ptl(i) = 0.0
            enddo
          do i=1,ns
            psl(i) = 0.0
            enddo
          do i=1,nt
            do j=1,ns
              arg = pbest-psave(i,j)
              arg = min(arg,72.0)
              psave(i,j) = exp(-arg)
              ptl(i)     = ptl(i) + psave(i,j)
              psl(j)     = psl(j) + psave(i,j)
              enddo
            enddo
          call getlimit(tlsave,ptl,nt,tlow,thig,tmedian)
          call getlimit(slsave,psl,ns,slow,shig,smedian)


          print*,'best = ',pbest,' chi = ',cbest,' chi/dof = ',cbest/(npt-1)
          print*,'maxlike sigma = ',sbest,' tau = ',tbest 
          print*,'median  sigma = ',smedian,' tau = ',tmedian

          cdof = cbest/float(npt-3)
          write(iout,'(i2,12(1x,f13.6),1x,a65)')iedge,avgmag,pbest,cbest,cdof,sbest,slow,shig,tbest,tlow,thig,tmedian,smedian,fname

           if (igrid.eq.1) then
             print*,'after interpolation '
             print*,'pbest ',pbest,' chibest ',cbest,' for ',npt
             print*,'tbest ',tbest,' sbest ',sbest
c             stop
             endif

          return
          end

       program BayesianNMA
c      Yeongjin Gwon, UNMC
c      May, 2021
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (iprint=3,maxlag=20,nthin=5,nrep=1000)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n)
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 x(nx,n),z(nz,n)
       real*8 xmean(nx),x2mean(nx),xsd(nx)
       real*8 zmean(nz),z2mean(nz),zsd(nz)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ylatent(n),LU(2,n)
       real*8 ratioP(n),ratioP70(n),ratioP100(n)
       real*8 aBL(n,4),aBU(n,4),maxBL(n),minBU(n)
       real*8 aLmuk(K),aUmuk(K)
       real*8 temp1,temp2,temp3,temp4,temp5
       real*8 cnt,cnt1,cnt2,cnt3,cnt4,tempcnt
       integer ids(n),nStudy(n),trt(n),collap(n)
       integer icount1,icount2,icount3
       integer icount,nthin,burn
c
       real*8 seqphi(nz+1,nrep)
       real*8 seqdelta(nT,nrep)
       real*8 seqbeta(nx+1,nrep)
       real*8 seqalpha2(nx+1,nrep)
       real*8 seqalpha3(nx+1,nrep)
       real*8 seqmu(K,nrep)
       real*8 seqy(n,nrep)
       real*8 seqy1(n,nrep),seqy2(n,nrep)
       real*8 seqy3(n,nrep),seqy4(n,nrep)
       real*8 seqp70(n,nrep)
       real*8 seqp100(n,nrep)
       real*8 seqpratio(n,nrep)
c
       real*8 ac(maxlag+1),acv(maxlag+1),seac(maxlag+1),
     1        xacmean
       real*8 alow(2),aupp(2)
       real*8 ahpd(nrep),conf
       real*8 stime(500)
       integer imean,iseopt
c......phi0 - phi1
       real*8 sphi(nz+1),s2phi(nz+1)
       real*8 ephi(nz+1),sephi(nz+1)
       real*8 philow(nz+1),phiupp(nz+1)
c......deta1 - delta13
       real*8 sdelta(nT),s2delta(nT)
       real*8 edelta(nT),sedelta(nT)
       real*8 deltalow(nT),deltaupp(nT)
c......beta0 - beta5
       real*8 sbeta(nx+1),s2beta(nx+1)
       real*8 ebeta(nx+1),sebeta(nx+1)
       real*8 betalow(nx+1),betaupp(nx+1)
c......alpha20 - alpha25
       real*8 salpha2(nx+1),s2alpha2(nx+1)
       real*8 ealpha2(nx+1),sealpha2(nx+1)
       real*8 alpha2low(nx+1),alpha2upp(nx+1)
c......alph30 - alpha35
       real*8 salpha3(nx+1),s2alpha3(nx+1)
       real*8 ealpha3(nx+1),sealpha3(nx+1)
       real*8 alpha3low(nx+1),alpha3upp(nx+1)
c......mu1 - mu17
       real*8 ssmu(K),s2mu(K)
       real*8 emu(K),semu(K)
c......DIC computation
       real*8 stheta(K),etheta(K)
       real*8 barDIC,DICbar,pD,DIC
       real*8 res1,res2
       real*8 comp1,comp2,sumK,sumT
       real*8 tempDIC1,tempDIC2
       real*8 tbarDIC(n),tDICbar(n),tpD(n),tDIC(n)
       real*8 tsumT(n)
c......WAIC computation
       real*8 barWAIC,wpD,WAIC
       real*8 tbarWAIC(n),tbarWAIC2(n),twpD(n),tWAIC(n)
c......Prior parameters
       integer now(3)
       integer nowb(3),nowe(3)
       integer ntoday(3),ntodayb(3)
       real etime
       real elapsed(2)
       real total
c
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecylatent/ylatent       
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /vecratioP/ratioP
       common /vecratioP70/ratioP70
       common /vecratioP100/ratioP100
c 
       common /vecphi/phi
       common /vecdelta/delta
       common /vecmu/mu 
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vectheta/theta
       common /vecb/b
       common /vecaBL/aBL
       common /vecaBU/aBU
       common /vecLU/LU
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk  
       common /dummy/index
c
       common /vecephi/ephi
       common /vecedelta/edelta
       common /vecemu/emu 
       common /vecebeta/ebeta
       common /vecealpha2/ealpha2
       common /vecealpha3/ealpha3
       common /vecetheta/etheta
       common /vecseqphi/seqphi
       common /vecseqdelta/seqdelta
       common /vecseqmu/seqmu 
       common /vecseqbeta/seqbeta
       common /vecseqalpha2/seqalpha2
       common /vecseqalpha3/seqalpha3
       common /vecseqy/seqy
       common /vecseqy1/seqy1
       common /vecseqy2/seqy2
       common /vecseqy3/seqy3
       common /vecseqy4/seqy4
       common /vecseqp70/seqp70
       common /vecseqp100/seqp100
       common /vecseqpratio/pratio
       common /vecymode/ymode
       common /vecy1mode/y1mode
       common /vecy2mode/y2mode
       common /vecy3mode/y3mode
       common /v3cy4mode/y4mode

c......Start the running time 
       call idate(ntoday)
       ntodayb(1)=ntoday(1)
       ntodayb(2)=ntoday(2)
       ntodayb(3)=ntoday(3)
       call itime(now)
       nowb(1)=now(1)
       nowb(2)=now(2)
       nowb(3)=now(3)

       do j=1,2
        do i=1,n
         LU(j,i)=0.0d0
        enddo
       enddo

c......Read the data
       open(unit=24,file='CrohnRawRead.out',
     1    access='sequential',status='unknown')
       open(unit=15,file='CrohnRawSecondary-2.txt',status='old')
       do i=1,n
       read(15,*) ids(i),nStudy(i),trt(i),collap(i),
     1          ss(i),n70(i),n100(i),nrem(i),(x(j,i),j=1,nx)
       write(24,700) ids(i),nStudy(i),trt(i),collap(i),
     1          ss(i),n70(i),n100(i),nrem(i),(x(j,i),j=1,nx)
  700  format(I2,1x,I2,1x,I2,1x,I2,1x,9f10.2)
       enddo
       close(15)
       close(24)

c......set z variable for phi
       temp1=0.0d0
       temp2=0.0d0
       temp3=0.0d0
       temp4=0.0d0
       cnt=0.0d0

c......Put all covariates
       do j1=1,nz
        do i1=1,n
         z(j1,i1)=0.0d0
        enddo
       enddo
c
       do i2=1,n
        cnt=x(2,i2)
        x(2,i2)=cnt/ss(i2)
       enddo
c
       do k1=1,K 
        temp1=(x(4,2*k1-1)+x(4,2*k1))/2.0d0
        z(1,2*k1-1)=temp1
        z(1,2*k1)=z(1,2*k1-1)
        temp2=(x(1,2*k1-1)+x(1,2*k1))/2.0d0
        z(2,2*k1-1)=temp2
        z(2,2*k1)=z(2,2*k1-1)
        temp3=(x(2,2*k1-1)+x(2,2*k1))/2.0d0
        z(3,2*k1-1)=temp3
        z(3,2*k1)=z(3,2*k1-1)
        temp4=(x(3,2*k1-1)+x(3,2*k1))/2.0d0
        z(4,2*k1-1)=temp4
        z(4,2*k1)=z(4,2*k1-1)
        temp5=(x(5,2*k1-1)+x(5,2*k1))/2.0d0
        z(5,2*k1-1)=temp5
        z(5,2*k1)=z(5,2*k1-1)
       enddo
c
c......Standardize trial-level covariates
       do j1=1,nz
        zmean(j1)=0.0d0
        z2mean(j1)=0.0d0
        do i1=1,n
         zmean(j1)=zmean(j1)+z(j1,i1)
         z2mean(j1)=z2mean(j1)+z(j1,i1)**2
        enddo
        zmean(j1)=zmean(j1)/real(n)
        zsd(j1)=dsqrt((z2mean(j1)-
     1     real(n)*zmean(j1)**2)/real(n-1))
       enddo
c
       do j2=1,nz
        do i2=1,n
         z(j2,i2)=(z(j2,i2)-zmean(j2))/zsd(j2)
        enddo
       enddo

c......Standardize arm-level covariates
       do j3=1,nx
        xmean(j3)=0.0d0
        x2mean(j3)=0.0d0
        do i3=1,n
         xmean(j3)=xmean(j3)+x(j3,i3)
         x2mean(j3)=x2mean(j3)+x(j3,i3)**2
        enddo
        xmean(j3)=xmean(j3)/real(n)
        xsd(j3)=dsqrt((x2mean(j3)-
     1     real(n)*xmean(j3)**2)/real(n-1))
       enddo
c
       do j4=1,nx
        do i4=1,n
         x(j4,i4)=(x(j4,i4)-xmean(j4))/xsd(j4)
        enddo
       enddo

c......Set initialization parameter values
       do j1=1,nz+1
        phi(j1)=0.0d0
        sphi(j1)=0.0d0
        s2phi(j1)=0.0d0
        ephi(j1)=0.0d0
        sephi(j1)=0.0d0
        philow(j1)=0.0d0
        phiupp(j1)=0.0d0
       enddo
       do j2=1,nT
        delta(j2)=0.0d0
        sdelta(j2)=0.0d0
        s2delta(j2)=0.0d0
        edelta(j2)=0.0d0
        sedelta(j2)=0.0d0
        deltalow(j2)=0.0d0
        deltaupp(j2)=0.0d0
       enddo
       do j3=1,K
        mu(j3)=0.0d0
        aLmuk(j3)=0.0d0
        aUmuk(j3)=0.0d0
        ssmu(j3)=0.0d0
        s2mu(j3)=0.0d0
        emu(j3)=0.0d0
        semu(j3)=0.0d0
       enddo
       do j4=1,nx+1
        beta(j4)=0.0d0
        sbeta(j4)=0.0d0
        s2beta(j4)=0.0d0
        ebeta(j4)=0.0d0
        sebeta(j4)=0.0d0
        betalow(j4)=0.0d0
        betaupp(j4)=0.0d0
       enddo
       do j5=1,nx+1
        alpha2(j5)=0.0d0
        salpha2(j5)=0.0d0
        s2alpha2(j5)=0.0d0
        ealpha2(j5)=0.0d0
        sealpha2(j5)=0.0d0
        alpha2low(j5)=0.0d0
        alpha2upp(j5)=0.0d0
        alpha3(j5)=0.0d0
        salpha3(j5)=0.0d0
        s2alpha3(j5)=0.0d0
        ealpha3(j5)=0.0d0
        sealpha3(j5)=0.0d0
        alpha3low(j5)=0.0d0
        alpha3upp(j5)=0.0d0
       enddo
       do j6=1,n
        maxBL(j6)=0.0d0
        minBU(j6)=0.0d0
        do j7=1,4
         b(j6,j7)=0.0d0
         aBL(j6,j7)=0.0d0
         aBU(j6,j7)=0.0d0
        enddo
       enddo
       do j9=1,n
        y1(j9)=0.0d0
        y2(j9)=0.0d0
        y3(j9)=0.0d0
        y4(j9)=0.0d0
        ylatent(j9)=0.0d0
       enddo

c......Compute lower and upper bounds for the latent counts
       temp1=0.0d0
       temp2=0.0d0
       temp3=0.0d0
       temp4=0.0d0
       tempcnt=0.0d0

c......Code bounds 
       do i=1,n
        if (collap(i) .eq. 0) then
         temp1=(n70(i)-n100(i))-nrem(i)
         if (temp1 .le. 0.0d0) then               
          LU(1,i)=nrem(i)-(n70(i)-n100(i))
          LU(2,i)=nrem(i)
         else
          LU(1,i)=0.0d0
          LU(2,i)=nrem(i)
         endif
        else if (collap(i) .eq. -1) then
         temp2=(ss(i)-n100(i))-nrem(i)
         if (temp2 .le. 0.0d0) then
          LU(1,i)=nrem(i)-(ss(i)-n100(i))
          LU(2,i)=nrem(i)
         else 
          LU(1,i)=0.0d0
          LU(2,i)=nrem(i)
         endif
        endif
       enddo

c......Set initial values of the parameters       
       phi(1)=-6.511
       phi(2)=-0.501
       beta(1)=0.650
       beta(2)=-0.073
       beta(3)=0.248
       beta(4)=0.090
       beta(5)=-0.044
       beta(6)=-0.173
       delta(1)=-0.443
       delta(2)=-0.401
       delta(3)=-2.162
       delta(4)=-0.718
       delta(5)=-0.752
       delta(6)=-1.002
       delta(7)=-1.233
       delta(8)=-1.532
       delta(9)=-1.225
       delta(10)=-1.135
       delta(11)=-0.895
       delta(12)=-0.950
       delta(13)=-1.030
       do i=1,K
        call rnset( iseed )       
        mu(i)=DRNNOF()
        call rnget( iseed )        
       enddo       
       alpha2(1)=-4.196
       alpha2(2)=-0.090
       alpha2(3)=-0.123
       alpha2(4)=-0.139
       alpha2(5)=-0.075
       alpha2(6)=0.137
       alpha3(1)=-0.422
       alpha3(2)=-0.119
       alpha3(3)=0.078
       alpha3(4)=0.420
       alpha3(5)=-0.085
       alpha3(6)=-0.069

c......Warm up Gibbs
       burn=5000
       iseed=987412487

       call rnset( iseed )
       icount=0
       do i1=1,burn
        call gibbs( iseed )
         write(*,*) '-------------------------------------'
         write(*,*) '         Warm up Gibbs i =',i1
         write(*,*) '-------------------------------------'
         write(*,900) (beta(j),j=1,nx+1)
 900     format('beta  =',10f9.3)
         write(*,901) (delta(j),j=1,nT)
 901     format('delta =',10f9.3)
         write(*,903) (alpha2(j),j=1,nx+1)
 903     format('alpha2=',10f9.3)
         write(*,904) (alpha3(j),j=1,nx+1)
 904     format('alpha3=',10f9.3)
         write(*,905) (mu(j),j=1,K)
 905     format('mu    =',10f9.3) 
         write (*,906) (phi(j),j=1,2)
 906     format('phi   =',10f9.3)
         write(*,912) (aLmuk(j),j=1,K)
 912     format('Lower muk=',10f10.3)
         write(*,913) (aUmuk(j),j=1,K)
 913     format('Upper muk=',10f10.3)
         write(*,914) (ylatent(j),j=1,n)
 914     format('Latent Count=',18f6.0)
         write(*,915) (LU(1,j),j=1,n)
 915     format('Lower Count =',18f6.0)
         write(*,916) (LU(2,j),j=1,n)
 916     format('Upper Count =',18f6.0)
         icount=icount+1
       enddo         
       call rnget( iseed )

c......Initializing variables to calculate barDIC	    
       barDIC=0.0d0
       barWAIC=0.0d0
       res1=0.0d0
       sumT=0.0d0
       do k1=1,n
        tbarDIC(k1)=0.0d0
        tbarWAIC(k1)=0.0d0
        tbarWAIC2(k1)=0.0d0
       enddo

c......Running gibbs sampling
       icount=0
       call rnset ( iseed )
        do 220 i2=1,nrep
         do j0=1,nthin
          call gibbs ( iseed )
         enddo

c......Print out estimated parameters
         write(*,*) '======================================='
         write(*,*) '         Gibbs iteration =',i2
         write(*,*) '======================================='
         write(*,900) (beta(j),j=1,nx+1)
         write(*,901) (delta(j),j=1,nT)
         write(*,903) (alpha2(j),j=1,nx+1)
         write(*,904) (alpha3(j),j=1,nx+1)
         write(*,905) (mu(j),j=1,K)
         write(*,906) (phi(j),j=1,2)
         write(*,914) (ylatent(j),j=1,n)
         write(*,915) (LU(1,j),j=1,n)
         write(*,916) (LU(2,j),j=1,n)         
         icount=icount+1
c
         do 230 j1=1,nz+1
          sphi(j1)=sphi(j1)+phi(j1)
          s2phi(j1)=s2phi(j1)+phi(j1)**2
          seqphi(j1,i2)=phi(j1)
 230     continue
c
         do 231 j2=1,nT
          sdelta(j2)=sdelta(j2)+delta(j2)
          s2delta(j2)=s2delta(j2)+delta(j2)**2
          seqdelta(j2,i2)=delta(j2)          
 231     continue
c          
         do 232 j3=1,nx+1
          sbeta(j3)=sbeta(j3)+beta(j3)
          s2beta(j3)=s2beta(j3)+beta(j3)**2
          seqbeta(j3,i2)=beta(j3)          
 232     continue
c          
         do 233 j4=1,nx+1
          salpha2(j4)=salpha2(j4)+alpha2(j4)
          s2alpha2(j4)=s2alpha2(j4)+alpha2(j4)**2
          seqalpha2(j4,i2)=alpha2(j4)
          salpha3(j4)=salpha3(j4)+alpha3(j4)
          s2alpha3(j4)=s2alpha3(j4)+alpha3(j4)**2
          seqalpha3(j4,i2)=alpha3(j4)
 233     continue
c          
         do 234 j5=1,K
          ssmu(j5)=ssmu(j5)+mu(j5)
          s2mu(j5)=s2mu(j5)+mu(j5)**2
          seqmu(j5,i2)=mu(j5)
 234     continue
c          
         do 236 j7=1,n
c          seqy1(j7,i2)=y1(j7)
c          seqy2(j7,i2)=y2(j7)
c          seqy3(j7,i2)=y3(j7)
c          seqy4(j7,i2)=y4(j7)
          seqy(j7,i2)=ylatent(j7)
 236     continue
c         do j8=1,n
c          seqp100(j8,i2)=ratioP100(j8)
c          seqpratio(j8,i2)=ratioP(j8)
c         enddo
 
c......Compute barDIC
         do k1=1,n
          tempDIC=0.0d0
          index=k1
          call DIC1(index,res1)
          tempDIC=res1
          barDIC=barDIC+tempDIC
          tbarDIC(k1)=tbarDIC(k1)+tempDIC
c         write(*,*) 'n=',k1,'tempDIC=',tempDIC
c
c......Compute arm-level barWAIC
          tempWAIC=0.0d0
          tempWAIC=dexp(-0.5d0*tempDIC)
          tbarWAIC(k1)=tbarWAIC(k1)+tempWAIC
c         write(*,*) 'n=',k1,'tempWAIC=',tbarWAIC(k1)
         enddo
 220   continue

c......End Gibbs 
       call rnget( iseed ) 
       do k1=1,n
        tbarDIC(k1)=tbarDIC(k1)/real(nrep)
        tbarWAIC(k1)=tbarWAIC(k1)/real(nrep)
        barWAIC=barWAIC+dlog(tbarWAIC(k1))
        tbarWAIC2(k1)=-2.0d0*dlog(tbarWAIC(k1))
       enddo
       barWAIC=-2.0d0*barWAIC
       barDIC=barDIC/real(nrep)

c......Obtain posterior estimates
c......Posterior mean and std
       do j1=1,2
        ephi(j1)=sphi(j1)/real(nrep)
        sephi(j1)=dsqrt((s2phi(j1)-real(nrep)*
     1     ephi(j1)**2)/real(nrep-1))
       enddo
c       
       do j2=1,nT
        edelta(j2)=sdelta(j2)/real(nrep)
        sedelta(j2)=dsqrt((s2delta(j2)-real(nrep)*
     1     edelta(j2)**2)/real(nrep-1))
       enddo
c       
       do j3=1,nx+1
        ebeta(j3)=sbeta(j3)/real(nrep)
        sebeta(j3)=dsqrt((s2beta(j3)-real(nrep)*
     1      ebeta(j3)**2)/real(nrep-1))
       enddo
c       
       do j4=1,nx+1
        ealpha2(j4)=salpha2(j4)/real(nrep)
        sealpha2(j4)=dsqrt((s2alpha2(j4)-real(nrep)*
     1      ealpha2(j4)**2)/real(nrep-1))
        ealpha3(j4)=salpha3(j4)/real(nrep)
        sealpha3(j4)=dsqrt((s2alpha3(j4)-real(nrep)*
     1      ealpha3(j4)**2)/real(nrep-1))
       enddo
c       
       do j5=1,K
        emu(j5)=ssmu(j5)/real(nrep)
        semu(j5)=dsqrt((s2mu(j5)-real(nrep)*
     1      emu(j5)**2)/real(nrep-1))
       enddo
c       

c......Compute DICbar and complexity of the model
        DICbar=0.0d0
        pD=0.0d0
        DIC=0.0d0
        wpD=0.0d0
        WAIC=0.0d0
        sumT=0.0d0
        res2=0.0d0
        comp1=0.0d0
        comp2=0.0d0

        do k2=1,n
         tsumT(k2)=0.0d0
         tDICbar(k2)=0.0d0
         tpD(k2)=0.0d0
         tDIC(k2)=0.0d0
         twpD(k2)=0.0d0
         tWAIC(k2)=0.0d0
        enddo

c......Compute arm-level DIC
        do k2=1,n
         tempDIC=0.0d0
         index=k2
         call DIC2(index,res2)
         tempDIC=res2
         sumT=sumT+tempDIC
         tsumT(k2)=tempDIC
         tDICbar(k2)=tsumT(k2)
c         write(*,*) 'Arm=',k2,'barDIC=',tbarDIC(k2),
c     1       'DICbar=',tDICbar(k2)
        enddo

        DICbar=sumT
        pD=barDIC-DICbar
        DIC=DICbar+2.0d0*pD
        do k2=1,n
         tpD(k2)=tbarDIC(k2)-tDICbar(k2)
         tDIC(k2)=tDICbar(k2)+2.0d0*tpD(k2)
        enddo
c
        comp1=-0.5d0*barWAIC
        comp2=-0.5d0*barDIC
        wpD=comp1-comp2
        WAIC=barWAIC+2.0d0*wpD
        do k3=1,n
         twpD(k3)=-0.5d0*tbarWAIC2(k3)+0.5d0*tbarDIC(k3)
         tWAIC(k3)=tbarWAIC2(k3)+2.0d0*twpD(k3)
        enddo
c

c......Obtain 95% HPD intervals
        imean=1.0d0
        isepot=1.0d0
        alphahpd=0.05d0

c......HPD intervals for phi0 and phi1
        do 300 jmac=1,2
          write(*,*) 'HPD for phi=',jmac
          do 301 i=1,nrep
            ahpd(i)=seqphi(jmac,i)
 301      continue
          write(*,*) 'autocorrelation for phi',jmac
          call dacf(nrep,ahpd,iprint,isepot,imean,xacmean,
     1           maxlag,acv,ac,seac)
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          phiupp(jmac)=aupp(1)
          philow(jmac)=alow(1)
 300    continue
 
c......HPD intervals for betas
        do jmac=1,nx+1
          write(*,*) 'HPD for beta j=',jmac
          do i=1,nrep
            ahpd(i)=seqbeta(jmac,i)
          enddo
          write(*,*) 'autocorrelation for beta',jmac
          call dacf(nrep,ahpd,iprint,isepot,imean,xacmean,
     1         maxlag,acv,ac,seac)
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          betaupp(jmac)=aupp(1)
          betalow(jmac)=alow(1)
        enddo

c......HPD intervals for delta's
        do jmac=1,nT
          write(*,*) 'HPD for delta j=',jmac
          do i=1,nrep
            ahpd(i)=seqdelta(jmac,i)
          enddo
          write(*,*) 'autocorrelation for delta',jmac
          call dacf(nrep,ahpd,iprint,isepot,imean,xacmean,
     1       maxlag,acv,ac,seac)
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          deltaupp(jmac)=aupp(1)
          deltalow(jmac)=alow(1)
        enddo

c......HPD intervals for alpha2
        do jmac=1,nx+1
          write(*,*) 'HPD for alpha2 j=',jmac
          do i=1,nrep
            ahpd(i)=seqalpha2(jmac,i)
          enddo
          write(*,*) 'autocorrelation for alpha2',jmac
          call dacf(nrep,ahpd,iprint,isepot,imean,xacmean,
     1       maxlag,acv,ac,seac)
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          alpha2upp(jmac)=aupp(1)
          alpha2low(jmac)=alow(1)
        enddo

c......HPD intervals for alpha3
        do jmac=1,nx+1
          write(*,*) 'HPD for alpha3 j=',jmac
          do i=1,nrep
            ahpd(i)=seqalpha3(jmac,i)
          enddo
          write(*,*) 'autocorrelation for alpha3',jmac
          call dacf(nrep,ahpd,iprint,isepot,imean,xacmean,
     1        maxlag,acv,ac,seac)
          call hpd(nrep,alphahpd,ahpd,alow,aupp)
          alpha3upp(jmac)=aupp(1)
          alpha3low(jmac)=alow(1)
        enddo

c......Generate MCMC samples
       open(unit=25,file='R1-Logit1P100-P5-M3.txt',
     1        access='sequential',status='unknown')
       do j=1,nrep
       write(25,500) seqphi(1,j),seqphi(2,j),
     1   seqdelta(1,j),seqdelta(2,j),seqdelta(3,j),
     2   seqdelta(4,j),seqdelta(5,j),seqdelta(6,j),
     3   seqdelta(7,j),seqdelta(8,j),seqdelta(9,j),
     4   seqdelta(10,j),seqdelta(11,j),seqdelta(12,j),
     5   seqdelta(13,j),
     6   seqbeta(1,j),seqbeta(2,j),seqbeta(3,j),
     7   seqbeta(4,j),seqbeta(5,j),seqbeta(6,j)
 500   format(' ',21f10.3)
       enddo
       close(25)
c
       open(unit=16,file='R1-Logit2P100-P5-M3.txt',
     1        access='sequential',status='unknown')
       do j=1,nrep
       write(16,600) seqalpha2(1,j),seqalpha2(2,j),seqalpha2(3,j),
     1  seqalpha2(4,j),seqalpha2(5,j),seqalpha2(6,j),
     2  seqalpha3(1,j),seqalpha3(2,j),seqalpha3(3,j),
     3  seqalpha3(4,j),seqalpha3(5,j),seqalpha3(6,j)
 600   format(' ',12f10.4)
       enddo
       close(16)
c
       open(unit=19,file='R1-Logit3P100-P5-M3.txt',
     1     access='sequential',status='unknown')
       do j=1,nrep
       write(19,650) seqmu(1,j),seqmu(2,j),seqmu(3,j),
     1   seqmu(4,j),seqmu(5,j),seqmu(6,j),seqmu(7,j),
     2   seqmu(8,j),seqmu(9,j),seqmu(10,j),seqmu(11,j),
     3   seqmu(12,j),seqmu(13,j),seqmu(14,j),
     4   seqmu(15,j),seqmu(16,j),seqmu(17,j),seqmu(18,j)
 650   format(' ',18f10.4)
       enddo
       close(19)
c       
c......Generate latent counts
       open(unit=20,file='R1-LatentCountP100_P5.txt',
     1    access='sequential',status='unknown')
       do j=1,nrep
        write(20,550) seqy(1,j),seqy(2,j),seqy(3,j),seqy(4,j),
     1  seqy(5,j),seqy(6,j),seqy(7,j),seqy(8,j),seqy(9,j),
     2  seqy(10,j),seqy(11,j),seqy(12,j),seqy(13,j),seqy(14,j),
     3  seqy(15,j),seqy(16,j),seqy(17,j),seqy(18,j),seqy(19,j),
     4  seqy(20,j),seqy(21,j),seqy(22,j),seqy(23,j),seqy(24,j),
     5  seqy(25,j),seqy(26,j),seqy(27,j),seqy(28,j),seqy(29,j),
     6  seqy(30,j),seqy(31,j),seqy(32,j),seqy(33,j),seqy(34,j),
     7  seqy(35,j),seqy(36,j)
 550    format(' ',36f6.1)
       enddo
       close(20)

c......Output for trial level DIC	   
       open(unit=19,file='R1_DIC_P100_P5_M3.out',
     1        access='sequential',status='unknown')
       write(19,*) '================================================='
       write(19,*) '    Decomposition to the trial level DIC'
       write(19,*) '================================================='
       write(19,*) ' K     barDIC       pD     DICbar      DIC       '
       do j1=1,n
       write(19,810) j1,tbarDIC(j1),tpD(j1),tDICbar(j1),tDIC(j1)
 810   format(1x,I2,1x,4f10.3)
       enddo
       write(19,*) '================================================='
       close(19)

c......Output for trial level WAIC   
       open(unit=18,file='R1_WAIC_P100_P5_M3.out',
     1        access='sequential',status='unknown')
       write(18,*) '================================================='
       write(18,*) '    Decomposition to the trial level WAIC'
       write(18,*) '================================================='
       write(18,*) ' K     barWAIC       pD     barDIC     WAIC      '
       do j1=1,n
       write(18,810) j1,tbarWAIC2(j1),twpD(j1),tDICbar(j1),tWAIC(j1)
 811   format(1x,I2,1x,4f10.3)
       enddo
       write(18,*) '================================================='
       close(18)

c......Output of posterior summaries	   
       open(unit=14,file='R1_OUTPUT_P100_P5_M3.out',
     1        access='sequential',status='unknown')
       write(14,*) 'Bayesian NMA model with trial level random effect'
       write(14,*) 'Relative treatment model'
       write(14,*) '--------------------------------------------------'
       write(14,*) '           Model Comparison using DIC'
       write(14,*) '--------------------------------------------------'
       write(14,*) 'barDIC=',barDIC,'   pD=',pD  
       write(14,*) 'DICbar=',DICbar,'   DIC=',DIC
       write(14,*) '--------------------------------------------------'
       write(14,*) '           Model Comparison using WAIC'
       write(14,*) '--------------------------------------------------'
       write(14,*) 'barWAIC=',barWAIC,'    pD=',wpD
       write(14,*) 'WAIC=',WAIC
       write(14,*) '--------------------------------------------------'
       write(14,*) 'number of Gibbs iterations =',nrep
       write(14,*) 'number of thins =',nthin
       write(14,*) 'number of burn-in sample =',burn
       write(14,*) '=================================================='
       write(14,*) '              Posterior Summary                   '
       write(14,*) '=================================================='
       write(14,*) '  j     ephi       sephi    philow      phiupp '
c       
       do j1=1,2
       write(14,999) j1,ephi(j1),sephi(j1),philow(j1),phiupp(j1)
 999   format(1x,I2,1x,4f10.3)
       enddo
       write(14,*) '--------------------------------------------------'
       write(14,*) '  j    edelta     sedelta  deltalow    deltaupp'
       do j2=1,nT
       write(14,1100) j2,edelta(j2),sedelta(j2),deltalow(j2),
     1              deltaupp(j2)
 1100  format(1x,I2,1x,4f10.3)
       enddo
       write(14,*) '--------------------------------------------------'
       write(14,*) '  j    ebeta      sebeta    betalow    betaupp '
       do j3=1,nx+1
       write(14,1200) j3,ebeta(j3),sebeta(j3),betalow(j3),betaupp(j3)
 1200  format(1x,I2,1x,4f10.3)
       enddo
       write(14,*) '--------------------------------------------------'
       write(14,*) '  j   ealpha2    sealpha2  alpha2low   alpha2upp'
       do j4=1,nx+1
       write(14,1300) j4,ealpha2(j4),sealpha2(j4),alpha2low(j4),
     1            alpha2upp(j4)
 1300  format(1x,I2,4f10.3)
       enddo
       write(14,*) '  j   ealpha3    sealpha3  alpha3low   alpha3upp'
       do j5=1,nx+1
       write(14,1300) j5,ealpha3(j5),sealpha3(j5),alpha3low(j5),
     1            alpha3upp(j5) 
       enddo
       write(14,*) '---------------------------------------------------'

       call idate(ntoday)
       call itime(now)
       write(14,*)
     1   'date: month/day/year and time: hour, minute, and second'
       write(14,*)
     1  'Begining at date=',ntodayb(2),'/',ntodayb(1),'/',ntodayb(3)
       write(14,*) 'time at',nowb(1),':',nowb(2),':',nowb(3)
       write(14,*)
     1   'Ending at date=',ntoday(2),'/',ntoday(1),'/',ntoday(3)
       write(14,*) 'time at',now(1),':',now(2),':',now(3)
       total = etime(elapsed)
       write(14, *) 'elapsed time in minutes'
       write(14,2216) total/60.0d0,elapsed(1)/60.0d0,elapsed(2)/60.0d0
 2216  format(1x,'end: total=',f12.4,' user=',f12.4,
     1         ' system=', f12.4)
       write(14,*) '--------------------------------------------------'
       close(14)

       stop
       end program BayesianNMA
c
       include 'gibbsp100.f'
       include 'optim1.f'
       include 'hpd.f'
       include 'tnorm.f'
       include 'gilks1.f'
c......DICs for P1 to P3 and P6       
c       include 'DICp100.f'
c......DICs for P4 to P5
       include 'DICp100-Latent2.f'

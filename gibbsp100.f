       subroutine gibbs( iseed )
c      Yeongjin Gwon, UNMC
c      January 19, 2021
       implicit real*8 (a-h,o-z) 
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)       
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 ratioP(n),ratioP70(n),ratioP100(n)
       real*8 ylatent(n)       
       real*8 aBL(n,4),aBU(n,4),maxBL(n),minBU(n)
       real*8 aLmuk(K),aUmuk(K)
c       
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecLU/LU       
       common /vecphi/phi
       common /vecdelta/delta
       common /vecmu/mu   
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecaBL/aBL
       common /vecaBU/aBU
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk
       common /vecylatent/ylatent 
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /vecraioP/ratioP
       common /vecratioP70/ratioP70
       common /vecratioP100/ratioP100

c......Generate ylatent
       call Gy( iseed )

c......Generate beta 
       call Gbeta( iseed )
c      write(*,*) 'Gbeta=',beta

c......Generate delta 
       call Gdelta( iseed )
c      write(*,*) 'Gdelta=',delta
   
c......Generate alpha2
       call Galpha2( iseed )
c      write(*,*) 'Galpha2=',alpha2
   
c......Generate alpha3
       call Galpha3( iseed )
c      write(*,*) 'Galpha3=',alpha3

c......Generate auxiliary variaaBLes	   
       call Gb( iseed )
c      write(*,*) 'Auxiliary=',b

c......Generate phi modeling variance
       call Gphi( iseed )
c      write(*,*) 'Gphi=',phi
 
c......Generate trial level random effects
       call GmuK( iseed )
c      write(*,*) 'GmuK=',mu

       return
       end

c......Generation of latent counts
       subroutine Gy( iseed )
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5,nmax=100)
       integer nStudy(n),trt(n),collap(n)
       real*8 y1(n),y2(n),y3(n),y4(n)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 ylatent(n),priorunif
       real*8 ratioP70(n),ratioP100(n),ratioP(n)
       real*8 ein1,ein2,ein3,c1,c2,c3
       real*8 tRTnt,muK,xbeta,zbeta,eta2,eta3
       real*8 rn1,rn2,ru1,ru2
       real*8 upp,bot,a1,b1,a2,b2,pw
       real*8 p12star,p2star,p3star
       real*8 p1,p2,p3,diff
       real*8 p70,p100,pt1,pt2,maxlf,maxlg
       real*8 priorprob1(nmax),priorprob2(nmax)
       real*8 numprob1(nmax),numprob2(nmax)
       real*8 lf(nmax),lg(nmax),pt12(nmax)
       real*8 temp1,temp2,temp3,temp4,tempcnt
       real*8 cnt1,cnt2,cntk1,cntk2
       real*8 sumC,sumC1,sumC2
       real*8 tempval,tempval1,tempval2
       real*8 sval1,sval2
       real*8 coef1(nmax),coef2(nmax)
       integer sumy,iseed
c
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecLU/LU 
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu 
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /vecratioP/ratioP
       common /vecratioP70/ratioP70
       common /vecratioP100/ratioP100
       common /vecpriorprob1/priorprob1
       common /vecpriorprob2/priorprob2
       common /vecpt12/pt12
       external DBINDF,DBINPR

c......Initialization of variables   
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       p1=0.0d0
       p2=0.0d0
       p3=0.0d0
       ru1=0.0d0
       ru2=0.0d0
       priorunif=0.0d0
       diff=0.0d0
       tempval=0.0d0
       tempcnt=0.0d0

       do j=1,n
        y1(j)=0.0d0
        y2(j)=0.0d0
        y3(j)=0.0d0
        y4(j)=0.0d0
       enddo

c......Derivation of each count       
       do i=1,n
        if (collap(i) .eq. 0) then
         y1(i)=ss(i)-n70(i)
         y2(i)=n70(i)-n100(i)-(nrem(i)-ylatent(i))
         y3(i)=n100(i)-ylatent(i)
         y4(i)=nrem(i)
        else if (collap(i) .eq. -1) then
         y1(i)=0.0d0
         y2(i)=ss(i)-n100(i)-nrem(i)+ylatent(i)
         y3(i)=n100(i)-ylatent(i)
         y4(i)=nrem(i)
        else if (collap(i) .eq. 1) then
         y1(i)=ss(i)-n70(i)
         y2(i)=0.0d0
         y3(i)=n70(i)-nrem(i)
         y4(i)=nrem(i)
        else
         y1(i)=ss(i)-n70(i)
         y2(i)=0.0d0
         y3(i)=0.0d0
         y4(i)=n70(i)
        endif
       enddo

       do j=1,n
        ratioP(j)=0.0d0
        ratioP70(j)=0.0d0
        ratioP100(j)=0.0d0
       enddo

c......Construction of likelihood function
       do 100 i=1,n
        if (trt(i) .gt. 0.0d0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)
c
        rn1=0.0d0
        rn2=0.0d0
        p12star=0.0d0
        p2star=0.0d0
        p3star=0.0d0
        p70=0.0d0
        p100=0.0d0
        prior70=0.0d0
        prior100=0.0d0
        upp=0.0d0
        bot=0.0d0
        pt1=0.0d0
        pt2=0.0d0
        maxlf=0.0d0
        maxlg=0.0d0
        temp1=0.0d0
        temp2=0.0d0
        temp3=0.0d0
        temp4=0.0d0
        cnt1=0.0d0
        cnt2=0.0d0
        cnkt1=0.0d0
        cnkt2=0.0d0
c        
        do j2=1,nmax
         numprob1(j2)=0.0d0
         numprob2(j2)=0.0d0
         priorprob1(j2)=0.0d0
         priorprob2(j2)=0.0d0
         lf(j2)=0.0d0
         lg(j2)=0.0d0
         pt12(j2)=0.0d0
         coef1(j2)=0.0d0
         coef2(j2)=0.0d0
        enddo

c......Calculate cumulative probabilities
        if (ein1 .le. 0.0d0) then
         p1=1.0d0/(1.0d0+dexp(ein1))
        else
         p1=dexp(-ein1)/(1.0d0+dexp(-ein1))
        endif
        if (ein2 .le. 0.0d0) then
         p2=1.0d0/(1.0d0+dexp(ein2))
        else
         p2=dexp(-ein2)/(1.0d0+dexp(-ein2))
        endif
        if (ein3 .le. 0.0d0) then
         p3=1.0d0/(1.0d0+dexp(ein3))
        else
         p3=dexp(-ein3)/(1.0d0+dexp(-ein3))
        endif

c......Determine different power parameters
        a1=1.0d0
c       a1=0.5d0
c       a1=1.5d0
c        b1=1.0d0
c        b1=0.5d0
        b1=1.0d0
c......Power parameter needs to be explored        
        pw=1.0d0
c        pw=0.5d0
c        pw=1.5d0

       if (collap(i) .eq. 0) then
c......Conditional probability on the latent count from 100 response
        temp1=LU(1,i)
        temp2=LU(2,i)
        p3star=p3-p2
        if (p3star .lt. 0.001d0) p3star=0.001d0
        p2star=p2-p1
        if (p2star .lt. 0.001d0) p2star=0.001d0
        p100=p2star/p3star

c......Step1- Calculate numerator
c......The method itself up to the remission
        do j1=temp1,temp2
          sumC=0.0d0
          sumC1=0.0d0
          sumC2=0.0d0
          cntk1=0.0d0
          cntk2=0.0d0

c......Multinomial coefficient C1 in the paper
c......back to original version (coef1 and coef2)
          coef1(j1+1)=n70(i)-n100(i)-nrem(i)+j1
c         coef1(j1+1)=j1
          if (coef1(j1+1) .eq. 0.0d0) coef1(j1+1)=1.0d0
          coef2(j1+1)=n100(i)-j1
c         coef2(j1+1)=n70(i)-nrem(i)-j1
          if (coef2(j1+1) .eq. 0.0d0) coef2(j1+1)=1.0d0
          do k1=1,coef1(j1+1)
           cntk1=k1
           sumC1=sumC1+dlog(cntk1)
          enddo
          do k2=1,coef2(j1+1)
           cntk2=k2
           sumC2=sumC2+dlog(cntk2)
          enddo
          sumC=-sumC1-sumC2
c         write(*,*) 'i=',i,'C1=',sumC1,'C2=',sumC2,'SumC=',sumC

c......50% chance uncertainty was assigned to construct prior
          temp=floor((temp2+temp1)/2.0d0)
          tempval=0.0d0
c         temp=floor(temp1+(temp2-temp1+1.0d0)*0.30d0)
c         temp=temp2
          diff=(p3-p2)-(p2-p1)            
          if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
            temp=floor((temp1+temp2)/2.0d0)
            if (j1 .lt. temp) then
             priorprob1(j1+1)=(temp-j1)*dlog((p3-p1)/2.0d0)
            else if (j1 .gt. temp) then
             priorprob1(j1+1)=(j1-temp)*dlog((p3-p1)/2.0d0)
            else
             tempval1=0.0d0
             tempval2=0.0d0
             sval1=0.0d0
             sval2=0.0d0
             do k1=temp1,temp-1
              sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
              tempval1=tempval1+sval1
             enddo
             do k2=temp+1,temp2
              sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
              tempval2=tempval2+sval2
             enddo
             tempval=dexp(tempval1)+dexp(tempval2)
             if (tempval .gt. 1.0d0) tempval=0.999d0
             priorprob1(j1+1)=dlog(1.0d0-tempval)
            endif
          else
            if ((p2-p1) .lt. (p3-p2)) then
             num=(p3-p2)**b1
             den=(p2-p1)**a1+(p3-p2)**b1
             prior100=(num/den)**pw
             if (j1 .le. temp) then
c             priorprob1(j1+1)=(temp2-j1)*dlog(prior100)
              priorprob1(j1+1)=j1*dlog(prior100)
             else
              priorprob1(j1+1)=0.0d0
             endif
            else if ((p2-p1) .gt. (p3-p2)) then
             num=(p2-p1)**a1
             den=(p2-p1)**a1+(p3-p2)**b1
             prior100=(num/den)**pw
             if (j1 .gt. temp) then
c            priorprob1(j1+1)=j1*dlog(prior100)
              priorprob1(j1+1)=(temp2-j1)*dlog(prior100)
             else
              priorprob1(j1+1)=0.0d0
             endif
            endif
          endif
c
          priorunif=-dlog(temp2-temp1+1.0d0)
          numprob1(j1+1)=j1*dlog(p100)
          lf(j1+1)=numprob1(j1+1)
     1       +priorprob1(j1+1)
c     1       +priorunif
     2       +sumC
        enddo
        maxlf=lf(temp1+1)

c......Step2- Calculate denominator
        do j2=temp1,temp2
         if (maxlf .lt. lf(j2+1)) maxlf=lf(j2+1)
        enddo

c......Step3- Calculate probability
        do j3=temp1,temp2
         pt1=pt1+dexp(lf(j3+1)-maxlf)
        enddo
        do j4=temp1,temp2
         pt2=dexp(lf(j4+1)-maxlf)
         pt12(j4+1)=pt2/pt1
        enddo
c       
c        write(*,*) 'Study ID=',i
c        write(*,201) pt2
c 201    format('Prob Num=',1f8.3)
c        write(*,202) pt1
c 202    format('Prob Dem =',1f8.3)
c        write(*,203) pt12
c 203    format('All ratio Probs =',10f10.3)
c        write(*,*) 'i=',i
c        write(*,204) priorprob1
c 204    format('Prior prob= =',f16.4)

c......Step4- Generation of latent count
        rn1=0.0d0
        call rnset( iseed )
         rn1=DRNUNF()
        call rnget( iseed )
c       
        cnt1=LU(1,i)
 110    if (rn1 .le. pt12(cnt1+1)) then
         ylatent(i)=cnt1
        else
         rn1=rn1-pt12(cnt1+1)
         cnt1=cnt1+1.0d0
         goto 110
        endif
c......Keep track of ratio and other probabilities
        ratioP100(i)=p100
        ratioP(i)=pt12(cnt1+1)
c       write(*,*) 'i=',i
c       write(*,41) ratioP(i)
c41     format('Ratio Prob=',1f8.3)

c......Extreme case 1: all remission is from 100 response (P1)
c        ylatent(i)=LU(2,i)
c......Extreme case 2: all remission is from 70 response (P2)
c         ylatent(i)=LU(1,i)
c......Equal allocation ratio 50-50 (P3)
c        ylatent(i)=floor((temp1+temp2)/2.0d0)
c......Discrete uniform (P4) and proposed P5
         ylatent(i)=cnt1
c......P6 to combine P1 to P3
c        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
c         ylatent(i)=floor((temp1+temp2)/2.0d0)
c        else
c         if ((p2-p1) .lt. (p3-p2)) then
c          ylatent(i)=LU(1,i)
c         else if ((p2-p1) .gt. (p3-p2)) then
c          ylatent(i)=LU(2,i)
c         endif
c        endif
c	   
       else if (collap(i) .eq. -1) then
c......Conditional probability on the latent count from <=100 response
        temp3=LU(1,i)
        temp4=LU(2,i)
        p12star=p2
        if (p12star .lt. 0.001d0) p12star=0.001d0
        p3star=p3-p2
        if (p3star .lt. 0.001d0) p3star=0.001d0
        p100=p12star/p3star
c
c......Step1- Calculate numerator
        do j1=temp3,temp4
         sumC=0.0d0
         sumC1=0.0d0
         sumC2=0.0d0
         cntk1=0.0d0
         cntk2=0.0d0

c......Multinomial coefficient C2 in the paper
c         tempcnt=nrem(i)-(ss(i)-n100(i))
c         if (tempcnt .lt. 0.0d0) then
c           coef1(j1+1)=ss(i)-n100(i)-nrem(i)+j1
c           if (coef1(j1+1) .eq. 0.0d0) coef1(j1+1)=1.0d0
c           coef2(j1+1)=n100(i)-j1
c           if (coef2(j1+1) .eq. 0.0d0) coef2(j1+1)=1.0d0
c         else
c           coef1(j1+1)=j1
c           if (coef1(j1+1) .eq. 0.0d0) coef1(j1+1)=1.0d0
c           coef2(j1+1)=ss(i)-nrem(i)-j1
c           if (coef2(j1+1) .eq. 0.0d0) coef2(j1+1)=1.0d0
c         endif
c         
         coef1(j1+1)=ss(i)-n100(i)-nrem(i)+j1
         coef2(j1+1)=n100(i)-j1
         do k1=1,coef1(j1+1)
          cntk1=k1
          sumC1=sumC1+dlog(cntk1)
         enddo
         do k2=1,coef2(j1+1)
          cntk2=k2
          sumC2=sumC2+dlog(cntk2)
         enddo
         sumC=-sumC1-sumC2

c......Truncated Power series prior distribution
          temp=floor((temp4+temp3)/2.0d0)
c         temp=floor(temp3+(temp4-temp3+1.0d0)*0.30d0)
c         temp=temp4
          tempval=0.0d0
          tempval1=0.0d0
          tempval2=0.0d0
          sval1=0.0d0
          sval2=0.0d0
          diff=(p3-p2)-(p2-p1)
          if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
            temp=floor((temp3+temp4)/2.0d0)
            if (j1 .lt. temp) then
             priorprob2(j1+1)=(temp-j1)*dlog((p3-p1)/2.0d0)
            else if (j1 .gt. temp) then
             priorprob2(j1+1)=(j1-temp)*dlog((p3-p1)/2.0d0)
            else
             do k1=temp1,temp-1
              sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
              tempval1=tempval1+sval1
             enddo
             do k2=temp+1,temp2
              sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
              tempval2=tempval2+sval2
             enddo
             tempval=dexp(tempval1)+dexp(tempval2)
             if (tempval .ge. 1.0d0) tempval=0.999d0
             priorprob2(j1+1)=dlog(1.0d0-tempval)
            endif
          else
            if ((p2-p1) .lt. (p3-p2)) then
             num=(p3-p2)**b1
             den=(p2-p1)**a1+(p3-p2)**b1
             prior100=(num/den)**pw
             if (j1 .le. temp) then
c             priorprob2(j1+1)=(temp4-j1)*dlog(prior100)
              priorprob2(j1+1)=j1*dlog(prior100)
             else
              priorprob2(j1+1)=0.0d0
             endif
            else if ((p2-p1) .gt. (p3-p2)) then
             num=(p2-p1)**a1
             den=(p2-p1)**a1+(p3-p2)**b1
             prior100=(num/den)**pw
             if (j1 .gt. temp) then
c             priorprob2(j1+1)=j1*dlog(prior100)
              priorprob2(j1+1)=(temp4-j1)*dlog(prior100)
             else
              priorprob2(j1+1)=0.0d0
             endif
            endif
          endif
c
          priorunif=-dlog(temp4-temp3+1.0d0)
          numprob2(j1+1)=j1*dlog(p100)
          lg(j1+1)=numprob2(j1+1)
     1       +priorprob2(j1+1)
c     1       +priorunif
     2       +sumC
        enddo
        maxlg=lg(temp3+1)

c......Step2- Calculate denominator
        do j2=temp3,temp4
         if (maxlg .lt. lg(j2+1)) maxlg=lg(j2+1)
        enddo

c......Step3- Calculate probability
        do j3=temp3,temp4
         pt1=pt1+dexp(lg(j3+1)-maxlg)
        enddo
        do j4=temp3,temp4
         pt2=dexp(lg(j4+1)-maxlg)
         pt12(j4+1)=pt2/pt1
        enddo
c       write(*,*) 'i=',i,'Multi coeffi=',sumC

c......Step4- Generation of latent count
        rn2=0.0d0
        call rnset( iseed )
         rn2=DRNUNF()
        call rnget( iseed )
c
        cnt2=LU(1,i)
 120    if (rn2 .le. pt12(cnt2+1)) then
         ylatent(i)=cnt2
        else
         rn2=rn2-pt12(cnt2+1)
         cnt2=cnt2+1.0d0
         goto 120
        endif

c......Keep track of ratio and other probabilities
        ratioP100(i)=p100        
        ratioP(i)=pt12(cnt2+1)
c       write(*,*) 'i=',i
c        write(*,43) ratioP(i)
c 43     format('Ratio Prob=',1f8.3)

c......Extreme case 1: all remission is from 100 reponse (P1)
c        ylatent(i)=LU(2,i)
c......Extreme case 2: all remission is from 70 response (P2)
c        ylatent(i)=LU(1,i)
c......Equal allocation 50-50 (P3)
c        ylatent(i)=floor((temp3+temp4)/2.0d0)
c......Discrete uniform (P4) and proposed P5
         ylatent(i)=cnt2
c......P6 to combine P1 to P3
c        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
c         ylatent(i)=floor((temp3+temp4)/2.0d0)
c        else
c         if ((p2-p1) .lt. (p3-p2)) then
c          ylatent(i)=LU(1,i)
c         else if ((p2-p1) .gt. (p3-p2)) then
c          ylatent(i)=LU(2,i)
c         endif
c       endif
c
       endif
 100   continue
       return
       end


c......Generation of betas
       subroutine Gbeta( iseed )
c      Yeongjin Gwon,UNMC
c      Jun 19, 2019
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 beta(nx+1)
       parameter (maxtries=25,im=2,NS=10,np=7,in=1000)
       real*8 emax,step
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1          xlb,xub,rwv(6*NS+15+in)
       integer ifault,iwv(NS+7+in),ifcount
       integer iseed
       logical lb,ub
       common /vecbeta/beta
       common /dummy/idum1   
       external eval1
c  
       emax=60.0d0
       do 500 j=1,nx+1
         idum1=j
         call Optbeta
         step=1.0d0 
         a(1)=beta(j)-step
         a(2)=beta(j)+step
520      call eval1(a(1), ha(1), hpa(1))
         call eval1(a(2), ha(2), hpa(2))
c        write(*,*) 'a(1)=',a(1),' ha(1)=',ha(1)
c        write(*,*) ' hpa(1)=',hpa(1)
c ---------- Get starting points for Gilks ----------
         if (hpa(2) .gt. 0.0) then
c......Both points to the left of the mode; push a(2) to the right
530        a(2)=a(2)+step
           call eval1(a(2), ha(2), hpa(2))
c          if a(2) still to the left of the mode repeat, else continue
             if (hpa(2) .gt. 0.0) then
               GOTO 530
             else
               a(1)=a(2)-step
             endif
           else
           if (hpa(1) .le. 0.0) then
c......Both points to the right of the mode
c......Insert new point to the left of a(1)
             a(2)=a(1)
             ha(2)=ha(1)
             hpa(2)=hpa(1)
540          a(1)=a(1)-step
             call eval1(a(1), ha(1), hpa(1))
             if (hpa(1) .le. 0.0) then
               GOTO 540
c            else
c              a(2)=a(1)+step
             endif
           endif
         endif
c -----------------------------------------------------------------------
         ifcount=0
550      call INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, eval1)
         if ((ifault .eq. 3) .or. (ifault .eq. 4)) GOTO 520
c
         call SAMPLE(iwv,rwv,eval1,val,ifault,iseed)
         if (ifault .eq. 5) then
c          print *, 'ifault =', ifault, ' trying again'
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           if(ifcount .le. maxtries) then
             GOTO 550
           else
c          print *, 'Too many tries resetting starting values'
             a(1)=0.0
             a(2)=0.0
             step=1.2d0*step
             GOTO 520
           endif
         endif
         if (ifault .eq. 7) then
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           if (ifcount .le. maxtries) then
             GOTO 550
           else
c            print *, 'Too many tries resetting starting values'
             a(1)=0.0
             a(2)=0.0
             step=1.2d0*step
             GOTO 520
           endif
         endif
         if (IFAULT .NE. 0) then
           print *, ' The sampling status is ', IFAULT
         endif
         beta(j)=val
 500   continue
       return
       end 
 

       subroutine eval1(z1,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       integer nStudy(n),trt(n),collap(n)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n),ylatent(n)
       real*8 sumLL,sumDL
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a   
       real*8 dllik1,dllik2,dllik3,dllik4  
       real*8 dllik2a,dllik2b,dllik2c
       real*8 dllik3a,dllik3b,dllik3c,dllik4a
       real*8 ein1,ein2,ein3,c1,c2,c3
       real*8 tRTnt,muK,xbeta,zbeta,eta2,eta3
       real*8 sy12,sy23,sy234
       real*8 z1,ahz,ahpz
       real*8 dX(n,nx+1),xt(n,nx+1)
       real*8 XtX(nx+1,nx+1),INVXtX(nx+1,nx+1)       
       integer iseed,idum1
c	 
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu 
       common /dummy/idum1
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       external DLINDS       
       beta(idum1)=z1

c......Initialization of variables
       sumLL=0.0d0
       sumDL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       dllik1=0.0d0
       dllik2=0.0d0
       dllik3=0.0d0
       dllik4=0.0d0
       dllik2a=0.0d0
       dllik2b=0.0d0
       dllik2c=0.0d0
       dllik3a=0.0d0
       dllik3b=0.0d0
       dllik3c=0.0d0
       dllik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0

c......Construction of likelihood function
       do 590 i=1,n
        if (trt(i) .gt. 0.0d0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Logliklihood construction
c......Fully observed categories
        if (collap(i) .eq. 0) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
           if (idum1 .eq. 1) then
             dllik1=dllik1+
     1         (dexp(ein1)/(1.0d0+dexp(ein1)))*y1(i)
           else
             dllik1=dllik1+x(idum1-1,i)*
     1         (dexp(ein1)/(1.0d0+dexp(ein1)))*y1(i)
           endif
          else
           llik1=llik1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
           if (idum1 .eq. 1) then
            dllik1=dllik1+
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))*y1(i)
           else
            dllik1=dllik1+x(idum1-1,i)*
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))*y1(i)
           endif
          endif
         endif
c        write(*,*) 'i=',i,'llik1=',llik1,'dllik1=',dllik1
c
         if (y2(i) .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
           if (idum1 .eq. 1) then
             dllik2b=-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1))
           else
             dllik2b=x(idum1-1,i)*
     1         (-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1)))
           endif
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
           if (idum1 .eq. 1) then
             dllik2b=-dexp(-ein1)/(1.0d0+dexp(-ein1))
           else
             dllik2b=x(idum1-1,i)*
     1          (-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
          if (ein2 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein2))
           if (idum1 .eq. 1) then
            dllik2c=dexp(ein2)/(1.0d0+dexp(ein2))
           else
            dllik2c=x(idum1-1,i)*
     1        (dexp(ein2)/(1.0d0+dexp(ein2)))
           endif
          else
           llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
           if (idum1 .eq. 1) then
            dllik2c=1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2))
           else
            dllik2c=x(idum1-1,i)*
     1        (1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           endif
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
          dllik2=dllik2+(dllik2b+dllik2c)*y2(i)
c         write(*,*) 'i=',i,'llik2=',llik2,'dllik2=',dllik2	
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
           if (idum1 .eq. 1) then
            dllik3b=-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2))
           else
            dllik3b=x(idum1-1,i)*
     1         (-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2)))
           endif
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
           if (idum1 .eq. 1) then
            dllik3b=-dexp(-ein2)/(1.0d0+dexp(-ein2))
           else
            dllik3b=x(idum1-1,i)*
     1        (-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           endif
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
           if (idum1 .eq. 1) then
            dllik3c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            dllik3c=x(idum1-1,i)*
     1        (dexp(ein3)/(1.0d0+dexp(ein3)))
           endif
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
           if (idum1 .eq. 1) then
            dllik3c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           else
            dllik3c=x(idum1-1,i)*
     1        (1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
          dllik3=dllik3+(dllik3b+dllik3c)*y3(i)
c         write(*,*) 'i=',i,'dllik3=',dllik3			  
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1       (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           endif
          else
           llik4=llik4+y4(i)*
     1       (-dlog(1.0d0+dexp(-ein3)))
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1       (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1       (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
         endif
c		  write(*,*) 'i=',i,llik1,dllik1,llik2,
c     1          dllik2,llik3,dllik3

c......Partially observed 1 and 2 are collapsed
        else if (collap(i) .eq. -1) then
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein2))*sy12
           if (idum1 .eq. 1) then
            dllik1=dllik1+sy12*
     1        (dexp(ein2)/(1.0d0+dexp(ein2)))
           else
            dllik1=dllik1+sy12*x(idum1-1,i)*
     1        (dexp(ein2)/(1.0d0+dexp(ein2)))
           endif 
          else
           llik1=llik1+(-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
           if (idum1 .eq. 1) then
            dllik1=dllik1+sy12*
     1        (1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           else
            dllik1=dllik1+sy12*x(idum1-1,i)*
     1        (1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           endif
          endif
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
           if (idum1 .eq. 1) then
            dllik3b=-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2))
           else
            dllik3b=x(idum1-1,i)*
     2       (-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2)))
           endif
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
           if (idum1 .eq. 1) then
            dllik3b=-dexp(-ein2)/(1.0d0+dexp(-ein2))
           else
            dllik3b=x(idum1-1,i)*
     1         (-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           endif
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
           if (idum1 .eq. 1) then
            dllik3c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            dllik3c=x(idum1-1,i)*
     1        (dexp(ein3)/(1.0d0+dexp(ein3)))
           endif
          else
            llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
           if (idum1 .eq. 1) then
            dllik3c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           else
            dllik3c=x(idum1-1,i)*
     1        (1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
          dllik3=dllik3+(dllik3b+dllik3c)*y3(i)
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+y4(i)*
     1        (ein3-dlog(1.0d0+dexp(ein3)))
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1      (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           endif
          else
           llik4=llik4+(-dlog(1.0d0+dexp(-ein3)))*y4(i)
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1       (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(i) .eq. 1) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
           if (idum1 .eq. 1) then
            dllik1=dllik1+dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           else
            dllik1=dllik1+x(idum1-1,i)*
     1        dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           endif
          else
           llik1=llik1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
           if (idum1 .eq. 1) then
            dllik1=dllik1+y1(i)*
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           else
            dllik1=dllik1+y1(i)*x(idum1-1,i)*
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
         endif
c         
         sy23=y2(i)+y3(i)
         if (sy23 .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein3-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
           if (idum1 .eq. 1) then
            dllik2b=-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1))
           else
            dllik2b=x(idum1-1,i)*
     1        (-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1)))
           endif
          else
            llik2b=-dlog(1.0d0+dexp(-ein1))
           if (idum1 .eq. 1) then
            dllik2b=-dexp(-ein1)/(1.0d0+dexp(-ein1))
           else
            dllik2b=x(idum1-1,i)*
     1        (-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
          if (ein3 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein3))
           if (idum1 .eq. 1) then
            dllik2c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            dllik2c=x(idum1-1,i)*
     1        (dexp(ein3)/(1.0d0+dexp(ein3)))
           endif 
          else
           llik2c=-ein3-dlog(1.0d0+dexp(-ein3))
           if (idum1 .eq. 1) then
            dllik2c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           else
            dllik2c=x(idum1-1,i)*
     1       (1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*sy23
          dllik2=dllik2+(dllik2b+dllik2c)*sy23
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+y4(i)*
     1       (ein3-dlog(1.0d0+dexp(ein3)))
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           endif
          else
           llik4=llik4+(-dlog(1.0d0+dexp(-ein3)))*y4(i)
           if (idum1 .eq. 1) then
            dllik4=dllik4+y4(i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           else
            dllik4=dllik4+y4(i)*x(idum1-1,i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif
         endif

c......partially observed 2, 3, and 4 collapsed
        else if (collap(i) .eq. 2) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
           if (idum1 .eq. 1) then
            dllik1=dllik1+
     1        dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           else
            dllik1=dllik1+x(idum1-1,i)*
     1        dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           endif
          else
           llik1=llik1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
           if (idum1 .eq. 1) then
            dllik1=dllik1+y1(i)*
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           else
            dllik1=dllik1+y1(i)*x(idum1-1,i)*
     1        (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
         endif
c         
         sy234=y2(i)+y3(i)+y4(i)
         if (sy234 .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik4=llik4+sy234*
     1       (ein1-dlog(1.0d0+dexp(ein1)))
           if (idum1 .eq. 1) then
            dllik4=dllik4+sy234*
     1       (-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1)))
           else
            dllik4=dllik4+sy234*x(idum1-1,i)*
     1         (-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1)))
           endif
          else
           llik4=llik4+sy234*
     1       (-dlog(1.0d0+dexp(-ein1)))
           if (idum1 .eq. 1) then
            dllik4=dllik4+sy234*
     1       (-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           else
            dllik4=dllik4+sy234*x(idum1-1,i)*
     1         (-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
         endif
        endif
c      write(*,*) 'i=',i,'sumLL=',sumLL,'sumDL=',sumDL	   
 590   continue
c
       sumLL=llik1+llik2+llik3+llik4
       sumDL=dllik1+dllik2+dllik3+dllik4
c       write(*,*) 'sumLL=',sumLL,'sumDL=',sumDL
       ahz=sumLL-1.0d0/10000.0d0*z1**2
       ahpz=sumDL-1.0d0/5000.0d0*z1
       return
       end

c......Generation of deltas
       subroutine Gdelta( iseed )
c      Yeongjin Gwon, UNMC
c      Jun 19, 2019
       implicit real*8 (a-h,o-z)
       parameter(n=36,K=18,nT=13,nx=5,nz=5)
       real*8 delta(nT)
c      use Gilks log-concave sampling
       parameter (maxtries=25,im=2,NS=10,np=7,in=1000)
       real*8  emax,step
       real*8  val,a(im+in),ha(im+in),hpa(im+in),
     1          xlb,xub,rwv(6*NS+15+in)
       integer ifault,iwv(NS+7+in),ifcount
       integer iseed
       logical lb,ub
       common /vecdelta/delta
       common /dummy/idum2
       external eval2
c	 
       emax=60.0d0
       do 300 j=1,nT
         idum2=j
         call Optdelta
         step=1.0d0
         a(1)=delta(j)-step
         a(2)=delta(j)+step
320      call eval2(a(1), ha(1), hpa(1))
         call eval2(a(2), ha(2), hpa(2))
c        write(*,*) 'a(1)=',a(1),' ha(1)=',ha(1)
c        write(*,*) ' hpa(1)=',hpa(1)
c ------ Get starting points for Gilks ------
         if (hpa(2) .gt. 0.0) then
c......Both points to the left of the mode; push a(2) to the right
330        a(2)=a(2)+step
           call eval2(a(2), ha(2), hpa(2))
c          if a(2) still to the left of the mode repeat, else continue
           if (hpa(2).gt.0.0) then
             GOTO 330
           else
             a(1)=a(2)-step
           endif
         else
           if (hpa(1) .le. 0.0) then
c......Both points to the right of the mode
c......Insert new point to the left of a(1)
             a(2)=a(1)
             ha(2)=ha(1)
             hpa(2)=hpa(1)
340          a(1)=a(1)-step
             call eval2(a(1), ha(1), hpa(1))
             if (hpa(1) .le. 0.0) then
               GOTO 340
c            else
c             a(2)=a(1)+step
             endif
           endif
         endif
c -------------------------------------------
         ifcount=0
350      call INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, eval2)
         if ((ifault .eq. 3) .or. (ifault .eq. 4)) GOTO 320
c
         call SAMPLE(iwv,rwv,eval2,val,ifault,iseed)
         if (ifault .eq. 5) then
c          print *, 'ifault =', ifault, ' trying again'
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           if (ifcount .le. maxtries) then
             GOTO 350
           else
c            print *, 'Too many tries resetting starting values'
             a(1)=0.0
             a(2)=0.0
             step=1.2d0*step
             GOTO 320
           endif
         endif
         if (ifault .eq. 7) then
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           if (ifcount .le. maxtries) then
             GOTO 350
           else
c            print *, 'Too many tries resetting starting values'
             a(1)=0.0
             a(2)=0.0
             step=1.2d0*step
             GOTO 320
           endif
         endif
         if (IFAULT .ne. 0) then
c          print *, ' The sampling status is ', IFAULT
         endif
         delta(j)=val
 300   continue
       return
       end


       subroutine eval2(z2,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n),ylatent(n)
       real*8 sumLL,sumDL,sumC
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a   
       real*8 dllik1,dllik2,dllik3,dllik4  
       real*8 dllik2a,dllik2b,dllik2c
       real*8 dllik3a,dllik3b,dllik3c,dllik4a
       real*8 allik,adlik
       real*8 tRTnt,muK,xbeta,zbeta,eta2,eta3   
       real*8 c1,c2,c3,ein1,ein2,ein3
       real*8 sy12,sy23,sy234
       real*8 z2,ahz,ahpz
       integer nStudy(n),trt(n),collap(n)
       integer iseed
c	 
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu   
       common /dummy/idum2
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       delta(idum2)=z2

c......Initialization of variables
       sumC=0.0d0
       sumLL=0.0d0
       sumDL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       dllik1=0.0d0
       dllik2=0.0d0
       dllik3=0.0d0
       dllik4=0.0d0
       dllik2a=0.0d0
       dllik2b=0.0d0
       dllik2c=0.0d0
       dllik3a=0.0d0
       dllik3b=0.0d0
       dllik3c=0.0d0
       dllik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
c
       do 390 i=1,n
        if (trt(i) .gt. 0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

        if (tRTnt .eq. delta(idum2)) then
c......Fully observed categories		
         if (collap(i) .eq. 0) then
          if (y1(i) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
            dllik1=dllik1+
     1        (dexp(ein1)/(1.0d0+dexp(ein1)))*y1(i)
           else
            llik1=llik1+
     1         (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
            dllik1=dllik1+
     1         (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))*y1(i)
           endif
          endif
c
          if (y2(i) .gt. 0.0d0) then
           llik2a=dlog(1.0d0-dexp(ein2-ein1))
           if (ein1 .le. 0.0d0) then
            llik2b=ein1-dlog(1.0d0+dexp(ein1))
            dllik2b=-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1))
           else
            llik2b=-dlog(1.0d0+dexp(-ein1))
            dllik2b=-dexp(-ein1)/(1.0d0+dexp(-ein1))
           endif
           if (ein2 .le. 0.0d0) then
            llik2c=-dlog(1.0d0+dexp(ein2))
            dllik2c=dexp(ein2)/(1.0d0+dexp(ein2))
           else
            llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
            dllik2c=1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2))
           endif
           llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
           dllik2=dllik2+(dllik2b+dllik2c)*y2(i)
          endif
c          
          if (y3(i) .gt. 0.0d0) then
           llik3a=dlog(1.0d0-dexp(ein3-ein2))
           if (ein2 .le. 0.0d0) then
            llik3b=ein2-dlog(1.0d0+dexp(ein2))
            dllik3b=-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2))
           else
            llik3b=-dlog(1.0d0+dexp(-ein2))
            dllik3b=-dexp(-ein2)/(1.0d0+dexp(-ein2))
           endif
           if (ein3 .le. 0.0d0) then
            llik3c=-dlog(1.0d0+dexp(ein3))
            dllik3c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
            dllik3c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           endif
           llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
           dllik3=dllik3+(dllik3b+dllik3c)*y3(i)
          endif
c          
          if (y4(i) .gt. 0.0d0) then
           if (ein3 .lt. 0.0d0) then
            llik4=llik4+y4(i)*
     1        (ein3-dlog(1.0d0+dexp(ein3)))
            dllik4=dllik4+y4(i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            llik4=llik4+y4(i)*
     1        (-dlog(1.0d0+dexp(-ein3)))
            dllik4=dllik4+y4(i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif

c......Partially observed 1 and 2 are collapsed
         else if (collap(i) .eq. -1) then
          sy12=y1(i)+y2(i)
          if (sy12 .gt. 0.0d0) then
           if (ein2 .le. 0.0d0) then
            llik1=llik1-dlog(1.0d0+dexp(ein2))*sy12
            dllik1=dllik1+sy12*
     1          (dexp(ein2)/(1.0d0+dexp(ein2)))
           else
            llik1=llik1+
     1         (-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
            dllik1=dllik1+sy12*
     1         (1.0d0-dexp(-ein2)/(1.0d0+dexp(-ein2)))
           endif
          endif
c          
          if (y3(i) .gt. 0.0d0) then
           llik3a=dlog(1.0d0-dexp(ein3-ein2))
           if (ein2 .le. 0.0d0) then
            llik3b=ein2-dlog(1.0d0+dexp(ein2))
            dllik3b=-1.0d0+dexp(ein2)/(1.0d0+dexp(ein2))
           else
            llik3b=-dlog(1.0d0+dexp(-ein2))
            dllik3b=-dexp(-ein2)/(1.0d0+dexp(-ein2))             
           endif
           if (ein3 .le. 0.0d0) then
            llik3c=-dlog(1.0d0+dexp(ein3))
            dllik3c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
            dllik3c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           endif
           llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
           dllik3=dllik3+(dllik3b+dllik3c)*y3(i)
          endif
c
          if (y4(i) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4=llik4+y4(i)*
     1        (ein3-dlog(1.0d0+dexp(ein3)))
            dllik4=dllik4+y4(i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            llik4=llik4+y4(i)*
     1        (-dlog(1.0d0+dexp(-ein3)))
            dllik4=dllik4+y4(i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif

c......Partially observed 2 and 3 are collapsed
         else if (collap(i) .eq. 1) then
          if (y1(i) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
            dllik1=dllik1+
     1           dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           else
            llik1=llik1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
            dllik1=dllik1+y1(i)*
     1           (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
c
          sy23=y2(i)+y3(i)
          if (sy23 .gt. 0.0d0) then
           llik3a=dlog(1.0d0-dexp(ein3-ein1))
           if (ein1 .le. 0.0d0) then
            llik3b=ein1-dlog(1.0d0+dexp(ein1))
            dllik3b=-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1))
           else
            llik3b=-dlog(1.0d0+dexp(-ein1))
            dllik3b=-dexp(-ein1)/(1.0d0+dexp(-ein1))
           endif
           if (ein3 .le. 0.0d0) then
            llik3c=-dlog(1.0d0+dexp(ein3))
            dllik3c=dexp(ein3)/(1.0d0+dexp(ein3))
           else
            llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
            dllik3c=1.0d0-dexp(-ein3)/(1.0d0+dexp(-ein3))
           endif
           llik3=llik3+(llik3a+llik3b+llik3c)*sy23
           dllik3=dllik3+(dllik3b+dllik3c)*sy23
          endif
c
          if (y4(i) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4=llik4+y4(i)*
     1        (ein3-dlog(1.0d0+dexp(ein3)))
            dllik4=dllik4+y4(i)*
     1        (-1.0d0+dexp(ein3)/(1.0d0+dexp(ein3)))
           else
            llik4=llik4+y4(i)*
     1        (-dlog(1.0d0+dexp(-ein3)))
            dllik4=dllik4+y4(i)*
     1        (-dexp(-ein3)/(1.0d0+dexp(-ein3)))
           endif
          endif

c......Partially observed 2, 3, and 4 are collapsed
         else if (collap(i) .eq. 2) then
          if (y1(i) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
            dllik1=dllik1+
     1           dexp(ein1)/(1.0d0+dexp(ein1))*y1(i)
           else
            llik1=llik1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
            dllik1=dllik1+y1(i)*
     1           (1.0d0-dexp(-ein1)/(1.0d0+dexp(-ein1)))
           endif
          endif
c
          sy234=y2(i)+y3(i)+y4(i)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik4=llik4+
     1          (ein1-dlog(1.0d0+dexp(ein1)))*sy234
            dllik4=dllik4+
     1          (-1.0d0+dexp(ein1)/(1.0d0+dexp(ein1)))*sy234
           else
            llik4=llik4-dlog(1.0d0+dexp(-ein1))*sy234
            dllik4=dllik4-
     1          (dexp(-ein1)/(1.0d0+dexp(-ein1)))*sy234
           endif
          endif
         endif
        endif
390    continue
c
       sumLL=llik1+llik2+llik3+llik4
       sumDL=dllik1+dllik2+dllik3+dllik4
       ahz=sumLL-1.0d0/10000.0d0*z2**2
       ahpz=sumDL-1.0d0/5000.0d0*z2
       return
       end

c......Generation of alpha2
       subroutine Galpha2( iseed )
c      Yeongjin Gwon, UNMC
c      Jan 19, 2019
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       integer nStudy(n),trt(n),collap(n)
       integer iseed
c	   
       real*8 start,step,ximax,plmax,step1
       real*8 e1,fl(3,1),cl(5,3),dl(5,1),al(3,3),el(3,1)
       real*8 sigmaa,xistar,r1(1),r2(1),ratio
       integer konvge,kcount,ifault1,numres,icount
       external DRNNOF,DRNUNF,a2Nloglik
c	        
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu   
       common /dummy/idum3

c......Optimal variance searching
       do 600 kk=1,nx+1
        idum3=kk
        ximax=0.0d0
        plmax=0.0d0
        ratio=0.0d0
        r1(1)=0.0d0
        r2(1)=0.0d0
        nopt=1
        reqmin=1.0d-10
        konvge=5
        kcount=1000
        step=0.20d0
        start=alpha2(kk)
        call nelmin(a2Nloglik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
c
c       write(*,*) 'i=',kk,'alpha2=',ximax
        step1=0.5d0
 655    do i=1,5
         e1 = real(i-3)
         cl(i,1)=(ximax+e1*step1)**2
         cl(i,2)=ximax+e1*step1
         cl(i,3)=1.0d0
         dl(i,1)=a2Nloglik(ximax+e1*step1)
        enddo
c
        do ni=1,5
         if (ni .ne. 3) then
          if (dl(ni,1) .le. plmax) then
           step1=step1*1.50d0
           GOTO 655
          endif
         endif
        enddo

c......Estimate coefficients using LSE
        al = matmul(transpose(cl),cl)
        el = matmul(transpose(cl),dl)
        call lin_sol_self(al,el,fl)
c
        sigmaa=dsqrt(1.0d0/(fl(1,1)*2.0d0))
        call rnset( iseed )
        r1(1)=DRNNOF()
        call rnget(iseed)
        xistar=sigmaa*r1(1)+ximax
        ratio = -a2Nloglik(xistar)+a2Nloglik(start)
     1          -0.50d0*(start-ximax)**2/sigmaa**2
     2          +0.50d0*(xistar-ximax)**2/sigmaa**2

c......Perform Metropolis-Hastings algorithm to sample alpha2
        if (ratio .ge. 0.0d0) then
         alpha2(kk)=xistar
        else
         call rnset(iseed)
         r2(1)=DRNUNF()
         call rnget(iseed)
         if (dlog(r2(1)) .le. ratio) then
          alpha2(kk)=xistar
         endif
        endif
 600   continue
       return
       end

c......Negative loglikelihood function of alpha2
       real*8 function a2Nloglik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n),ylatent(n)
       real*8 sumLL,sumC
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 ein1,ein2,ein3  
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 xi,c1,c2,c3,eta2,eta3
       real*8 sy12,sy23,sy234      
       integer nStudy(n),trt(n),collap(n)
       integer iseed
c   
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu
       common /dummy/idum3
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       external DLINDS
       alpha2(idum3)=xi
c
       sumC=0.0d0
       sumLL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik1a=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
c       
       do 690 i=1,n
        if (trt(i) .gt. 0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c  
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Logliklihood construction
c......Fully observed categories
        if (collap(i) .eq. 0) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
          else
           llik1=llik1+
     1        (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c         llik1=llik1+llik1a*y1(i)
         endif
c          
         if (y2(i) .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein2 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein2))
          else
           llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
         endif
c          
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1        (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1        (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 1 and 2 are collapsed
        else if (collap(i) .eq. -1) then
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein2)))*sy12
          else
           llik1=llik1+
     1       (-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
          endif
c         llik1=llik1+llik1a*sy12
         endif
c          
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
            llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(i) .eq. 1) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein1)))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c         llik1=llik1+llik1a*y1(i)
         endif
c          
         sy23=y2(i)+y3(i)
         if (sy23 .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein3-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein3 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein3))
          else
           llik2c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*sy23
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2,3,and 4 are collapsed
        else if (collap(i) .eq. 2) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1+
     1      (-dlog(1.0d0+dexp(ein1)))*y1(i)
          else
           llik1=llik1+
     1      (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c         llik1=llik1+llik1a*y1(i)
         endif
c          
         sy234=y2(i)+y3(i)+y4(i)
         if (sy234 .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik4=llik4+
     1       (ein1-dlog(1.0d0+dexp(ein1)))*sy234
          else
           llik4=llik4+
     1       (dlog(1.0d0+dexp(-ein1)))*sy234
          endif
c         llik4=llik4+llik4a*sy234
         endif
        endif 
 690   continue
c
       sumLL=llik1+llik2+llik3+llik4
c       prior=-1.0d0/6.0d0*xi**2
       prior=-1.0d0/10000.0d0*xi**2
       a2Nloglik=-sumLL-prior
       return
       end

c......Generation of alpha3
       subroutine Galpha3( iseed )
c      Yeongjin Gwon, UNMC
c      Jan 19, 2019
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       integer nStudy(n),trt(n),collap(n)
       integer iseed,idum4
c	   
       real*8 start,step,ximax,plmax,step1
       real*8 e1,fl(3,1),cl(5,3),dl(5,1),al(3,3),el(3,1)
       real*8 sigmaa,xistar,r1(1),r2(1),ratio
       integer konvge,kcount,ifault1,numres,icount
       external DRNNOF,DRNUNF,a3Nloglik
c	        
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu   
       common /dummy/idum4

c......Optimal variance searching
       do 700 kk=1,nx+1
        idum4=kk
        ximax=0.0d0
        plmax=0.0d0
        ratio=0.0d0
        r1(1)=0.0d0
        r2(1)=0.0d0
        nopt=1
        reqmin=1.0d-10
        konvge=5
        kcount=1000
        step=0.20d0
        start=alpha3(kk)
        call nelmin(a3Nloglik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
c
        step1=0.5d0
 755    do i=1,5
         e1 = real(i-3)
         cl(i,1)=(ximax+e1*step1)**2
         cl(i,2)=ximax+e1*step1
         cl(i,3)=1.0d0
         dl(i,1)=a3Nloglik(ximax+e1*step1)
        enddo
c
        do ni=1,5
         if (ni .ne. 3) then
          if(dl(ni,1) .le. plmax) then
           step1=step1*1.5d0
           GOTO 755
          endif
         endif
        enddo

c......Estimate coefficients using LSE
        al = matmul(transpose(cl),cl)
        el = matmul(transpose(cl),dl)
        call lin_sol_self(al,el,fl)
c
        sigmaa=dsqrt(1.0d0/(fl(1,1)*2.0d0))
        call rnset( iseed )
        r1(1)=DRNNOF()
        call rnget( iseed )
        xistar=sigmaa*r1(1)+ximax
        ratio = -a3Nloglik(xistar)+a3Nloglik(start)
     1          -0.50d0*(start-ximax)**2/sigmaa**2
     2          +0.50d0*(xistar-ximax)**2/sigmaa**2

c......Perform Metropolis-Hastings algorithm to sample alpha3
        if (ratio .ge. 0.0d0) then
         alpha3(kk)=xistar
        else
         call rnset(iseed)
         r2(1)=DRNUNF()
         call rnget(iseed)
         if (dlog(r2(1)) .le. ratio) then
          alpha3(kk)=xistar
         endif
        endif
 700   continue
       return
       end

c......Negative loglikelihood function of alpha3
       real*8 function a3Nloglik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n),ylatent(n)
       real*8 sumLL,sumC
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 ein1,ein2,ein3
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 xi,c1,c2,c3,eta2,eta3
       real*8 sy12,sy23,sy234
       integer nStudy(n),trt(n),collap(n)
       integer idum4
c   
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu   
       common /dummy/idum4
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       alpha3(idum4)=xi
c
       sumC=0.0d0
       sumLL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik1a=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
c	 
       do 790 i=1,n
        if (trt(i) .gt. 0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c  
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Logliklihood construction
c......Fully observed categories
        if (collap(i) .eq. 0) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein1)))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c         llik1=llik1+llik1a*y1(i)
         endif
c          
         if (y2(i) .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein2 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein2))
          else
           llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
         endif
c          
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 1 and 2 are collapsed
        else if (collap(i) .eq. -1) then
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein2)))*sy12
          else
           llik1=llik1+
     1       (-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
          endif
c         llik1=llik1+llik1a*sy12
         endif
c          
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(i) .eq. 1) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein1)))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c          llik1=llik1+llik1a*y1(i)
         endif
c          
         sy23=y2(i)+y3(i)
         if (sy23 .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein3-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein3 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein3))
          else
           llik2c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*sy23
         endif
c          
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           llik4=llik4+
     1       (ein3-dlog(1.0d0+dexp(ein3)))*y4(i)
          else
           llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein3)))*y4(i)
          endif
c         llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2,3,and 4 are collapsed
        else if (collap(i) .eq. 2) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1+
     1       (-dlog(1.0d0+dexp(ein1)))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
c         llik1=llik1+llik1a*y1(i)
         endif
c          
         sy234=y2(i)+y3(i)+y4(i)
         if (sy234 .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
            llik4=llik4+
     1       (ein1-dlog(1.0d0+dexp(ein1)))*sy234
          else
            llik4=llik4+
     1       (-dlog(1.0d0+dexp(-ein1)))*sy234
          endif
c         llik4=llik4+llik4a*sy234
         endif
        endif 
 790   continue
c 
       sumLL=llik1+llik2+llik3+llik4
c       prior=-1.0d0/6.0d0*xi**2
       prior=-1.0d0/10000.0d0*xi**2
       a3Nloglik=-sumLL-prior
       return
       end


c......Generating auxiliary variables bounds
       subroutine Gb( iseed )
c      Yeongjin, UNMC
c      June 14, 2019
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       integer nStudy(n),trt(n),collap(n)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1),b(n,4)
       real*8 ss(n),n70(n),n100(n),nrem(n),ylatent(n)
       real*8 aBL(n,4),aBU(n,4),maxBL(n),minBU(n)
       real*8 sumC
       real*8 aLmuk(K),aUmuk(K)
       real*8 lein,ein1,ein2,ein3,ein4   
       real*8 c1,c2,c3
       real*8 eb1,eb2,eb3,eb4
       real*8 eb2a,eb2b,eb2c,eb3a,eb3b,eb3c,eb4a
       real*8 sy12,sy23,sy234
       real*8 tRTnt,xbeta,zbeta,eta2,eta3   
       external Glow,DRNUNF
c
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecb/b
       common /vecmu/mu
       common /vecaBL/aBL
       common /vecaBU/aBU
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk
       common /vecylatent/ylatent
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
c
       sumC=0.0d0
       eb1=0.0d0
       eb2=0.0d0
       eb3=0.0d0
       eb4=0.0d0
       eb2a=0.0d0
       eb2b=0.0d0
       eb2c=0.0d0
       eb3a=0.0d0
       eb3b=0.0d0
       eb3c=0.0d0 
       eb4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0   
       lein=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       ein4=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
c 
       do 300 i=1,n
        if (trt(i) .gt. 0.0d0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif 
c	
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        if (i .eq. 35) then
         xbeta=x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        endif
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
        lein=tRTnt+xbeta
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Lower and upper bound calculation using developed Theorem 3-2
c......Fully observed categories	
        if (collap(i) .eq. 0) then
         call rnset( iseed )
          rn1=DRNUNF()
         call rnget( iseed )
         call rnset( iseed )
          rn2=DRNUNF()
         call rnget( iseed )
         call rnset( iseed )
          rn3=DRNUNF()
         call rnget( iseed )
         call rnset( iseed )
          rn4=DRNUNF()
         call rnget( iseed )
c
         eb1=dlog(rn1)
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           eb1=eb1-dlog(1.0d0+dexp(ein1))*y1(i)   
          else
           eb1=eb1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
         b(i,1)=eb1
c
         eb2=dlog(rn2)
         if (y2(i) .gt. 0.0d0) then
          eb2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           eb2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           eb2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein2 .le. 0.0d0) then
           eb2c=-dlog(1.0d0+dexp(ein2))
          else
           eb2c=-ein2-dlog(1.0d0+dexp(-ein2))
          endif
          eb2=eb2+(eb2a+eb2b+eb2c)*y2(i)
         endif
         b(i,2)=eb2
c
         eb3=dlog(rn3)
         if (y3(i) .gt. 0.0d0) then
          eb3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           eb3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           eb3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           eb3c=-dlog(1.0d0+dexp(ein3))
          else
           eb3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          eb3=eb3+(eb3a+eb3b+eb3c)*y3(i)
         endif
         b(i,3)=eb3
c
         eb4=dlog(rn4)
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           eb4=eb4+y4(i)*(ein3-dlog(1.0d0+dexp(ein3)))
          else
           eb4=eb4+y4(i)*(-dlog(1.0d0+dexp(-ein3)))
          endif
         endif
         b(i,4)=eb4

c......The 1st category only needs to have lower bound
         aBL(i,1)=-100000.0d0
         aBU(i,1)=100000.0d0
         if (y1(i) .gt. 0.0d0) then
          eb1=b(i,1)/y1(i)
          aBL(i,1)=eb1-dlog(1.0d0-dexp(eb1))-lein-c1
         endif

c.......The 2nd category needs two bounds		 		 
         aBL(i,2)=-100000.0d0
         aBU(i,2)=100000.0d0
         if (y2(i) .gt. 0.0d0) then
          eb2=b(i,2)/y2(i)
          ein1=Glow(c2-c1,eb2)
          aBL(i,2)=-(c2+c1)/2.0d0-dexp(ein1)-lein
          aBU(i,2)=-(c2+c1)/2.0d0+dexp(ein1)-lein
         endif

c.......The 3rd category needs two bounds		 
         aBL(i,3)=-100000.0d0
         aBU(i,3)=100000.0d0
         if (y3(i) .gt. 0.0d0) then
          eb3=b(i,3)/y3(i)
          ein2=Glow(c3-c2,eb3)
          aBL(i,3)=-(c3+c2)/2.0d0-dexp(ein2)-lein
          aBU(i,3)=-(c3+c2)/2.0d0+dexp(ein2)-lein
         endif

c......The 4th category only needs to have upper bound
         aBL(i,4)=-100000.0d0
         aBU(i,4)=100000.0d0
         if (y4(i) .gt. 0.0d0) then
          eb4=b(i,4)/y4(i)
          aBU(i,4)=-eb4+dlog(1.0d0-dexp(eb4))-lein-c3
         endif

c......Partially observed 1 and 2 collapsing		  
        else if (collap(i) .eq. -1) then
         call rnset( iseed )
          rn1=DRNUNF()
         call rnget( iseed )
         call rnset( iseed )
          rn3=DRNUNF()
         call rnget( iseed )
         call rnset( iseed )
          rn4=DRNUNF()
         call rnget( iseed )
c		 
         eb1=dlog(rn1)
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           eb1=eb1-dlog(1.0d0+dexp(ein2))*sy12
          else
           eb1=eb1+(-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
          endif
         endif
         b(i,1)=eb1
c		 
         eb3=dlog(rn3)
         if (y3(i) .gt. 0.0d0) then
          eb3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           eb3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           eb3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           eb3c=-dlog(1.0d0+dexp(ein3))
          else
           eb3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          eb3=eb3+(eb3a+eb3b+eb3c)*y3(i)
         endif
         b(i,3)=eb3
c
         eb4=dlog(rn4)
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .le. 0.0d0) then
           eb4=eb4+y4(i)*(ein3-dlog(1.0d0+dexp(ein3)))
          else
           eb4=eb4+y4(i)*(-dlog(1.0d0+dexp(-ein3)))
          endif
         endif
         b(i,4)=eb4

c......The collapsing 1st and 2nd needs only lower bound
         aBL(i,1)=-100000.0d0
         aBU(i,1)=100000.0d0
         if (sy12 .gt. 0.0d0) then
          eb1=b(i,1)/sy12
          aBL(i,1)=eb1-dlog(1.0d0-dexp(eb1))-lein-c2
         endif
         aBL(i,2)=aBL(i,1)
         aBU(i,2)=aBU(i,1)

c......The 3rd category needs to have lower and upper bounds
         aBL(i,3)=-100000.0d0
         aBU(i,3)=100000.0d0
         if (y3(i) .gt. 0.0d0) then
          eb3=b(i,3)/y3(i)
          ein2=Glow(c3-c2,eb3)
          aBL(i,3)=-(c3+c2)/2.0d0-dexp(ein2)-lein
          aBU(i,3)=-(c3+c2)/2.0d0+dexp(ein2)-lein
         endif

c......The 4th category needs upper bound
         aBL(i,4)=-100000.0d0
         aBU(i,4)=100000.0d0
         if (y4(i) .gt. 0.0d0) then
          eb4=b(i,4)/y4(i)
          aBU(i,4)=-eb4+dlog(1.0d0-dexp(eb4))-lein-c3
         endif

c......Partially observed 2 and 3 collapsing		   
        else if (collap(i) .eq. 1) then
          call rnset( iseed )
           rn1=DRNUNF()
          call rnget( iseed )
          call rnset( iseed )
           rn2=DRNUNF()
          call rnget( iseed )
          call rnset( iseed )
           rn4=DRNUNF()
          call rnget( iseed )
c		 
          eb1=dlog(rn1)
          if (y1(i) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            eb1=eb1-dlog(1.0d0+dexp(ein1))*y1(i)
           else
            eb1=eb1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
           endif
          endif
          b(i,1)=eb1
c		 
          eb2=dlog(rn2)
          sy23=y2(i)+y3(i)
          if (sy23 .gt. 0.0d0) then
           eb2a=dlog(1.0d0-dexp(ein3-ein1))
           if (ein1 .le. 0.0d0) then
            eb2b=ein1-dlog(1.0d0+dexp(ein1))
           else
            eb2b=-dlog(1.0d0+dexp(-ein1))
           endif
           if (ein3 .le. 0.0d0) then
            eb2c=-dlog(1.0d0+dexp(ein3))
           else
            eb2c=-ein3-dlog(1.0d0+dexp(-ein3))
           endif
           eb2=eb2+(eb2a+eb2b+eb2c)*sy23
          endif
          b(i,2)=eb2
c		 
          eb4=dlog(rn4)
          if (y4(i) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            eb4=eb4+y4(i)*(ein3-dlog(1.0d0+dexp(ein3)))
           else
            eb4=eb4+y4(i)*(-dlog(1.0d0+dexp(-ein3)))
           endif
          endif
          b(i,4)=eb4

c......The 1st category needs only lower bound
          aBL(i,1)=-100000.0d0
          aBU(i,1)=100000.0d0
          if (y1(i) .gt. 0.0d0) then
           eb1=b(i,1)/y1(i)
           aBL(i,1)=eb1-dlog(1.0d0-dexp(eb1))-lein-c1
          endif

c......The collapsing 2nd and 3rd categories need to have lower and upper bound
          aBL(i,2)=-100000.0d0
          aBU(i,2)=100000.0d0
          if (sy23 .gt. 0.0d0) then
           eb2=b(i,2)/sy23
           ein1=Glow(c3-c1,eb2)
           aBL(i,2)=-(c3+c1)/2.0d0-dexp(ein1)-lein
           aBU(i,2)=-(c3+c1)/2.0d0+dexp(ein1)-lein
          endif
          aBL(i,3)=aBL(i,2)
          aBU(i,3)=aBU(i,2)

c......The 4th category needs only upper bound
          aBL(i,4)=-100000.0d0
          aBU(i,4)=100000.0d0
          if (y4(i) .gt. 0.0d0) then
           eb4=b(i,4)/y4(i)
           aBU(i,4)=-eb4+dlog(1.0d0-dexp(eb4))-lein-c3
          endif

c......Partially observed 2,3,and 4 are collapsed	   
        else if (collap(i) .eq. 2) then
          call rnset( iseed )
           rn1=DRNUNF()
          call rnget( iseed )
          call rnset( iseed )
           rn4=DRNUNF()
          call rnget( iseed )
c		 
          eb1=dlog(rn1)
          if (y1(i) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            eb1=eb1-dlog(1.0d0+dexp(ein1))*y1(i)
           else
            eb1=eb1+(-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
           endif
          endif
          b(i,1)=eb1
c		 
          eb4=dlog(rn4)
          sy234=y2(i)+y3(i)+y4(i)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            eb4=eb4+(ein1-dlog(1.0d0+dexp(ein1)))*sy234
           else
            eb4=eb4-dlog(1.0d0+dexp(-ein1))*sy234
           endif
          endif
          b(i,4)=eb4
c	 
c......The 1st category needs only lower bound
          aBL(i,1)=-100000.0d0
          aBU(i,1)=100000.0d0
          if (y1(i) .gt. 0.0d0) then
           eb1=b(i,1)/y1(i)
           aBL(i,1)=eb1-dlog(1.0d0-dexp(eb1))-lein-c1
          endif

c......The collapsing 2,3,and 4th categories need to have lower and upper bound
          aBL(i,4)=-100000.0d0
          aBU(i,4)=100000.0d0
          if (sy234 .gt. 0.0d0) then
           eb4=b(i,4)/sy234
           aBU(i,4)=-eb4+dlog(1.0d0-dexp(eb4))-lein-c1
          endif
          aBL(i,3)=aBL(i,4)
          aBU(i,3)=aBU(i,4)
          aBL(i,2)=aBL(i,4)
          aBU(i,2)=aBU(i,4)
        endif

c......Find maximum of lower and minimum of upper bounds
        eBL=aBL(i,1)
        eBU=aBU(i,1)
        do j=2,4
         if (eBL .lt. aBL(i,j)) eBL=aBL(i,j)
         if (eBU .gt. aBU(i,j)) eBU=aBU(i,j)
        enddo
        maxBL(i)=eBL
        minBU(i)=eBU
 300   continue   
c......Checking direction for maximum lower and minimum upper bounds
c       do i=1,n
c        write(*,1000) i,maxBL(i),minBU(i)
c 1000   format(1x,I2,1x,2f10.8)
c       enddo   
       return
       end


c......Optimization to find lower and upper bounds	   
       real*8 function Glow(a1,a2)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 mode1,temp,a1,a2
       real*8 start,step,ximax,plmax,reqmin
       integer konvge,kcount,ifault1,numres,icount,nopt
       external alowlik
       common /vecdcut/dcut
       common /vecdbdd/dbdd

c......Set up initial variables
       dcut=a1
       dbdd=a2
       plmax=0.0d0
       nopt=1
       reqmin=1.0d-10
       konvge=5
       kcount=1000
       step=0.2d0
       start=dcut
       call nelmin(alowlik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
       Glow=ximax
c      write(*,*) 'Optimal Low=',ximax
       return
       end
  
c......Objective function for lowrer bound calculation  
       real*8 function alowlik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 temp1,temp2,sumL,sumT,ein,ein2
       integer iseed
       common /vecdcut/dcut
       common /vecdbdd/dbdd
c
       temp1=0.0d0
       temp2=0.0d0
       sumL=0.0d0
       sumT=0.0d0
       ein1=0.0d0
       ein2=0.0d0
c       
       temp1=dcut/2.0d0-dexp(xi)
       temp2=-dcut/2.0d0-dexp(xi)
       temp1=-temp1
       temp2=-temp2
c......Previous code	   
c       sumL=temp1+dlog(1.0d0-dexp(temp2-temp1))
c     1    -dlog(1.0d0+dexp(temp2))
c     1    -dbdd	 
c       if (temp1 .le. 0.0d0) then
c         ein=-dlog(1.0d0+dexp(temp1))
c       else
c         ein=-temp1-dlog(1.0d0+dexp(-temp1))
c       endif
c       ein2=sumL+ein
c......Current code
       sumL=dlog(1.0d0-dexp(temp1-temp2))
     1       -dbdd 
       if (temp1 .le. 0.0d0) then
        ein1=-dlog(1.0d0+dexp(temp1))
       else
        ein1=-temp1-dlog(1.0d0+dexp(-temp1))
       endif
       if (temp2 .le. 0.0d0) then
        ein2=temp2-dlog(1.0d0+dexp(temp2))
       else
        ein2=-dlog(1.0d0+dexp(-temp2))
       endif
       sumT=sumL+ein1+ein2
c......Revised code
c       sumL=temp1+dlog(1.0d0-dexp(temp2-temp1))
c     1    -dlog(1.0d0+dexp(temp1))
c     2    -dlog(1.0d0+dexp(temp2))
c     3    -dbdd
c       sumT=sumL
       if (sumT .lt. 0.0d0) sumT=-sumT
       alowlik=sumT
       return
       end


c......Generation of phi
       subroutine Gphi( iseed )
c      Yeongjin Gwon, UNMC
c      Jan 30, 2019
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 phi(nz+1)
       real*8 x(nx,n),z(nz,n)
c	   
       real*8 start,step,ximax,plmax,step1,reqmin
       real*8 e1,fl(3,1),cl(5,3),dl(5,1),al(3,3),el(3,1)
       real*8 sigmaa,xistar,r1(1),r2(1),ratio
       real*8 aBL(n,4),aBU(n,4),maxBL(n),minBU(n)
       real*8 aLmuk(K),aUmuk(K)
       real*8 emax,emin,low,upp
       integer iseed
       integer konvge,kcount,ifault1,numres,icount,nopt
c
       common /vecx/x
       common /vecz/z
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk
       common /vecphi/phi
       common /dummy/idum0
       external DRNNOF,DRNUNF,pNloglik
c
       emax=0.0d0
       emin=0.0d0

c......Optimal variance searching
       do 100 kk=1,2
        idum0=kk
        ximax=0.0d0
        plmax=0.0d0
        ratio=0.0d0
        r1(1)=0.0d0
        r2(1)=0.0d0
        nopt=1
        reqmin=1.0d-10
        konvge=5
        kcount=1000
        step=0.2d0
        start=phi(kk)
        call nelmin(pNloglik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
c
c       write(*,*) 'value =',plmax
        step1=0.5d0
 155    do i=1,5
         e1 = real(i-3)
         cl(i,1)=(ximax+e1*step1)**2
         cl(i,2)=ximax+e1*step1
         cl(i,3)=1.0d0
         dl(i,1)=pNloglik(ximax+e1*step1)
        enddo
c
        do i=1,5
         if (i .ne. 3) then
          if(dl(i,1) .le. plmax) then
           step1=step1*1.2d0
           GOTO 155
          endif
         endif
        enddo

c......Estimate coefficients using LSE
        al = matmul(transpose(cl),cl)
        el = matmul(transpose(cl),dl)
        call lin_sol_self(al,el,fl)
        sigmaa=dsqrt(1.0d0/(fl(1,1)*2.0d0))
c       write(*,*) 'Value of a =',fl(1,1)
c       write(*,*) 'New sigma =',sigmaa
        call rnset( iseed )
        r1(1)=DRNNOF()
        call rnget( iseed )
        xistar=sigmaa*r1(1)+ximax
        ratio = -pNloglik(xistar)+pNloglik(start)
     1          -0.50d0*(start-ximax)**2/sigmaa**2
     2          +0.50d0*(xistar-ximax)**2/sigmaa**2

c......Perform Metropolis-Hastings algorithm to sample from phi
        if (ratio .ge. 0.0d0) then
         phi(kk)=xistar
        else
         call rnset( iseed )
         r2(1)=DRNUNF()
         call rnget( iseed )
         if (dlog(r2(1)) .le. ratio) then
          phi(kk)=xistar
         endif
        endif
 100   continue
       return
       end

c......Negative loglikelihood function of phi
       real*8 function pNloglik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 x(nx,n),z(nz,n)
       real*8 phi(nz+1)
       real*8 sum1,sum2,sum3,sum4
       real*8 sumLL,sumC
       real*8 adif,aprod,emax,emin
       real*8 lik1,lik2,lik3,lik4
       real*8 alow,alow1,alik2
       real*8 ein,ein1,ein2,ein3,ein4,ein5
       real*8 eein2,eein3
       real*8 xi,c1,c2,c3,eta2,eta3
       real*8 low,upp,eps,temp1,temp2,low1,upp1   
       real*8 aBL(n,4),aBU(n,4),maxBL(n),minBU(n)
       real*8 aLmuk(K),aUmuk(K)
       integer iseed
c   
       common /vecx/x
       common /vecz/z
       common /vecphi/phi
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecaBL/aBL
       common /vecaBU/aBU  
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk
       common /dummy/idum0
       common /vecisim/isim
       external DRNNOF,DRNUNF,DNORDF
c       
       phi(idum0)=xi
       sumC=0.0d0
       sumLL=0.0d0  
       sum1=0.0d0
       sum2=0.0d0
       sum3=0.0d0
       sum4=0.0d0
       emax=0.0d0
       emin=0.0d0
       ein=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       ein4=0.0d0
       ein5=0.0d0
       eein2=0.0d0
       eein3=0.0d0
       temp=0.0d0

c......Calculate normal cdf of random effects	
       do 1000 k1=1,K
        emax=maxBL(2*k1-1)
        if (emax .lt. maxBL(2*k1)) emax=maxBL(2*k1)
        emin=minBU(2*k1-1)
        if (emin .gt. minBU(2*k1)) emin=minBU(2*k1)
        temp=phi(1)+phi(2)*z(1,2*k1)
        low=emax/dexp(temp)
        upp=emin/dexp(temp)
c......Stopping criterion
        if (upp .le. low) then
         write(*,*) 'in pNloglik ERROR ','low=',low,'upp=',upp
         temp=upp
         upp=low
         low=temp
        endif
        adif=upp-low
        aprod=upp*low

c......Calculation loglikelihood 
c......Opposite sign (affect the speed)
        if (aprod .lt. 0.0d0) then
c          ein2=DNORDF(upp)-DNORDF(low)	
c......Above coding is too slow in generating phi
          if (upp .ge. 5.3d0) then
           eein2=1.0d0
          else 
           eein2=DNORDF(upp)
          endif
          if (low .le. -5.3d0) then 
           eein3=0.0d0
          else
           eein3=DNORDF(low)
          endif
          ein1=eein2-eein3
          sum1=sum1+dlog(ein1)
c         write(*,*) 'Likelihood 1=',sum1,'LB=',low,'UB=',upp		 

c......Same sign		 
        else
        if ((adif .gt. 0.0d0) .and. 
     1       (adif .le. 0.01d0)) then
          ein2=-((upp+low)/2.0d0)**2/2.0d0+dlog(adif)
          sum2=sum2+ein2
c         write(*,*) 'Likelihood 2=',sum2	  
        else
          if ((upp .gt. 0.0d0) .and. 
     1         (upp .lt. 5.3d0)) then
            ein3=DNORDF(upp)-DNORDF(low)
c           write(*,*) '1: upp=',upp,' low=',low
c           write(*,*) 'DNORDF: ein2=',ein2
            sum3=sum3+dlog(ein3)
          else if (upp .ge. 5.3d0) then
            if (low .lt. 5.2d0) then
              ein3=1.0d0-DNORDF(low)
c             ein3=DNORDF(upp)-DNORDF(low)
              sum3=sum3+dlog(ein3)
c             write(*,*) '2: upp=',upp,' low=',low
c             write(*,*) 'DNORDF: ein2=',ein2
            else
              ein=low
              ein3=-ein**2/2.0d0
     1          +dlog(1.0d0/ein-1.0d0/ein**3+3.0d0/ein**5)
              sum3=sum3+ein3
            endif
c           write(*,*) 'Likelihood 3=',sum3	   
          else if (upp .lt. 0.0d0) then
            upp1=-low
            low1=-upp
            if (upp1 .lt. 5.3d0) then
              ein4=DNORDF(upp1)-DNORDF(low1)
c             write(*,*) '3: upp1=',upp1,' low1=',low1
c             write(*,*) 'DNORDF: ein2=',ein2
              sum4=sum4+dlog(ein4)
            else if (upp1 .gt. 5.3d0) then
              if (low1 .lt. 5.2d0) then
                ein4=1.0d0-DNORDF(low1)
c               ein4=DNORDF(upp1)-DNORDF(low1)
c               write(*,*) '4: upp1=',upp1,' low1=',low1
c               write(*,*) 'DNORDF: ein2=',ein2
                sum4=sum4+dlog(ein4)
              else
                ein=low1
                ein4=-ein**2/2.0d0
     1             +dlog(1.0d0/ein-1.0d0/ein**3+3.0d0/ein**5)
                sum4=sum4+ein4
              endif
            endif
c           write(*,*) 'Likelihood 4=',sum4		   
          endif
        endif
       endif
 1000  continue
c
       sumLL=sum1+sum2+sum3+sum4
c      ein5=-1.0d0/9.0d0*xi**2
       ein5=-1.0d0/10000.0d0*xi**2
       pNloglik=-sumLL-ein5
       return
       end


c......Generation of trial level random effects	   
       subroutine GmuK( iseed )
c      Yeongjin Gwon, UNMC
c      Jan 14, 2019
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 aLmuk(K),aUmuk(K),mu(K) 
       real*8 phi(nz+1),emax,emin,z(nz,n)
       real*8 maxBL(n),minBU(n)
       real*8 aleft,aright
       integer iseed     
       logical la,lb 
c
       common /vecz/z
       common /vecphi/phi
       common /vecmaxBL/maxBL
       common /vecminBU/minBU
       common /vecmu/mu
       common /vecaLmuk/aLmuk
       common /vecaUmuk/aUmuk
       external YTUVN
c	   
       temp=0.0d0
       la=.false.
       lb=.false.
       call rnset( iseed )
       do 3000 j=1,K
        emax=maxBL(2*j-1)
        if (emax .lt. maxBL(2*j)) emax=maxBL(2*j)
        emin=minBU(2*j-1)
        if (emin .gt. minBU(2*j)) emin=minBU(2*j)
        temp=phi(1)+phi(2)*z(1,2*j)
        aLmuk(j)=emax/dexp(temp)
        aUmuk(j)=emin/dexp(temp)
        aleft=aLmuk(j)
        aright=aUmuk(j)
        mu(j)=YTUVN(aleft,aright,la,lb,iseed)
 3000  continue
       call rnget( iseed )
       return
       end   


c......Optimization to get maximum points	   
       subroutine Optdelta
c      Yeongjin Gwon, UNMC
c      Jan 19, 2019
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 delta(nT)
       real*8 start,step,ximax,plmax,reqmin
       integer konvge,kcount,ifault1,numres,icount,nopt
       integer iseed
       external DRNNOF,DRNUNF,dNloglik
c	        
       common /vecdelta/delta
       common /dummy/idum2

c......Optimal variance searching
       ximax=0.0d0
       plmax=0.0d0
       nopt=1
       reqmin=1.0d-10
       konvge=5
       kcount=1000
       step=0.2d0
       start=delta(idum2)
       call nelmin(dNloglik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
c
c      write(*,*) 'Optimization =',ximax
       return
       end

c......Negative loglikelihood function of deltas
       real*8 function dNloglik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1)
       real*8 sumLL
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 ein1,ein2,ein3
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 xi,c1,c2,c3,eta2,eta3
       integer nStudy(n),trt(n),collap(n)
       integer iseed
c
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3   
       common /vecmu/mu  
       common /dummy/idum2
       delta(idum2)=xi
c  
       sumLL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
c
       do 1000 i=1,n
        if (trt(i) .gt. 0.0d0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c		
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Logliklihood construction
c......Fully observed categories
        if (collap(i) .eq. 0) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)   
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
c         
         if (y2(i) .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein2 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein2))
          else
           llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
           llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 1 and 2 are collapsed
        else if (collap(i) .eq. -1) then
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein2))*sy12
          else
           llik1=llik1+
     1       (-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
          endif
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
           llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(i) .eq. 1) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
          else
           llik1=llik1+
     1      (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
         sy23=y2(i)+y3(i)
c         
         if (sy23 .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein3-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein3 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein3))
          else
           llik2c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*sy23
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
           llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2,3,and 4 are collapsed
        else if (collap(i) .eq. 2) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
          else
           llik1=llik1+
     1        (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
c         
         sy234=y2(i)+y3(i)+y4(i)
         if (sy234 .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik4=llik4+
     1        (ein1-dlog(1.0d0+dexp(ein1)))*sy234
          else
           llik4=llik4-dlog(1.0d0+dexp(-ein1))*sy234
          endif
         endif
        endif
 1000  continue
c 
       sumLL=llik1+llik2+llik3+llik4
       prior=-1.0d0/10000.0d0*xi**2
       dNloglik=-sumLL-prior
       return
       end


       subroutine Optbeta
c      Yeongjin Gwon, UNMC
c      Jan 19, 2019
       use linear_operators
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 beta(nx+1)
       real*8 start,step,ximax,plmax,reqmin
       integer konvge,kcount,ifault1,numres,icount,nopt
       integer iseed
       external DRNNOF,DRNUNF,bNloglik
c	   
       common /vecbeta/beta   
       common /dummy/idum1

c......Optimal variance searching
       ximax=0.0d0
       plmax=0.0d0
       nopt=1
       reqmin=1.0d-10
       konvge=5
       kcount=1000
       step=0.2d0
       start=beta(idum1)
       call nelmin(bNloglik,nopt,start,ximax,plmax,reqmin,
     1         step,konvge,kcount,icount,numres,ifault1)
c
c      write(*,*) 'Optimization =',ximax
       return
       end

c......Negative loglikelihood function of beta
       real*8 function bNloglik(xi)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 y1(n),y2(n),y3(n),y4(n),x(nx,n),z(nz,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1)
       real*8 sumLL,sumC
       real*8 llik1,llik2,llik3,llik4  
       real*8 llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 ein1,ein2,ein3
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 xi,c1,c2,c3,eta2,eta3
       integer nStudy(n),trt(n),collap(n)
       integer iseed
c   
       common /vecy1/y1
       common /vecy2/y2
       common /vecy3/y3
       common /vecy4/y4
       common /vecx/x
       common /vecz/z
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3
       common /vecmu/mu 
       common /dummy/idum1
       beta(idum1)=xi
c  
       sumC=0.0d0
       sumLL=0.0d0
       llik1=0.0d0
       llik2=0.0d0
       llik3=0.0d0
       llik4=0.0d0
       llik2a=0.0d0
       llik2b=0.0d0
       llik2c=0.0d0
       llik3a=0.0d0
       llik3b=0.0d0
       llik3c=0.0d0
       llik4a=0.0d0
       xbeta=0.0d0
       zbeta=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       muK=0.0d0
       tRTnt=0.0d0
c
       do 2000 i=1,n
        if (trt(i) .gt. 0.0d0) then
         tRTnt=delta(trt(i))
        else
         tRTnt=0.0d0
        endif
c		
        muK=mu(nStudy(i))
        zbeta=phi(1)+z(1,i)*phi(2)
        xbeta=beta(1)+x(1,i)*beta(2)
     1     +x(2,i)*beta(3)+x(3,i)*beta(4)
     2     +x(4,i)*beta(5)+x(5,i)*beta(6)
        eta2=alpha2(1)+x(1,i)*alpha2(2)
     1     +x(2,i)*alpha2(3)+x(3,i)*alpha2(4)
     2     +x(4,i)*alpha2(5)+x(5,i)*alpha2(6)
        eta3=alpha3(1)+x(1,i)*alpha3(2)
     1     +x(2,i)*alpha3(3)+x(3,i)*alpha3(4)
     2     +x(4,i)*alpha3(5)+x(5,i)*alpha3(6)
        c1=0.0d0
        c2=c1+dexp(eta2)
        c3=c2+dexp(eta3)
c		
        ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
        ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
        ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)

c......Logliklihood construction
c......Fully observed categories
        if (collap(i) .eq. 0) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)   
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
c         
         if (y2(i) .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein2-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein2 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein2))
          else
           llik2c=-ein2-dlog(1.0d0+dexp(-ein2))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*y2(i)
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
           llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 1 and 2 are collapsed
        else if (collap(i) .eq. -1) then
         sy12=y1(i)+y2(i)
         if (sy12 .gt. 0.0d0) then
          if (ein2 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein2))*sy12
          else
           llik1=llik1+
     1       (-ein2-dlog(1.0d0+dexp(-ein2)))*sy12
          endif
         endif
c         
         if (y3(i) .gt. 0.0d0) then
          llik3a=dlog(1.0d0-dexp(ein3-ein2))
          if (ein2 .le. 0.0d0) then
           llik3b=ein2-dlog(1.0d0+dexp(ein2))
          else
           llik3b=-dlog(1.0d0+dexp(-ein2))
          endif
          if (ein3 .le. 0.0d0) then
           llik3c=-dlog(1.0d0+dexp(ein3))
          else
           llik3c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik3=llik3+(llik3a+llik3b+llik3c)*y3(i)
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
           llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(i) .eq. 1) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
c         
         sy23=y2(i)+y3(i)
         if (sy23 .gt. 0.0d0) then
          llik2a=dlog(1.0d0-dexp(ein3-ein1))
          if (ein1 .le. 0.0d0) then
           llik2b=ein1-dlog(1.0d0+dexp(ein1))
          else
           llik2b=-dlog(1.0d0+dexp(-ein1))
          endif
          if (ein3 .le. 0.0d0) then
           llik2c=-dlog(1.0d0+dexp(ein3))
          else
            llik2c=-ein3-dlog(1.0d0+dexp(-ein3))
          endif
          llik2=llik2+(llik2a+llik2b+llik2c)*sy23
         endif
c         
         if (y4(i) .gt. 0.0d0) then
          if (ein3 .lt. 0.0d0) then
           llik4a=ein3-dlog(1.0d0+dexp(ein3))
          else
            llik4a=-dlog(1.0d0+dexp(-ein3))
          endif
          llik4=llik4+llik4a*y4(i)
         endif

c......Partially observed 2,3,and 4 are collapsed
        else if (collap(i) .eq. 2) then
         if (y1(i) .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik1=llik1-dlog(1.0d0+dexp(ein1))*y1(i)
          else
           llik1=llik1+
     1       (-ein1-dlog(1.0d0+dexp(-ein1)))*y1(i)
          endif
         endif
c         
         sy234=y2(i)+y3(i)+y4(i)
         if (sy234 .gt. 0.0d0) then
          if (ein1 .le. 0.0d0) then
           llik4=llik4+(ein1-dlog(1.0d0+dexp(ein1)))*sy234
          else
           llik4=llik4-dlog(1.0d0+dexp(-ein1))*sy234
          endif
         endif
        endif
 2000  continue
c
       sumLL=llik1+llik2+llik3+llik4
       prior=-1.0d0/10000.0d0*xi**2
       bNloglik=-sumLL-prior
       return
       end

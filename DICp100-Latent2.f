c......Calculation DIC using version DIC7 
c......DIC for missing data models Bayesian analysis (2006)
c
c......Calculate barDIC: average of posterior deviance function
       subroutine DIC1(xi,res1)
       parameter (n=36,K=18,nT=13,nx=5,nz=5,nmax=500)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1)
       real*8 theta(K),prior
       real*8 llik,llik1,llik2,llik3,llik4  
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 c1,c2,c3,eta2,eta3
       real*8 p1,p2,p3
       real*8 xi,ein1,ein2,ein3
       real*8 sumT,sumK,sumL,sumP,einK,mcoeff
       real*8 ldist,lldist
       real*8 PI,valueC,fval,res1
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 n1(n),n2(n),n3(n),n4(n)
       real*8 ymode(n),ylatent(n)
       real*8 p100(nmax),tempn100
       real*8 prior100,cnst1,cnst2
       real*8 num,den,a1,b1,pw
       real*8 tempval1,tempval2
       real*8 sval1,sval2
       real*8 aa,bb,ntilda,diff,tempcnt
c       
       integer l1,l2,l3
       integer sumy,sy12,sy23,sy234
       integer ids(n),nStudy(n),trt(n),collap(n)
       external gcoeff1,DICnorm1
c	
       common /vecx/x
       common /vecz/z  
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecmu/mu      
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3 
       common /vecLU/LU
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /dummy/index
       common /vecntilda/ntilda
c       
       PI=const('PI')
       sumy=0
       l1=0
       l2=0
       l3=0
       sumT=0.0d0
       sumL=0.0d0
       sumK=0.0d0
       sumP=0.0d0
       ldist=0.0d0
       lldist=0.0d0
       ntilda=0.0d0
       mcoeff=0.0d0
       aa=0.0d0
       bb=0.0d0
       diff=0.0d0
c       
       einK=0.0d0   
       fval=0.0d0
       res1=0.0d0
       num=0.0d0
       den=0.0d0
c
       do jj=1,n
        n1(jj)=0.0d0
        n2(jj)=0.0d0
        n3(jj)=0.0d0
        n4(jj)=0.0d0
       enddo
       do i1=1,nmax
        p100(i1)=0.0d0
       enddo

c......Variable initialization
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
       llik=0.0d0
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
       muK=0.0d0
       tRTnt=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       p1=0.0d0
       p2=0.0d0
       p3=0.0d0
       prior100=0.0d0
       cnst1=0.0d0
       cnst2=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       einK=0.0d0
       tempcnt=0.0d0
c
       id=index
       if (trt(id) .gt. 0) then
        tRTnt=delta(trt(id))
       else
        tRTnt=0.0d0
       endif
c 
       muK=mu(nStudy(id))
       zbeta=phi(1)+z(1,id)*phi(2)
       xbeta=beta(1)+x(1,id)*beta(2)
     1     +x(2,id)*beta(3)+x(3,id)*beta(4)
     2     +x(4,id)*beta(5)+x(5,id)*beta(6)
       if (id .eq. 35) then
        xbeta=x(1,id)*beta(2)
     1     +x(2,id)*beta(3)+x(3,id)*beta(4)
     2     +x(4,id)*beta(5)+x(5,id)*beta(6)
       endif
       eta2=alpha2(1)+x(1,id)*alpha2(2)
     1     +x(2,id)*alpha2(3)+x(3,id)*alpha2(4)
     2     +x(4,id)*alpha2(5)+x(5,id)*alpha2(6)
       eta3=alpha3(1)+x(1,id)*alpha3(2)
     1     +x(2,id)*alpha3(3)+x(3,id)*alpha3(4)
     2     +x(4,id)*alpha3(5)+x(5,id)*alpha3(6)
       c1=0.0d0
       c2=c1+dexp(eta2)
       c3=c2+dexp(eta3)
c	   
       ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
       ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
       ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)
c
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

c......Lower and upper bounds for latent counts
       aa=LU(1,id)
       bb=LU(2,id)

c......Initial setting 
       a1=1.0d0
       b1=1.0d0
       pw=1.0d0

c......To use normalized distribution
       do 50 kk=aa,bb
        temp=floor((aa+bb)/2.0d0)
c       temp=bb
        if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
         diff=(p3-p2)-(p2-p1)
         if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
          temp=floor((aa+bb)/2.0d0)
          if (kk .lt. temp) then
c......We need (p2+p3)/2, this is p3-p1                  
           lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
          else if (kk .gt. temp) then
           lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
          else
           tempval1=0.0d0
           tempval2=0.0d0
           sval1=0.0d0
           sval2=0.0d0
           do k1=aa,temp-1.0d0
            sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
            tempval1=tempval1+sval1
           enddo
           do k2=temp+1.0d0,bb
            sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
            tempval2=tempval2+sval2
           enddo
           lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
          endif
         else
          if ((p2-p1) .lt. (p3-p2)) then
           num=(p3-p2)**b1
           den=(p2-p1)**a1+(p3-p2)**b1
           if (kk .le. temp) then
            lldist=kk*pw*(dlog(num)-dlog(den))
           else
            lldist=0.0d0
           endif
          else if ((p2-p1) .gt. (p3-p2)) then
           num=(p2-p1)**a1
           den=(p2-p1)**a1+(p3-p2)**b1
           if (kk .gt. temp) then
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
           else
            lldist=0.0d0
           endif
          endif
         endif
        endif
        sumP=sumP+dexp(lldist)
 50    continue

c......Summing latent counts over all possible values
       do 100 kk=aa,bb
        mcoeff=0.0d0
        sumK=0.0d0
        sumL=0.0d0
        sumC=0.0d0
        ntilda=kk

c......Define cell counts
        if (collap(id) .eq. 0) then
         n1(id)=ss(id)-n70(id)
         n2(id)=n70(id)-n100(id)-(nrem(id)-ntilda)
         n3(id)=n100(id)-ntilda
         n4(id)=nrem(id)
        else if (collap(id) .eq. -1) then
         n1(id)=0.0d0
         n2(id)=ss(id)-n100(id)-nrem(id)+ntilda
         n3(id)=n100(id)-ntilda
         n4(id)=nrem(id)
        else if (collap(id) .eq. 1) then
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=n70(id)-nrem(id)
         n4(id)=nrem(id)
        else
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=0.0d0
         n4(id)=n70(id)
        endif

c......Define multinomial coefficients per latent count
        mcoeff=gcoeff1(id)

c......Define log-likelihood on each latent count
c......Fully observed categories		 
        if (collap(id) .eq. 0) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c         
          if (n2(id) .gt. 0.0d0) then
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
           llik2=(llik2a+llik2b+llik2c)*n2(id)
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 1 and 2 are collapsed
        else if (collap(id) .eq. -1) then
          sy12=n1(id)+n2(id)
          if (sy12 .gt. 0.0d0) then
           if (ein2 .le. 0.0d0) then
            llik2a=-dlog(1.0d0+dexp(ein2))
           else
            llik2a=-ein2-dlog(1.0d0+dexp(-ein2))
           endif
           llik2=llik2a*sy12
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif 

c         if (id .eq. 18) then
c          write(*,*) 'i=',id,'llik2=',llik2,'llik3=',llik3
c         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(id) .eq. 1) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy23=n2(id)+n3(id)
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
           llik2=(llik2a+llik2b+llik2c)*sy23
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 2,3,and 4 are collapsed
        else if (collap(id) .eq. 2) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy234=n2(id)+n3(id)+n4(id)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik4a=ein1-dlog(1.0d0+dexp(ein1))
           else
            llik4a=-dlog(1.0d0+dexp(-ein1))
           endif
           llik4=llik4a*sy234
          endif
        endif 

c......Customized construction
       temp=floor((aa+bb)/2.0d0)
c      temp=floor(aa+(bb-aa+1.0d0)*0.70d0)
c      temp=bb
       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
        diff=(p3-p2)-(p2-p1)
        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
         temp=floor((aa+bb)/2.0d0)
         if (kk .lt. temp) then
          lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
         else if (kk .gt. temp) then
          lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
         else
          tempval1=0.0d0
          tempval2=0.0d0
          sval1=0.0d0
          sval2=0.0d0
          do k1=aa,temp-1.0d0
           sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
           tempval1=tempval1+sval1
          enddo
          do k2=temp+1.0d0,bb
           sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
           tempval2=tempval2+sval2
          enddo
          lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
         endif
        else 
         if ((p2-p1) .lt. (p3-p2)) then
          num=(p3-p2)**b1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .le. temp) then
           lldist=kk*pw*(dlog(num)-dlog(den))
c          lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         else if ((p2-p1) .gt. (p3-p2)) then
          num=(p2-p1)**a1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .gt. temp) then               
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
c           lldist=kk*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         endif
        endif
       endif
   
c......Discrete Uniform distribution (P4)
c       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
c        lldist=-dlog(bb-aa+1.0d0)
c       else 
c        lldist=0.0d0
c       endif

c......Define a deviance function
       sumK=llik1+llik2+llik3+llik4
       sumC=DICnorm1(id)
       sumL=dexp(sumK+mcoeff+lldist-sumC)
       sumT=sumT+sumL
 100   continue
c
c      write(*,*) 'id=',id,'sumT=',sumT
       sumK=DICnorm1(id)
c      write(*,*) 'id=',id,'sumT=',sumT,'val=',sumK
       fval=dlog(sumT)
       res1=-2.0d0*fval
       return
       end


c......Compute normalizing constant for correct DIC
       real*8 function DICnorm1(id)
       implicit real*8 (a-h,o-z) 
       parameter (n=36,K=18,nT=13,nx=5,nz=5,nmax=500)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 phi(nz+1),delta(nT),beta(nx+1),mu(K)
       real*8 alpha2(nx+1),alpha3(nx+1)
       real*8 theta(K),prior
       real*8 llik,llik1,llik2,llik3,llik4  
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 c1,c2,c3,eta2,eta3
       real*8 p1,p2,p3
       real*8 xi,ein1,ein2,ein3
       real*8 sumT,sumK,sumL,sumP,einK,mcoeff
       real*8 ldist,lldist
       real*8 PI,valueC,fval,res1
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 n1(n),n2(n),n3(n),n4(n)
       real*8 ymode(n),ylatent(n)
       real*8 p100(nmax),tempn100
       real*8 prior100,cnst1,cnst2
       real*8 num,den,a1,b1,pw
       real*8 tempval1,tempval2
       real*8 aa,bb,ntilda,diff,tempcnt
c       
       integer l1,l2,l3
       integer sumy,sy12,sy23,sy234
       integer ids(n),nStudy(n),trt(n),collap(n)
       external gcoeff1
c	
       common /vecx/x
       common /vecz/z  
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecphi/phi
       common /vecdelta/delta
       common /vecmu/mu      
       common /vecbeta/beta
       common /vecalpha2/alpha2
       common /vecalpha3/alpha3 
       common /vecLU/LU
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /dummy/index
c
       sumy=0
       l1=0
       l2=0
       l3=0
       sumT=0.0d0
       sumP=0.0d0
       ldist=0.0d0
       lldist=0.0d0
       ntilda=0.0d0
       aa=0.0d0
       bb=0.0d0
       diff=0.0d0
c       
       einK=0.0d0   
       fval=0.0d0
       res1=0.0d0
       num=0.0d0
       den=0.0d0
c
       do jj=1,n
        n1(jj)=0.0d0
        n2(jj)=0.0d0
        n3(jj)=0.0d0
        n4(jj)=0.0d0
       enddo
       do i1=1,nmax
        p100(i1)=0.0d0
       enddo

c......Variable initialization
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
       llik=0.0d0
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
       muK=0.0d0
       tRTnt=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       p1=0.0d0
       p2=0.0d0
       p3=0.0d0
       prior100=0.0d0
       cnst1=0.0d0
       cnst2=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       einK=0.0d0
       tempcnt=0.0d0
c
       id=index
       if (trt(id) .gt. 0) then
        tRTnt=delta(trt(id))
       else
        tRTnt=0.0d0
       endif
c 
       muK=mu(nStudy(id))
       zbeta=phi(1)+z(1,id)*phi(2)
       xbeta=beta(1)+x(1,id)*beta(2)
     1     +x(2,id)*beta(3)+x(3,id)*beta(4)
     2     +x(4,id)*beta(5)+x(5,id)*beta(6)
       if (id .eq. 35) then
        xbeta=x(1,id)*beta(2)
     1     +x(2,id)*beta(3)+x(3,id)*beta(4)
     2     +x(4,id)*beta(5)+x(5,id)*beta(6)
       endif
       eta2=alpha2(1)+x(1,id)*alpha2(2)
     1     +x(2,id)*alpha2(3)+x(3,id)*alpha2(4)
     2     +x(4,id)*alpha2(5)+x(5,id)*alpha2(6)
       eta3=alpha3(1)+x(1,id)*alpha3(2)
     1     +x(2,id)*alpha3(3)+x(3,id)*alpha3(4)
     2     +x(4,id)*alpha3(5)+x(5,id)*alpha3(6)
       c1=0.0d0
       c2=c1+dexp(eta2)
       c3=c2+dexp(eta3)
c	   
       ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
       ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
       ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)
c
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

c......Lower and upper bounds for latent counts
       aa=LU(1,id)
       bb=LU(2,id)

c......Initial setting 
       a1=1.0d0
       b1=1.0d0
       pw=1.0d0

c......Summing latent counts over all possible values
       do 100 kk=aa,bb
        mcoeff=0.0d0
        sumK=0.0d0
        sumL=0.0d0
        sumC=0.0d0
        sumD=0.0d0

c......Define cell counts
        if (collap(id) .eq. 0) then
         n1(id)=ss(id)-n70(id)
         n2(id)=n70(id)-n100(id)-(nrem(id)-kk)
         n3(id)=n100(id)-kk
         n4(id)=nrem(id)
        else if (collap(id) .eq. -1) then
         n1(id)=0.0d0
         n2(id)=ss(id)-n100(id)-nrem(id)+kk
         n3(id)=n100(id)-kk
         n4(id)=nrem(id)
        else if (collap(id) .eq. 1) then
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=n70(id)-nrem(id)
         n4(id)=nrem(id)
        else
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=0.0d0
         n4(id)=n70(id)
        endif

c......Compute multinomial coefficients
        mcoeff=gcoeff1(id)

c......Compute likelihoods
        if (collap(id) .eq. 0) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c         
          if (n2(id) .gt. 0.0d0) then
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
           llik2=(llik2a+llik2b+llik2c)*n2(id)
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 1 and 2 are collapsed
        else if (collap(id) .eq. -1) then
          sy12=n1(id)+n2(id)
          if (sy12 .gt. 0.0d0) then
           if (ein2 .le. 0.0d0) then
            llik2a=-dlog(1.0d0+dexp(ein2))
           else
            llik2a=-ein2-dlog(1.0d0+dexp(-ein2))
           endif
           llik2=llik2a*sy12
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif 

c         if (id .eq. 18) then
c          write(*,*) 'i=',id,'llik2=',llik2,'llik3=',llik3
c         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(id) .eq. 1) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy23=n2(id)+n3(id)
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
           llik2=(llik2a+llik2b+llik2c)*sy23
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 2,3,and 4 are collapsed
        else if (collap(id) .eq. 2) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy234=n2(id)+n3(id)+n4(id)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik4a=ein1-dlog(1.0d0+dexp(ein1))
           else
            llik4a=-dlog(1.0d0+dexp(-ein1))
           endif
           llik4=llik4a*sy234
          endif
        endif

c......Compute distribution of the latent counts
       temp=floor((aa+bb)/2.0d0)
c      temp=bb
       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
        diff=(p3-p2)-(p2-p1)
        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
         temp=floor((aa+bb)/2.0d0)
         if (kk .lt. temp) then
          lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
         else if (kk .gt. temp) then
          lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
         else
          tempval1=0.0d0
          tempval2=0.0d0
          sval1=0.0d0
          sval2=0.0d0
          do k1=aa,temp-1.0d0
           sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
           tempval1=tempval1+sval1
          enddo
          do k2=temp+1.0d0,bb
           sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
           tempval2=tempval2+sval2
          enddo
          lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
         endif
        else 
         if ((p2-p1) .lt. (p3-p2)) then
          num=(p3-p2)**b1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .le. temp) then
           lldist=kk*pw*(dlog(num)-dlog(den))
c          lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         else if ((p2-p1) .gt. (p3-p2)) then
          num=(p2-p1)**a1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .gt. temp) then               
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
c           lldist=kk*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         endif
        endif
       endif  

c......Discrete uniform distribution (P4)
c       if ((collap(id) .eq. -1) .or. (collap(id) .eq. 0)) then
c        lldist=-dlog(bb-aa+1.0d0)
c       else
c        lldist=0.0d0
c       endif

c......Compute the sum for normalization
       sumC=mcoeff
       sumL=llik1+llik2+llik3+llik4
       sumD=lldist
       sumK=sumC+sumL+sumD
       sumT=sumT+sumK
 100   continue
c
c      write(*,*) 'id=',id,'sumT=',sumT
       DICnorm1=dexp(sumT)
       return
       end   



c......Calculation DICbar: Goodness-of-fit
       subroutine DIC2(xi,res2) 
       parameter (n=36,K=18,nT=13,nx=5,nz=5,nmax=500)
       real*8 n1(n),n2(n),n3(n),n4(n)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 ephi(nz+1),edelta(nT),ebeta(nx+1)
       real*8 ealpha2(nx+1),ealpha3(nx+1),emu(K)
       real*8 etheta(K),prior,einK
       real*8 llik,llik1,llik2,llik3,llik4
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 tRTnt,emuK,xbeta,zbeta
       real*8 c1,c2,c3,eta2,eta3
       real*8 xi,ein1,ein2,ein3
       real*8 sumT,sumL,sumK,sumP
       real*8 PI,valueC,fval,res2
       real*8 temp1,temp2
       real*8 p100(nmax),tempn100
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 cnst1,cnst2,prior100
       real*8 p1,p2,p3,num,den
       real*8 a1,b1,pw
       real*8 ldist,lldist,mcoeff
       real*8 tempval1,tempval2
       real*8 aa,bb,ntilda,diff,tempcnt
c       
       integer sumy,sy12,sy23,sy234
       integer l1,l2,l3
       integer ids(n),nStudy(n),trt(n),collap(n)
       external gcoeff1,DICnorm2
c	
       common /vecx/x
       common /vecz/z  
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecephi/ephi
       common /vecedelta/edelta
       common /vecemu/emu
       common /vecebeta/ebeta
       common /vecealpha2/ealpha2
       common /vecealpha3/ealpha3
       common /vecLU/LU
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /dummy/index  
       common /vecntilda/ntilda
c
       PI=const('PI')
       sumy=0
       sumT=0.0d0
       sumL=0.0d0
       sumK=0.0d0
       sumP=0.0d0
       fval=0.0d0
       res2=0.0d0
       tempn100=0.0d0
       l1=0
       l2=0
       l3=0
       mcoeff=0.0d0
       cnst1=0.0d0
       cnst2=0.0d0
       p1=0.0d0
       p2=0.0d0
       p3=0.0d0
       num=0.0d0
       den=0.0d0
c
       aa=0.0d0
       bb=0.0d0       
       ldist=0.0d0
       lldist=0.0d0
       diff=0.0d0
c
       sy12=0
       sy23=0
       sy234=0
       llik1=0.0d0
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
       emuK=0.0d0
       tRTnt=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       einK=0.0d0
       temp1=0.0d0
       temp2=0.0d0
       tempcnt=0.0d0
c	   
       do jj=1,n
        n1(jj)=0.0d0
        n2(jj)=0.0d0
        n3(jj)=0.0d0
        n4(jj)=0.0d0
       enddo
       do i1=1,nmax
        p100(i1)=0.0d0
       enddo
c
       id=index
       if (trt(id) .gt. 0) then
        tRTnt=edelta(trt(id))
       else
        tRTnt=0.0d0
       endif
c 
       emuK=emu(nStudy(id))
       zbeta=ephi(1)+z(1,id)*ephi(2)
       xbeta=ebeta(1)+x(1,id)*ebeta(2)
     1       +x(2,id)*ebeta(3)+x(3,id)*ebeta(4)
     2       +x(4,id)*ebeta(5)+x(5,id)*ebeta(6)
       if (id .eq. 35) then
        xbeta=x(1,id)*ebeta(2)
     1       +x(2,id)*ebeta(3)+x(3,id)*ebeta(4)
     2       +x(4,id)*ebeta(5)+x(5,id)*ebeta(6)
       endif
       eta2=ealpha2(1)+x(1,id)*ealpha2(2)
     1       +x(2,id)*ealpha2(3)+x(3,id)*ealpha2(4)
     2       +x(4,id)*ealpha2(5)+x(5,id)*ealpha2(6)
       eta3=ealpha3(1)+x(1,id)*ealpha3(2)
     1       +x(2,id)*ealpha3(3)+x(3,id)*ealpha3(4)
     2       +x(4,id)*ealpha3(5)+x(5,id)*ealpha3(6)
       c1=0.0d0
       c2=c1+dexp(eta2)
       c3=c2+dexp(eta3)
c		
       ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*emuK)
       ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*emuK)
       ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*emuK)
c
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

       aa=LU(1,id)
       bb=LU(2,id)

c......Initial setting
       a1=1.0d0
       b1=1.0d0
       pw=1.0d0

c......To use normalized distribution
       do kk=aa,bb
        temp=floor((aa+bb)/2.0d0)
c       temp=bb
        if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
         diff=(p3-p2)-(p2-p1)
         if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
          temp=floor((aa+bb)/2.0d0)
          if (kk .lt. temp) then
           lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
          else if (kk .gt. temp) then
           lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
          else
           tempval1=0.0d0
           tempval2=0.0d0
           sval1=0.0d0
           sval2=0.0d0
           do k1=aa,temp-1.0d0
            sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
            tempval1=tempval1+sval1
           enddo
           do k2=temp+1.0d0,bb
            sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
            tempval2=tempval2+sval2
           enddo
           lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
          endif
         else
          if ((p2-p1) .lt. (p3-p2)) then
           num=(p3-p2)**b1
           den=(p2-p1)**a1+(p3-p2)**b1
           if (kk .le. temp) then
            lldist=kk*pw*(dlog(num)-dlog(den))
           else
            lldist=0.0d0
           endif
          else if ((p2-p1) .gt. (p3-p2)) then
           num=(p2-p1)**a1
           den=(p2-p1)**a1+(p3-p2)**b1
           if (kk .gt. temp) then
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
           else
            lldist=0.0d0
           endif
          endif
         endif
        endif
        sumP=sumP+dexp(lldist)
       enddo

c......Summing over latent counts
       do 200 kk=aa,bb
        mcoeff=0.0d0
        sumL=0.0d0
        sumK=0.0d0
        ntilda=kk

c......Define cell counts
        if (collap(id) .eq. 0) then
         n1(id)=ss(id)-n70(id)
         n2(id)=n70(id)-n100(id)-(nrem(id)-ntilda)
         n3(id)=n100(id)-ntilda
         n4(id)=nrem(id)
        else if (collap(id) .eq. -1) then
         n1(id)=0.0d0
         n2(id)=ss(id)-n100(id)-nrem(id)+ntilda
         n3(id)=n100(id)-ntilda
         n4(id)=nrem(id)
        else if (collap(id) .eq. 1) then
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=n70(id)-nrem(id)
         n4(id)=nrem(id)
        else
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=0.0d0
         n4(id)=n70(id)
        endif

c......Define multinomial coefficients per latent count
        mcoeff=gcoeff1(id)

c......Define log-likelihood on each latent count
c......Fully observed categories		 
        if (collap(id) .eq. 0) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c         
          if (n2(id) .gt. 0.0d0) then
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
           llik2=(llik2a+llik2b+llik2c)*n2(id)
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 1 and 2 are collapsed
        else if (collap(id) .eq. -1) then
          sy12=n1(id)+n2(id)
          if (sy12 .gt. 0.0d0) then
           if (ein2 .le. 0.0d0) then
            llik2a=-dlog(1.0d0+dexp(ein2))
           else
            llik2a=-ein2-dlog(1.0d0+dexp(-ein2))
           endif
           llik2=llik2a*sy12
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif 

c......Partially observed 2 and 3 are collapsed
        else if (collap(id) .eq. 1) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy23=n2(id)+n3(id)
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
           llik2=(llik2a+llik2b+llik2c)*sy23
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 2,3,and 4 are collapsed
        else if (collap(id) .eq. 2) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy234=n2(id)+n3(id)+n4(id)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik4a=ein1-dlog(1.0d0+dexp(ein1))
           else
            llik4a=-dlog(1.0d0+dexp(-ein1))
           endif
           llik4=llik4a*sy234
          endif
        endif 

c......Define a distribution of the latent counts
c......Define a distribution of the latent counts
       temp=floor((aa+bb)/2.0d0)
c      temp=floor(aa+(bb-aa+1.0d0)*0.70d0)
c      temp=bb
       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
        diff=(p3-p2)-(p2-p1)
        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
         temp=floor((aa+bb)/2.0d0)
         if (kk .lt. temp) then
          lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
         else if (kk .gt. temp) then
          lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
         else
          tempval1=0.0d0
          tempval2=0.0d0
          sval1=0.0d0
          sval2=0.0d0
          do k1=aa,temp-1.0d0
           sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
           tempval1=tempval1+sval1
          enddo
          do k2=temp+1.0d0,bb
           sval2=(k2-temp1)*dlog((p3-p1)/2.0d0)
           tempval2=tempval2+sval2
          enddo
          lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
         endif
        else 
         if ((p2-p1) .lt. (p3-p2)) then
          num=(p3-p2)**b1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .le. temp) then
           lldist=kk*pw*(dlog(num)-dlog(den))
c          lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         else if ((p2-p1) .gt. (p3-p2)) then
          num=(p2-p1)**a1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .gt. temp) then               
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
c           lldist=kk*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         endif
        endif
       endif
c		
c......Discrete Uniform distribution (P4)
c       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
c        lldist=-dlog(bb-aa+1.0d0)
c       else
c        lldist=0.0d0
c       endif

c......Define a deviance function
       sumK=llik1+llik2+llik3+llik4
       sumC=DICnorm2(id)
       sumL=dexp(sumK+mcoeff+lldist-sumC)
       sumT=sumT+sumL
 200   continue   
c      
       fval=dlog(sumT)
       res2=-2.0d0*fval 
       return
       end


c......Compute normalizing constant for correct DIC
       real*8 function DICnorm2(id)
       implicit real*8 (a-h,o-z) 
       parameter (n=36,K=18,nT=13,nx=5,nz=5,nmax=500)
       real*8 x(nx,n),z(nz,n),LU(2,n)
       real*8 ephi(nz+1),edelta(nT),ebeta(nx+1),emu(K)
       real*8 ealpha2(nx+1),ealpha3(nx+1)
       real*8 theta(K),prior
       real*8 llik,llik1,llik2,llik3,llik4  
       real*8 llik1a,llik2a,llik2b,llik2c
       real*8 llik3a,llik3b,llik3c,llik4a
       real*8 tRTnt,muK,xbeta,zbeta
       real*8 c1,c2,c3,eta2,eta3
       real*8 p1,p2,p3
       real*8 xi,ein1,ein2,ein3
       real*8 sumT,sumK,sumL,sumP,einK,mcoeff
       real*8 ldist,lldist
       real*8 PI,valueC,fval,res1
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 n1(n),n2(n),n3(n),n4(n)
       real*8 ymode(n),ylatent(n)
       real*8 p100(nmax),tempn100
       real*8 prior100,cnst1,cnst2
       real*8 num,den,a1,b1,pw
       real*8 tempval1,tempval2
       real*8 aa,bb,ntilda,diff,tempcnt
c       
       integer l1,l2,l3
       integer sumy,sy12,sy23,sy234
       integer ids(n),nStudy(n),trt(n),collap(n)
       external gcoeff1
c	
       common /vecx/x
       common /vecz/z  
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecephi/ephi
       common /vecedelta/edelta
       common /vecemu/emu
       common /vecebeta/ebeta
       common /vecealpha2/ealpha2
       common /vecealpha3/ealpha3 
       common /vecLU/LU
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /dummy/index
c
       sumy=0
       l1=0
       l2=0
       l3=0
       sumT=0.0d0
       sumP=0.0d0
       ldist=0.0d0
       lldist=0.0d0
       ntilda=0.0d0
       aa=0.0d0
       bb=0.0d0
       diff=0.0d0
c       
       einK=0.0d0   
       fval=0.0d0
       res1=0.0d0
       num=0.0d0
       den=0.0d0
c
       do jj=1,n
        n1(jj)=0.0d0
        n2(jj)=0.0d0
        n3(jj)=0.0d0
        n4(jj)=0.0d0
       enddo
       do i1=1,nmax
        p100(i1)=0.0d0
       enddo

c......Variable initialization
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
       llik=0.0d0
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
       muK=0.0d0
       tRTnt=0.0d0
       eta2=0.0d0 
       eta3=0.0d0
       c1=0.0d0
       c2=0.0d0
       c3=0.0d0
       p1=0.0d0
       p2=0.0d0
       p3=0.0d0
       prior100=0.0d0
       cnst1=0.0d0
       cnst2=0.0d0
       ein1=0.0d0
       ein2=0.0d0
       ein3=0.0d0
       einK=0.0d0
       tempcnt=0.0d0
c
       id=index
       if (trt(id) .gt. 0) then
        tRTnt=edelta(trt(id))
       else
        tRTnt=0.0d0
       endif
c 
       muK=emu(nStudy(id))
       zbeta=ephi(1)+z(1,id)*ephi(2)
       xbeta=ebeta(1)+x(1,id)*ebeta(2)
     1     +x(2,id)*ebeta(3)+x(3,id)*ebeta(4)
     2     +x(4,id)*ebeta(5)+x(5,id)*ebeta(6)
       if (id .eq. 35) then
        xbeta=x(1,id)*ebeta(2)
     1     +x(2,id)*ebeta(3)+x(3,id)*ebeta(4)
     2     +x(4,id)*ebeta(5)+x(5,id)*ebeta(6)
       endif
       eta2=ealpha2(1)+x(1,id)*ealpha2(2)
     1     +x(2,id)*ealpha2(3)+x(3,id)*ealpha2(4)
     2     +x(4,id)*ealpha2(5)+x(5,id)*ealpha2(6)
       eta3=ealpha3(1)+x(1,id)*ealpha3(2)
     1     +x(2,id)*ealpha3(3)+x(3,id)*ealpha3(4)
     2     +x(4,id)*ealpha3(5)+x(5,id)*ealpha3(6)
       c1=0.0d0
       c2=c1+dexp(eta2)
       c3=c2+dexp(eta3)
c	   
       ein1=-(c1+tRTnt+xbeta+dexp(zbeta)*muK)
       ein2=-(c2+tRTnt+xbeta+dexp(zbeta)*muK)
       ein3=-(c3+tRTnt+xbeta+dexp(zbeta)*muK)
c
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

c......Lower and upper bounds for latent counts
       aa=LU(1,id)
       bb=LU(2,id)

c......Initial setting 
       a1=1.0d0
       b1=1.0d0
       pw=1.0d0

c......Summing latent counts over all possible values
       do 100 kk=aa,bb
        mcoeff=0.0d0
        sumK=0.0d0
        sumL=0.0d0
        sumC=0.0d0
        sumD=0.0d0

c......Define cell counts
        if (collap(id) .eq. 0) then
         n1(id)=ss(id)-n70(id)
         n2(id)=n70(id)-n100(id)-(nrem(id)-kk)
         n3(id)=n100(id)-kk
         n4(id)=nrem(id)
        else if (collap(id) .eq. -1) then
         n1(id)=0.0d0
         n2(id)=ss(id)-n100(id)-nrem(id)+kk
         n3(id)=n100(id)-kk
         n4(id)=nrem(id)
        else if (collap(id) .eq. 1) then
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=n70(id)-nrem(id)
         n4(id)=nrem(id)
        else
         n1(id)=ss(id)-n70(id)
         n2(id)=0.0d0
         n3(id)=0.0d0
         n4(id)=n70(id)
        endif

c......Compute multinomial coefficients
        mcoeff=gcoeff1(id)

c......Compute likelihoods
        if (collap(id) .eq. 0) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c         
          if (n2(id) .gt. 0.0d0) then
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
           llik2=(llik2a+llik2b+llik2c)*n2(id)
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 1 and 2 are collapsed
        else if (collap(id) .eq. -1) then
          sy12=n1(id)+n2(id)
          if (sy12 .gt. 0.0d0) then
           if (ein2 .le. 0.0d0) then
            llik2a=-dlog(1.0d0+dexp(ein2))
           else
            llik2a=-ein2-dlog(1.0d0+dexp(-ein2))
           endif
           llik2=llik2a*sy12
          endif
c
          if (n3(id) .gt. 0.0d0) then
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
           llik3=(llik3a+llik3b+llik3c)*n3(id)
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif 

c         if (id .eq. 18) then
c          write(*,*) 'i=',id,'llik2=',llik2,'llik3=',llik3
c         endif

c......Partially observed 2 and 3 are collapsed
        else if (collap(id) .eq. 1) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy23=n2(id)+n3(id)
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
           llik2=(llik2a+llik2b+llik2c)*sy23
          endif
c
          if (n4(id) .gt. 0.0d0) then
           if (ein3 .le. 0.0d0) then
            llik4a=ein3-dlog(1.0d0+dexp(ein3))
           else
            llik4a=-dlog(1.0d0+dexp(-ein3))
           endif
           llik4=llik4a*n4(id)
          endif
   
c......Partially observed 2,3,and 4 are collapsed
        else if (collap(id) .eq. 2) then
          if (n1(id) .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik1a=-dlog(1.0d0+dexp(ein1))
           else
            llik1a=-ein1-dlog(1.0d0+dexp(-ein1))
           endif
           llik1=llik1a*n1(id)
          endif
c
          sy234=n2(id)+n3(id)+n4(id)
          if (sy234 .gt. 0.0d0) then
           if (ein1 .le. 0.0d0) then
            llik4a=ein1-dlog(1.0d0+dexp(ein1))
           else
            llik4a=-dlog(1.0d0+dexp(-ein1))
           endif
           llik4=llik4a*sy234
          endif
        endif

c......Compute distribution of the latent counts
       temp=floor((aa+bb)/2.0d0)
c      temp=bb
       if ((collap(id) .eq. 0) .or. (collap(id) .eq. -1)) then
        diff=(p3-p2)-(p2-p1)
        if ((diff .ge. -0.05d0) .and. (diff .le. 0.05d0)) then
         temp=floor((aa+bb)/2.0d0)
         if (kk .lt. temp) then
          lldist=(temp-kk)*dlog((p3-p1)/2.0d0)
         else if (kk .gt. temp) then
          lldist=(kk-temp)*dlog((p3-p1)/2.0d0)
         else
          tempval1=0.0d0
          tempval2=0.0d0
          sval1=0.0d0
          sval2=0.0d0
          do k1=aa,temp-1.0d0
           sval1=(temp-k1)*dlog((p3-p1)/2.0d0)
           tempval1=tempval1+sval1
          enddo
          do k2=temp+1.0d0,bb
           sval2=(k2-temp)*dlog((p3-p1)/2.0d0)
           tempval2=tempval2+sval2
          enddo
          lldist=dlog(1.0d0-dexp(tempval1)-dexp(tempval2))
         endif
        else 
         if ((p2-p1) .lt. (p3-p2)) then
          num=(p3-p2)**b1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .le. temp) then
           lldist=kk*pw*(dlog(num)-dlog(den))
c          lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         else if ((p2-p1) .gt. (p3-p2)) then
          num=(p2-p1)**a1
          if (num .lt. 0.001d0) num=0.001d0
          den=(p2-p1)**a1+(p3-p2)**b1
          if (den .lt. 0.001d0) den=0.001d0
          if (kk .gt. temp) then               
            lldist=(bb-kk)*pw*(dlog(num)-dlog(den))
c           lldist=kk*pw*(dlog(num)-dlog(den))
          else
           lldist=0.0d0
          endif
         endif
        endif
       endif  

c......Discrete uniform distribution (P4)
c       if ((collap(id) .eq. -1) .or. (collap(id) .eq. 0)) then
c        lldist=-dlog(bb-aa+1.0d0)
c       else
c        lldist=0.0d0
c       endif

c......Compute the sum for normalization
       sumC=mcoeff
       sumL=llik1+llik2+llik3+llik4
       sumD=lldist
       sumK=sumC+sumL+sumD
       sumT=sumT+sumK
 100   continue
c
c      write(*,*) 'id=',id,'sumT=',sumT
       DICnorm2=dexp(sumT)
       return
       end   

c......Multinomial coefficient for barDIC
       real*8 function gcoeff1(a1)
       implicit real*8 (a-h,o-z)
       parameter (n=36,K=18,nT=13,nx=5,nz=5)
       real*8 n1(n),n2(n),n3(n),n4(n)
       real*8 ein1,ein2,ein3
       real*8 cnst1,cnst2,cnst3,cnst4,cnst5
       real*8 temp1,temp2,temp3,temp4,temp5
       real*8 emp1,emp2,emp3,emp4,emp5
       real*8 sumy,sumC,coef,comp,einK
       real*8 sy12,sy23,sy234
       real*8 ss(n),n70(n),n100(n),nrem(n)
       real*8 a1,a2,ntilda,tempcnt
c
       integer ind1,ind2
       integer id,nlatent
       integer ids(n),nStudy(n),trt(n),collap(n)
c	
       common /vecids/ids
       common /vecnStudy/nStudy
       common /vectrt/trt
       common /veccollap/collap
       common /vecss/ss
       common /vecn70/n70
       common /vecn100/n100
       common /vecnrem/nrem
       common /dummy/index
       common /vecntilda/ntilda
c   
       cnst1=0.0d0
       cnst2=0.0d0
       cnst3=0.0d0
       cnst4=0.0d0
       cnst5=0.0d0
       temp1=0.0d0
       temp2=0.0d0
       temp3=0.0d0
       temp4=0.0d0
       temp5=0.0d0
       emp1=0.0d0
       emp2=0.0d0
       emp3=0.0d0
       emp4=0.0d0
       emp5=0.0d0
       sumy=0.0d0
       sumC=0.0d0
       coef=0.0d0
       einK=0.0d0
       sy12=0.0d0
       sy23=0.0d0
       sy234=0.0d0
       tempcnt=0.0d0

       do jj=1,n
        n1(jj)=0.0d0
        n2(jj)=0.0d0
        n3(jj)=0.0d0
        n4(jj)=0.0d0
       enddo

c......Compute multinomial coefficients
       id=index
       nlatent=ntilda
       if (collap(id) .eq. 0) then
        n1(id)=ss(id)-n70(id)
        n2(id)=n70(id)-n100(id)-(nrem(id)-nlatent)
        n3(id)=n100(id)-nlatent
        n4(id)=nrem(id)
       else if (collap(id) .eq. -1) then
        n1(id)=0.0d0
        n2(id)=ss(id)-n100(id)-nrem(id)+nlatent
        n3(id)=n100(id)-nlatent
        n4(id)=nrem(id)
       else if (collap(id) .eq. 1) then
        n1(id)=ss(id)-n70(id)
        n2(id)=0.0d0
        n3(id)=n70(id)-nrem(id)
        n4(id)=nrem(id)
       else 
        n1(id)=ss(id)-n70(id)
        n2(id)=0.0d0
        n3(id)=0.0d0
        n4(id)=n70(id)
       endif

c......Fully observed categories
       sumy=n1(id)+n2(id)+n3(id)+n4(id)
       if (collap(id) .eq. 0) then
         if (n1(id) .gt. 0.0d0) then
           emp1=n1(id)
           do j1=1,emp1
            temp1=real(j1)
            cnst1=cnst1+dlog(temp1)
           enddo
         endif
c
         if (n2(id) .gt. 0.0d0) then
           emp2=n2(id)
           do j2=1,emp2
            temp2=real(j2)
            cnst2=cnst2+dlog(temp2)
           enddo
         endif
c
         if (n3(id) .gt. 0.0d0) then
           emp3=n3(id)
           do j3=1,emp3
            temp3=real(j3)
            cnst3=cnst3+dlog(temp3)
           enddo
         endif
c
         if (n4(id) .gt. 0.0d0) then
           emp4=n4(id)
           do j4=1,emp4
            temp4=real(j4)
            cnst4=cnst4+dlog(temp4)
           enddo
         endif
c   
         do j5=1,sumy
           temp5=real(j5)
           cnst5=cnst5+dlog(temp5)
         enddo

c......Partially observed 1 and 2 are collapsed
       else if (collap(id) .eq. -1) then
         sy12=n1(id)+n2(id)
         do j1=1,sy12
          temp1=real(j1)
          cnst1=cnst1+dlog(temp1)
         enddo   
         cnst2=0.0d0
c
         if (n3(id) .gt. 0.0d0) then
          emp3=n3(id)
          do j3=1,emp3
           temp3=real(j3)
           cnst3=cnst3+dlog(temp3)
          enddo
         endif
c
         if (n4(id) .gt. 0.0d0) then
          emp4=n4(id)
          do j4=1,emp4
            temp4=real(j4)
            cnst4=cnst4+dlog(temp4)
          enddo
         endif
c
         do j5=1,sumy
          temp5=real(j5)
          cnst5=cnst5+dlog(temp5)
         enddo

c......Partially observed 2 and 3 are collapsed
       else if (collap(id) .eq. 1) then
         if (n1(id) .gt. 0.0d0) then
           emp1=n1(id)
           do j1=1,emp1
            temp1=real(j1)
            cnst1=cnst1+dlog(temp1)
           enddo
         endif
c
         sy23=n2(id)+n3(id)
         do j2=1,sy23
           temp2=real(j2)
           cnst2=cnst2+dlog(temp2)
         enddo
         cnst3=0.0d0
c           
         if (n4(id) .gt. 0.0d0) then
           emp4=n4(id)
           do j4=1,emp4
            temp4=real(j4)
            cnst4=cnst4+dlog(temp4)
           enddo
         endif
c           
         do j5=1,sumy
           temp5=real(j5)
           cnst5=cnst5+dlog(temp5)
         enddo

c......Partially observed 2,3,and 4 are collapsed
       else if (collap(id) .eq. 2) then
         if (n1(id) .gt. 0.0d0) then
           emp1=n1(id)
           do j1=1,emp1
            temp1=real(j1)
            cnst1=cnst1+dlog(temp1)
           enddo
         endif
c         
         sy234=n2(id)+n3(id)+n4(id)
         do j2=1,sy234
           temp2=real(j2)
           cnst2=cnst2+dlog(temp2)
         enddo
         cnst3=0.0d0
         cnst4=0.0d0
c		  
         do j5=1,sumy
           temp5=real(j5)
           cnst5=cnst5+dlog(temp5)
         enddo
       endif
c
       coef=cnst5-(cnst1+cnst2+cnst3+cnst4)
       gcoeff1=coef
       return
       end

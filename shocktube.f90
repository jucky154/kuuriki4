program shocktube
  implicit none
  !for grid
  integer::mx
  real::dx,x(-5:505),xmin,xmax
  !cmpcnd
  real::cfl,dt,time,eps,g,roul,ul,pl,rour,ur,pr,rgas,cvgas,ecp
  real::atmpa,platm,pratm,tl,tr
  real::ustar
  integer::nlast,n,i
  real::now
  !sulute
  real::p12
  real::rou(-5:505),u(-5:505),p(-5:505),e(-5:505)
  real::roue(-5:505),ue(-5:505),pe(-5:505),ee(-5:505)
  !solver
  real::q(3,-5:505),qold(3,-5:505),flux(3,-5:505)

  !set grid
  mx=401
  xmin=-2.0
  xmax=2.0
  dx=(xmax-xmin)/float(mx-1)
  do i=-2,mx+3
     x(i)=xmin+dx*float(i-1)
  end do

  atmpa=1.013e+05
  eps=1.0e-06
  g=1.4
  rgas=8.314e+03/28.96
  cvgas=rgas/(g-1.0)
  
  ul=0
  platm=5
  pl=platm*atmpa
  tl=300
  roul=pl/(rgas*tl)

  ! you must choose pr < pl
  ur=0
  pratm=1
  pr=pratm*atmpa
  tr=300
  rour=pr/(rgas*tr)

  time=0.001
  cfl=0.5
  ecp=0.125
  
  ustar=sqrt(g*rgas*tl)

  dt=cfl*dx/sqrt(g*rgas*amax1(tl,tr))
  nlast=int(time/dt)
  time=dt*float(nlast)
  
  call calc_p12(g,tr,tl,pr,pl,p12)
  print*,p12
  
  do i=1,mx
     if(x(i).lt.0.0)then
        rou(i)=roul
        u(i)=ul
        p(i)=pl
        e(i)=p(i)/((g-1.0)*rou(i))
     else
        rou(i)=rour
        u(i)=ur
        p(i)=pr
        e(i)=p(i)/((g-1.0)*rou(i))
     end if
     q(1,i)=rou(i)
     q(2,i)=rou(i)*u(i)
     q(3,i)=rou(i)*(e(i)+0.5*u(i)*u(i))
  end do

  do n=1,nlast
     now=float(n)*dt
     do i=1,mx
        qold(1,i)=q(1,i)
        qold(2,i)=q(2,i)
        qold(3,i)=q(3,i)
     end do

     call calflx(g,mx,q,ecp,flux)
     
     do i=4,mx-3
        q(1,i)=qold(1,i)-(dt/dx)*(flux(1,i)-flux(1,i-1))
        q(2,i)=qold(2,i)-(dt/dx)*(flux(2,i)-flux(2,i-1))
        q(3,i)=qold(3,i)-(dt/dx)*(flux(3,i)-flux(3,i-1))
        rou(i)=q(1,i)
        u(i)=q(2,i)/rou(i)
        p(i)=(g-1.0)*(q(3,i)-0.5*rou(i)*u(i)*u(i))
        e(i)=p(i)/((g-1.0)*rou(i))
     end do
     
     call exact(x,roue,ue,pe,ee,mx,roul,rour,pr,pl,p12,g,rgas,tr,tl,ur,ul,now)
     do i=1,mx
        print*,x(i),now,p(i)/atmpa,pe(i)/atmpa
        !print*,x(i),now,e(i)/cvgas,ee(i)/cvgas
     end do
     print*," " !ここはtのほうをつかわないときは消すこと
  end do
  

  do i=1,mx
     !print*,x(i),rou(i),roue(i)
     !print*,x(i),u(i),ue(i)
     !print*,x(i),p(i)/atmpa,pe(i)/atmpa
     !print*,x(i),e(i)/cvgas,ee(i)/cvgas
  end do
  do i=1,mx
  
  end do
  stop
     
contains
  subroutine calflx(g,mx,q,ecp,flux)
    implicit none
    integer,intent(in)::mx
    real,intent(in)::q(3,-5:505),g,ecp
    real,intent(out)::flux(3,-5:505)
    real::q1lt,q1rt,q2lt,q2rt,q3lt,q3rt
    real::dq1,dq2,dq3,rlt,ult,plt,hlt
    real::rrt,urt,prt,hrt
    real::ubar,hbar,abar
    real::ram(3,-5:505),veck(3,3,-5:505),alfa(3,-5:505)
    real::b1,b2,ecpx,ff1,ff2,sss,qqq,ph1,ph2,ph3
    real::tvd1,tvd2,tvd3
    real::f1lt,f2lt,f3lt,f1rt,f2rt,f3rt
    integer::i

    do i=1,mx-1
       q1lt=q(1,i)
       q1rt=q(1,i+1)
       q2lt=q(2,i)
       q2rt=q(2,i+1)
       q3lt=q(3,i)
       q3rt=q(3,i+1)

       dq1=q(1,i+1)-q(1,i)
       dq2=q(2,i+1)-q(2,i)
       dq3=q(3,i+1)-q(3,i)

       rlt=q1lt
       ult=q2lt/rlt
       plt=(g-1.0)*(q3lt-0.5*rlt*ult**2)
       hlt=(q3lt+plt)/rlt
       rrt=q1rt
       urt=q2rt/q1rt
       prt=(g-1.0)*(q3rt-0.5*rrt*urt**2)
       hrt=(q3rt+prt)/rrt
       
       ubar=(sqrt(rlt)*ult+sqrt(rrt)*urt)/(sqrt(rlt)+sqrt(rrt))
       hbar=(sqrt(rlt)*hlt+sqrt(rrt)*hrt)/(sqrt(rlt)+sqrt(rrt))

       abar=sqrt( amax1( (g-1.0)*(hbar-0.5*ubar**2) &
            &          , amin1(g*plt/rlt, g*prt/rrt) ) )

       ram(1,i)=ubar-abar
       veck(1,1,i)=1.0
       veck(2,1,i)=ubar-abar
       veck(3,1,i)=hbar-ubar*abar
       
       ram(2,i)=ubar
       veck(1,2,i)=1.0
       veck(2,2,i)=ubar
       veck(3,2,i)=0.5*ubar**2
       
       ram(3,i)=ubar+abar
       veck(1,3,i)=1.0
       veck(2,3,i)=ubar+abar
       veck(3,3,i)=hbar+ubar*abar

       b2=(g-1.0)/abar**2
       b1=b2*ubar**2/2.0
       
       alfa(1,i)=0.5*(b1+ubar/abar)*dq1+0.5*(-b2*ubar-1.0/abar)*dq2+0.5*b2*dq3
       alfa(2,i)=(1.0-b1)*dq1+b2*ubar*dq2-b2*dq3
       alfa(3,i)=0.5*(b1-ubar/abar)*dq1+0.5*(-b2*ubar+1.0/abar)*dq2+0.5*b2*dq3

    end do

    do i=3,mx-3
       ecpx=ecp*amax1(abs(ram(1,i)),abs(ram(2,i)),abs(ram(3,i)))
       
       ff1=(dt/dx)*ram(1,i)**2
       ff2=abs(ram(1,i))

       if(ff2.lt.ecpx) ff2=(ff2**2+ecpx**2)*0.5/ecpx
       sss =sign(1.0, alfa(1,i))
       qqq =sss*amax1(0.0, amin1(sss*2.0*alfa(1,i-1)&
            & ,sss*2.0*alfa(1,i ) &
            & ,sss*2.0*alfa(1,i+1)&
            & ,sss*0.5*(alfa(1,i-1)+alfa(1,i+1))))
       ph1 =-ff1*qqq-ff2*(alfa(1,i)-qqq)

       ff1=(dt/dx)*ram(2,i)**2
       ff2=abs(ram(2,i))
       if(ff2.lt.ecpx) ff2=(ff2**2+ecpx**2)*0.5/ecpx
       sss =sign(1.0, alfa(2,i))
       qqq =sss*amax1(0.0, amin1(sss*2.0*alfa(2,i-1)&
            & , sss*2.0*alfa(2,i)&
            & , sss*2.0*alfa(2,i+1)&
            & , sss*0.5*(alfa(2,i-1)+alfa(2,i+1))))
       ph2 =-ff1*qqq-ff2*(alfa(2,i)-qqq)

       ff1=(dt/dx)*ram(3,i)**2
       ff2=abs(ram(3,i))
       if(ff2.lt.ecpx) ff2=(ff2**2+ecpx**2)*0.5/ecpx
       sss =sign(1.0, alfa(3,i))
       qqq =sss*amax1(0.0, amin1(sss*2.0*alfa(3,i-1)&
            &, sss*2.0*alfa(3,i)&
            &, sss*2.0*alfa(3,i+1)&
            &, sss*0.5*(alfa(3,i-1)+alfa(3,i+1))))
       ph3 =-ff1*qqq-ff2*(alfa(3,i)-qqq)

       tvd1=veck(1,1,i)*ph1+veck(1,2,i)*ph2+veck(1,3,i)*ph3
       tvd2=veck(2,1,i)*ph1+veck(2,2,i)*ph2+veck(2,3,i)*ph3
       tvd3=veck(3,1,i)*ph1+veck(3,2,i)*ph2+veck(3,3,i)*ph3
       q1lt=q(1,i )
       q1rt=q(1,i+1)
       q2lt=q(2,i )
       q2rt=q(2,i+1)
       q3lt=q(3,i )
       q3rt=q(3,i+1)
       rlt=q1lt
       ult=q2lt/rlt
       plt=(g-1.0)*(q3lt-0.5*rlt*ult**2)
       hlt=(q3lt+plt)/rlt
       rrt=q1rt
       urt=q2rt/q1rt
       prt=(g-1.0)*(q3rt-0.5*rrt*urt**2)
       hrt=(q3rt+prt)/rrt
       f1lt=rlt*ult
       f2lt=rlt*ult**2+plt
       f3lt=ult*(q3lt+plt)
       f1rt=rrt*urt
       f2rt=rrt*urt**2+prt
       f3rt=urt*(q3rt+prt)
       flux(1,i)=0.5*(f1lt+f1rt+tvd1)
       flux(2,i)=0.5*(f2lt+f2rt+tvd2)
       flux(3,i)=0.5*(f3lt+f3rt+tvd3)
    end do

    return
  end subroutine calflx

  subroutine calc_p12(g,tr,tl,pr,pl,p12)
    real,intent(in)::g,tr,tl,pr,pl
    real,intent(out)::p12
    real::f1,f2,p12_a,p12_b,p12_c

    p12_a=1.0
    p12_b=pl/pr
    eps=p12_b-p12_a
    do while(eps > 0.00001)
       p12_c=(p12_a+p12_b)/2
       f1=p12_a*((1-((g-1)*sqrt(tr/tl)*(p12_a-1))/(sqrt(2*g*(2*g+(g+1)*(p12_a-1))))))**(-2*g/(g-1))-pl/pr
       f2=p12_c*((1-((g-1)*sqrt(tr/tl)*(p12_c-1))/(sqrt(2*g*(2*g+(g+1)*(p12_c-1))))))**(-2*g/(g-1))-pl/pr
       if(f1*f2<0.0)then
          p12_b=p12_c
       else
          p12_a=p12_c
       end if
       eps=p12_b-p12_a
    end do
    p12=p12_c       
  end subroutine calc_p12
  

  subroutine exact(x,roue,ue,pe,ee,mx,roul,rour,pr,pl,p12,g,rgas,tr,tl,ur,ul,now)
    implicit none
    integer::i
    integer,intent(in)::mx
    real,intent(in)::x(-5:505)
    real,intent(out)::roue(-5:505),ue(-5:505),pe(-5:505),ee(-5:505)
    real::tmp,ptmp
    real,intent(in)::roul,rour,p12,g,rgas,tr,tl,ur,ul,now,pr,pl
    real::ms,us,al,a3,u2

    ms=sqrt((g+1)/(2*g)*(p12+(g-1)/(g+1)))
    us=ms*sqrt(g*rgas*tr)
    u2=2*sqrt(g*rgas*tr)/(g+1) * (ms-1/ms)
    al=sqrt(g*rgas*tr)
    a3=al-(g-1)/2*u2
    
    do i=1,mx
       if(x(i)>us*now)ue(i)=ur
       if(x(i)<=us*now .and. x(i)>(u2-a3)*now)ue(i)=u2
       if(x(i)<=(u2-a3)*now .and. x(i)>-al*now)then
          ue(i)=2*al/(g+1)*(1+x(i)/(al*now))
       end if
       if(x(i)<=-al*now)ue(i)=ul
    end do

    do i=1,mx
       if(x(i)<=-al*now)pe(i)=pl
       if(x(i)<=us*now .and. x(i)>-al*now)then
          tmp=(1-(g-1)/2*ue(i)/al)**(-2*g/(g-1))
          pe(i)=pl/tmp
       end if
       if(x(i)>us*now)pe(i)=pr
    end do

    do i=1,mx
       if(x(i)>us*now)then
          roue(i)=rour
       end if
       if(x(i)<=us*now .and. x(i)>u2*now)then
          tmp=pe(i)/pr
          roue(i)=((g+1)/(g-1)*tmp+1)/(tmp+(g+1)/(g-1))*rour
       end if
       if(x(i)<=-al*now)roue(i)=roul
       if(x(i)<=u2*now .and. x(i)>-al*now)then
          roue(i)=(pe(i)/pl)**(1/g)*roul
       end if
    end do

    do i=1,mx
       tmp=pe(i)
       pe(i)=tmp
       ee(i)=pe(i)/((g-1.0)*roue(i))
    end do
  
    return
  end subroutine exact
end program shocktube

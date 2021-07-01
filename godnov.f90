program godnov
  implicit none
  integer::mx,nlast,n
  integer::i
  real::x(-1:110)
  real::vecUr(2),vecUl(2),tmp1vecU(2),tmp2vecU(2)
  real::tmpf(2)
  real::ui(-1:110),vi(-1:110),ue(-1:110),ve(-1:110)
  real::us(-1:110),vs(-1:110)
  real::w1(-1:110),w2(-1:110)
  real::A(2,2),R(2,2),R_1(2,2),lamda(2,2)
  real::flux_u(-1:110),flux_v(-1:110)
  real::work_u(-1:110),work_v(-1:110)
  real::cfl,dx,dt,lamda1,lamda2,b,c,R_det
  real::wr(2),wl(2),tmp12U(2)
  real::RAR_1(2,2)
  
  
  A(1,1)=3.0
  A(2,1)=2.0
  A(1,2)=1.0
  A(2,2)=4.0

  vecUl(1)=1.0
  vecUl(2)=1.0
  vecUr(1)=2.0
  vecUr(2)=-1.0

  
  mx=61
  nlast=50
  dx=0.02
  dt=0.002
  
  !set grid
  do i=-1,mx+1
     x(i)=dx*float(i-11)
  end do

  !set initial condition
  do i=-1,mx+1
     if(x(i)>0.0)then
        ui(i)=vecUr(1)
        vi(i)=vecUr(2)
     else
        ui(i)=vecUl(1)
        vi(i)=vecUl(2)
     end if
  end do

  !solve eigenvalue
  b=-(A(1,1)+A(2,2))
  c=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  lamda1=(-b+sqrt(b*b-4.0*c))/2.0
  lamda2=(-b-sqrt(b*b-4.0*c))/2.0

  lamda(1,1)=abs(lamda1)
  lamda(2,1)=0.0
  lamda(1,2)=0.0
  lamda(2,2)=abs(lamda2)
  
  R(1,1)=A(1,2)
  R(2,1)=lamda1-A(1,1)
  R(1,2)=A(1,2)
  R(2,2)=lamda2-A(1,1)

  R_det=R(1,1)*R(2,2)-R(1,2)*R(2,1)

  R_1(1,1)=R(2,2)/R_det
  R_1(2,1)=-R(2,1)/R_det
  R_1(1,2)=-R(1,2)/R_det
  R_1(2,2)=R(1,1)/R_det

  !set exact solution
  wr=matmul(R_1,vecUr)
  wl=matmul(R_1,vecUl)

  do i=-1,mx+1
     if(x(i) > lamda1*dt*nlast)then
        w1(i)=wr(1)
     else
        w1(i)=wl(1)
     end if

     if(x(i)>lamda2*dt*nlast)then
        w2(i)=wr(2)
     else
        w2(i)=wl(2)
     end if
  end do
 
  do i=-1,mx+1
     ue(i)=R(1,1)*w1(i)+R(1,2)*w2(i)
     ve(i)=R(2,1)*w1(i)+R(2,2)*w2(i)
  end do

  !solve scheme
  RAR_1=matmul(R,matmul(lamda,R_1))

  do i=-1,mx+1
     us(i)=ui(i)
     vs(i)=vi(i)
  end do

  do n=1,nlast
     do i=0,mx-1
        tmp1vecU(1)=us(i)
        tmp1vecU(2)=vs(i)
        tmp2vecU(1)=us(i+1)
        tmp2vecU(2)=vs(i+1)
        tmp12U(1)=tmp2vecU(1)-tmp1vecU(1)
        tmp12U(2)=tmp2vecU(2)-tmp1vecU(2)
        tmpf=0.5*(matmul(A,tmp1vecU)+matmul(A,tmp2vecU)-matmul(RAR_1,tmp12U))
        flux_u(i)=tmpf(1)
        flux_v(i)=tmpf(2)
     end do

     do i=1,mx-1
        work_u(i)=us(i)-(dt/dx)*(flux_u(i)-flux_u(i-1))
        work_v(i)=vs(i)-(dt/dx)*(flux_v(i)-flux_v(i-1))
     end do

     do i=1,mx-1
        us(i)=work_u(i)
        vs(i)=work_v(i)
     end do
  end do

  !print
  do i=0,mx
     print*,x(i),ue(i),ve(i),us(i),vs(i)
  end do
  
end program godnov


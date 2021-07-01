program tvd
  implicit none
  integer::mx,nlast,n
  integer::i
  real::x(-1:100),yi(-1:100)
  real::ye(-1:100),ys3(-1:100)
  real::ys1(-1:100),ys2(-1:100)
  real::flux(-1:100)
  real::work(-1:100)
  real::c,ecp,cfl,dx,dt,travel
  real::dm,d0,dp,s,q,ac,f,sb1,sb2
  
  c=1.0
  ecp=0.1
  mx=51
  nlast=30
  cfl=0.5
  dx=0.02
  dt=cfl*dx
  travel=dt*c*float(nlast)

  !set grid
  do i=-1,mx+1
     x(i)=dx*float(i-1)
  end do

  !set initial condition
  do i=-1,mx+1
     if(x(i).lt.0.5)then
        yi(i)=1.0
     else
        yi(i)=0.0
     end if
  end do

  !set exact solution
  do i=-1,mx+1
     if(x(i).lt.0.5+travel)then
        ye(i)=1.0
     else
        ye(i)=0.0
     end if
  end do

  !solve scheme minmod1
  do i=-1,mx+1
     ys1(i)=yi(i)
  end do

  do n=1,nlast
     do i=0,mx-1
        dm=ys1(i)-ys1(i-1)
        d0=ys1(i+1)-ys1(i)
        dp=ys1(i+2)-ys1(i+1)
        s=sign(1.0,dm)
        q=s*amax1(0.0,amin1(s*dm,s*d0,s*dp))
        ac=abs(c)
        if(ac.lt.ecp) ac=(c*c+ecp*ecp)*0.5/ecp
        f=-(dt*c**2/dx*q+ac*(d0-q))
        flux(i)=0.5*(c*ys1(i)+c*ys1(i+1)+f)
     end do

     do i=1,mx-1
        work(i)=ys1(i)-(dt/dx)*(flux(i)-flux(i-1))
     end do

     do i=1,mx-1
        ys1(i)=work(i)
     end do
  end do

  !solve scheme by minmod2
  do i=-1,mx+1
     ys2(i)=yi(i)
  end do

  do n=1,nlast
     do i=0,mx-1
        dm=ys2(i)-ys2(i-1)
        d0=ys2(i+1)-ys2(i)
        dp=ys2(i+2)-ys2(i+1)
        s=sign(1.0,dm)
        q=s*amax1(0.0,amin1(2.0*s*dm,2.0*s*d0,2.0*s*dp,0.5*s*(dm+dp)))
        ac=abs(c)
        if(ac.lt.ecp) ac=(c*c+ecp*ecp)*0.5/ecp
        f=-(dt*c**2/dx*q+ac*(d0-q))
        flux(i)=0.5*(c*ys2(i)+c*ys2(i+1)+f)
     end do

     do i=1,mx-1
        work(i)=ys2(i)-cfl*(flux(i)-flux(i-1))
     end do

     do i=1,mx-1
        ys2(i)=work(i)
     end do
  end do

  !solve scheme by superbee
  do i=-1,mx+1
     ys3(i)=yi(i)
  end do

  do n=1,nlast
     do i=0,mx-1
        dm=ys3(i)-ys3(i-1)
        d0=ys3(i+1)-ys3(i)
        dp=ys3(i+2)-ys3(i+1)
        s=sign(1.0,dm)
        sb1=s*amax1(0.0,amin1(2.0*abs(d0),s*dp),amin1(abs(d0),2.0*s*d0))
        s=sign(1.0,d0)
        sb2=s*amax1(0.0,amin1(2.0*abs(d0),s*dp),amin1(abs(d0),2.0*s*d0))
        q=sb1+sb2-d0
        ac=abs(c)
        if(ac.lt.ecp) ac=(c*c+ecp*ecp)*0.5/ecp
        f=-(dt*c**2/dx*q+ac*(d0-q))
        flux(i)=0.5*(c*ys3(i)+c*ys3(i+1)+f)
     end do

     do i=1,mx-1
        work(i)=ys3(i)-cfl*(flux(i)-flux(i-1))
     end do

     do i=1,mx-1
        ys3(i)=work(i)
     end do
  end do

  do i=0,mx
     print*,x(i),yi(i),ye(i),ys1(i),ys2(i),ys3(i)
  end do
end program tvd

  

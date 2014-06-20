% sor.m -- Program to solve Laplace's equation in
% generalized coordinates using SOR.  Solve for
% potential flow over an airfoil.  Use an O-grid.
%
% Converted from sor.f
%
%    ni = number of points around the airfoil
%    nj = number of radial points
%
      parameter(ni=100,nj=75)

      real psi(ni,nj), x(ni,nj), y(ni,nj),...
        xxi(ni,nj), xeta(ni,nj), yxi(ni,nj), yeta(ni,nj), jinv(ni,nj),...
        a(ni,nj), b(ni,nj), c(ni,nj), d(ni,nj), e(ni,nj),...
        u(ni,nj), v(ni,nj), cp(ni,nj)

%...Read the grid and set up the metrics

      call read_grid(ni,nj,x,y)
      call get_metrics(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv)
      call get_coefs(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,a,b,c,d,e,cp)

%...Choose what type of initial conditions to use
      fprintf(' Start SOR iteration from:\n');
      fprintf('    [1] uniform flow initial conditions\n');
      fprintf('    [2] restart file (from a previous run)\n');
      fprintf(' Enter choice:\n');
      read(5,*)irestart

      if(irestart==1)
         fprintf(' Enter angle of attack (deg):');
         read(5,*)alpha
         %pi = acos(-1.0); uneeded since matlab has built in pi()
         alpha = alpha*pi()/180;
         call initial(ni,nj,x,y,psi,alpha)
         itno_rst = 0;
         open(unit=9,file='history.dat',status='unknown')
      elseif(irestart==2)
         call restart(ni,nj,itno_rst,psi,alpha)
         open(unit=9,file='history.dat',status='old',access='append')
      else
         fprintf('*** ERROR on input ***\n');
         stop
      end

%...Enter overrelaxation parameter and number of iterations

      fprintf(' Enter SOR overrelaxation parameter omega:\n');
      read(5,*)omega

      if((omega<=0.0)||(omega>=2.0))
         fprintf('*** ERROR omega out of range ***\n');
         stop
      elseif(omega<1.0)
         fprintf('*   WARNING:  using UNDER-relaxation *\n');
      end

      fprintf(' Enter number of SOR iterations:\n');
      read(5,*)itmax

%***
%***Start SOR iteration loop
%***

for itno = 1:itmax
    
%........Sweep over the interior points
    
for j = 2:nj-1
    for i = 2:ni-1
        
        psi(i,j) = (omega/(2.0*(a(i,j)+b(i,j))))*...
            ( 0.25*c(i,j)*psi(i-1,j-1)...
            + (b(i,j) - 0.5*e(i,j))*psi(i,j-1)...
            - 0.25*c(i,j)*psi(i+1,j-1)...
            + (a(i,j) - 0.5*d(i,j))*psi(i-1,j)...
            + (a(i,j) + 0.5*d(i,j))*psi(i+1,j)...
            - 0.25*c(i,j)*psi(i-1,j+1)...
            + (b(i,j) + 0.5*e(i,j))*psi(i,j+1)...
            + 0.25*c(i,j)*psi(i+1,j+1) )...
            + (1.0 - omega)*psi(i,j);
    end
end
    
%........Sweep over the periodic boundary
    
for j = 2:nj-1
    psi(1,j) = (omega/(2.0*(a(1,j)+b(1,j))))*...
        ( 0.25*c(1,j)*psi(ni-1,j-1)...
        + (b(1,j) - 0.5*e(1,j))*psi(1,j-1)...
        - 0.25*c(1,j)*psi(2,j-1)...
        + (a(1,j) - 0.5*d(1,j))*psi(ni-1,j)...
        + (a(1,j) + 0.5*d(1,j))*psi(2,j)...
        - 0.25*c(1,j)*psi(ni-1,j+1)...
        + (b(1,j) + 0.5*e(1,j))*psi(1,j+1)...
        + 0.25*c(1,j)*psi(2,j+1) )...
        + (1.0 - omega)*psi(1,j);
    psi(ni,j) = psi(1,j);
end

%........Set the BC on the airfoil surface and the outer boundary

for i = 1:ni
    %psi(i,1) = psi(1,2)
    psi(i,1) = (4.0*psi(1,2) - psi(1,3))/3.0;
    psi(i,nj) = cos(alpha)*y(i,nj) - sin(alpha)*x(i,nj);
end

%........Compute the velocity on the airfoil surface and the outer boundary

u(1,1) = -xeta(1,1)*jinv(1,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + xxi(1,1)*jinv(1,1)...
    *(-1.5*psi(1,1) + 2.0*psi(1,2) - 0.5*psi(1,3));
v(1,1) = -yeta(1,1)*jinv(1,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + yxi(1,1)*jinv(1,1)...
    *(-1.5*psi(1,1) + 2.0*psi(1,2) - 0.5*psi(1,3));
for i = 2:ni-1
    u(i,1) = -xeta(i,1)*jinv(i,1)*0.5*(psi(i+1,1) - psi(i-1,1))...
        + xxi(i,1)*jinv(i,1)...
        *(-1.5*psi(i,1) + 2.0*psi(i,2) - 0.5*psi(i,3));
    v(i,1) = -yeta(i,1)*jinv(i,1)*0.5*(psi(i+1,1) - psi(i-1,1))...
        + yxi(i,1)*jinv(i,1)...
        *(-1.5*psi(i,1) + 2.0*psi(i,2) - 0.5*psi(i,3));
end
u(ni,1) = -xeta(ni,1)*jinv(ni,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + xxi(ni,1)*jinv(ni,1)...
    *(-1.5*psi(ni,1) + 2.0*psi(ni,2) - 0.5*psi(ni,3));
v(ni,1) = -yeta(ni,1)*jinv(ni,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + yxi(ni,1)*jinv(ni,1)...
    *(-1.5*psi(ni,1) + 2.0*psi(ni,2) - 0.5*psi(ni,3));

%........Compute the pressure coefficient on the airfoil surface

for i = 1:ni
    cp(i,1) = 1.0 - (u(i,1)^2 + v(i,1)^2);
end

%........Integrate Cp to get the lift coefficient, cl

         cn = 0.0;
         ca = 0.0;
         cmle = 0.0;
         for i = 1:ni-1
            cn = cn - 0.5*(cp(i+1,1) + cp(i,1))*(x(i+1,1) - x(i,1));
            ca = ca + 0.5*(cp(i+1,1) + cp(i,1))*(y(i+1,1) - y(i,1));
            cmle = cmle + 0.5*(cp(i+1,1)*x(i+1,1) + cp(i,1)*x(i,1))...
              *(x(i+1,1) - x(i,1));
         end
         cl = cn*cos(alpha) - ca*sin(alpha);
         cd = cn*sin(alpha) + ca*cos(alpha);

         fprintf(' iteration:  ',itno_rst + itno,...
           ',  Cl = ',cl,',  Cd = ',cd,', CmLE = ',cmle); % don't know if this statement will work
         write(9,101)itno_rst+itno,cl,cd,cmle
101      format(3x,i7,3(3x,e15.7))

end
%***
%***End of SOR loop
%***

%...Compute velocity and Cp throughout the flow field

      call velocity(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,...
        psi,alpha,u,v);
    for i = 1:ni
        for j = 1:nj
            cp(i,j) = 1.0 - (u(i,j)^2 + v(i,j)^2);
        end
    end

%...Write out the solution, velocity and Cp

      call output(ni,nj,itno_rst+itmax,alpha,psi,u,v,cp)

      stop
      end

%********************************************************************
%   read_grid
%********************************************************************
      subroutine read_grid(ni,nj,x,y)

      real x(ni,nj), y(ni,nj)

      open(unit=7,file='100_75_x.dat')
      open(unit=8,file='100_75_y.dat')

      for j = 1:nj
         read(7,*)(x(i,j),i=1,ni)
         read(8,*)(y(i,j),i=1,ni)
      end

%...Check to make sure we are at the end of the file

      read(unit=7,end=101,fmt=*)dummy
      read(unit=8,end=101,fmt=*)dummy
      fprintf('*** Mismatch between ni,nj and size of grid ***\n');
      stop

101   continue

      close(7)
      close(8)

      return
      end

%********************************************************************
%   get_metrics
%********************************************************************
      subroutine get_metrics(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv)

      real x(ni,nj), y(ni,nj), xxi(ni,nj), xeta(ni,nj), yxi(ni,nj),...
        yeta(ni,nj), jinv(ni,nj);

%...Interior points

for j = 2:nj-1
    for i = 2:ni-1
        xxi(i,j) = 0.5*(x(i+1,j) - x(i-1,j));
        yxi(i,j) = 0.5*(y(i+1,j) - y(i-1,j));
        xeta(i,j) = 0.5*(x(i,j+1) - x(i,j-1));
        yeta(i,j) = 0.5*(y(i,j+1) - y(i,j-1));
    end
end

%...Boundary points:

%.....airfoil surface and outer boundary

for i = 2:ni-1
    xxi(i,1) = 0.5*(x(i+1,1) - x(i-1,1));
    yxi(i,1) = 0.5*(y(i+1,1) - y(i-1,1));
    xeta(i,1) = -1.5*x(i,1) + 2.0*x(i,2) - 0.5*x(i,3);
    yeta(i,1) = -1.5*y(i,1) + 2.0*y(i,2) - 0.5*y(i,3);
end
for i = 2:ni-1
    xxi(i,nj) = 0.5*(x(i+1,nj) - x(i-1,nj));
    yxi(i,nj) = 0.5*(y(i+1,nj) - y(i-1,nj));
    xeta(i,nj) = 1.5*x(i,nj) - 2.0*x(i,nj-1) + 0.5*x(i,nj-2);
    yeta(i,nj) = 1.5*y(i,nj) - 2.0*y(i,nj-1) + 0.5*y(i,nj-2);
end

%.....periodic boundaries

for j = 2:nj-1
    xxi(1,j) = 0.5*(x(2,j) - x(ni-1,j));
    yxi(1,j) = 0.5*(y(2,j) - y(ni-1,j));
    xeta(1,j) = 0.5*(x(1,j+1) - x(1,j-1));
    yeta(1,j) = 0.5*(y(1,j+1) - y(1,j-1));
end
for j = 2:nj-1
    xxi(ni,j) = 0.5*(x(2,j) - x(ni-1,j));
    yxi(ni,j) = 0.5*(y(2,j) - y(ni-1,j));
    xeta(ni,j) = 0.5*(x(ni,j+1) - x(ni,j-1));
    yeta(ni,j) = 0.5*(y(ni,j+1) - y(ni,j-1));
end

%.....corners

      xxi(1,1) = 0.5*(x(2,1) - x(ni-1,1));
      yxi(1,1) = 0.5*(y(2,1) - y(ni-1,1));
      xeta(1,1) = -1.5*x(1,1) + 2.0*x(1,2) - 0.5*x(1,3);
      yeta(1,1) = -1.5*y(1,1) + 2.0*y(1,2) - 0.5*y(1,3);
      xxi(ni,1) = xxi(1,1);
      yxi(ni,1) = yxi(1,1);
      xeta(ni,1) = xeta(1,1);
      yeta(ni,1) = yeta(1,1);

      xxi(1,nj) = 0.5*(x(2,nj) - x(ni-1,nj));
      yxi(1,nj) = 0.5*(y(2,nj) - y(ni-1,nj));
      xeta(1,nj) = 1.5*x(1,nj) - 2.0*x(1,nj-1) + 0.5*x(1,nj-2);
      yeta(1,nj) = 1.5*y(1,nj) - 2.0*y(1,nj-1) + 0.5*y(1,nj-2);
      xxi(ni,nj) = xxi(1,nj);
      yxi(ni,nj) = yxi(1,nj);
      xeta(ni,nj) = xeta(1,nj);
      yeta(ni,nj) = yeta(1,nj);

%...Compute 1/J where J is the determinant of the Jacobian

for j = 1:nj
    for i = 1:ni
        jinv(i,j) = 1./(xxi(i,j)*yeta(i,j) - xeta(i,j)*yxi(i,j));
    end
end

      return
      end

%********************************************************************
%   get_coefs
%********************************************************************

      subroutine get_coefs(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,...
        a,b,c,d,e,tmp)

      real x(ni,nj), y(ni,nj), xxi(ni,nj), xeta(ni,nj), yxi(ni,nj),...
        yeta(ni,nj), jinv(ni,nj), a(ni,nj), b(ni,nj), c(ni,nj),...
        d(ni,nj), e(ni,nj), tmp(ni,nj);

%...Temporarilly store 2nd derivatives of x and y in a-tmp

      call get_2nd_derivs(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,...
        a,b,c,d,e,tmp);

%...Compute coefficients in the Laplacian

for j = 1:nj
    for i = 1:ni
        r1 = -(yeta(i,j)*jinv(i,j))^2*a(i,j)...
            + 2.0*yxi(i,j)*yeta(i,j)*jinv(i,j)^2*b(i,j)...
            - (yxi(i,j)*jinv(i,j))^2*c(i,j);
        r2 = -(yeta(i,j)*jinv(i,j))^2*d(i,j)...
            + 2.0*yxi(i,j)*yeta(i,j)*jinv(i,j)^2*e(i,j)...
            - (yxi(i,j)*jinv(i,j))^2*tmp(i,j);
        d1 = yeta(i,j)*jinv(i,j)*r1 - xeta(i,j)*jinv(i,j)*r2;
        e1 = -yxi(i,j)*jinv(i,j)*r1 + xxi(i,j)*jinv(i,j)*r2;
        
        s1 = -(xeta(i,j)*jinv(i,j))^2*a(i,j)...
            + 2.0*xxi(i,j)*xeta(i,j)*jinv(i,j)^2*b(i,j)...
            - (xxi(i,j)*jinv(i,j))^2*c(i,j);
        s2 = -(xeta(i,j)*jinv(i,j))^2*d(i,j)...
            + 2.0*xxi(i,j)*xeta(i,j)*jinv(i,j)^2*e(i,j)...
            - (xxi(i,j)*jinv(i,j))^2*tmp(i,j);
        d2 = yeta(i,j)*jinv(i,j)*s1 - xeta(i,j)*jinv(i,j)*s2;
        e2 = -yxi(i,j)*jinv(i,j)*s1 + xxi(i,j)*jinv(i,j)*s2;
        
        e(i,j) = e1 + e2;
        d(i,j) = d1 + d2;
        c(i,j) = -2.0*(yxi(i,j)*yeta(i,j) + xxi(i,j)*xeta(i,j))...
            *jinv(i,j)^2;
        b(i,j) = (yxi(i,j)*jinv(i,j))^2 + (xxi(i,j)*jinv(i,j))^2;
        a(i,j) = (yeta(i,j)*jinv(i,j))^2 + (xeta(i,j)*jinv(i,j))^2;
    end
end

      return
      end

%********************************************************************
%   get_2nd_derivs
%********************************************************************

      subroutine get_2nd_derivs(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,...
        a,b,c,d,e,tmp);

      real x(ni,nj), y(ni,nj), xxi(ni,nj), xeta(ni,nj), yxi(ni,nj),...
        yeta(ni,nj), jinv(ni,nj), a(ni,nj), b(ni,nj), c(ni,nj),...
        d(ni,nj), e(ni,nj), tmp(ni,nj);

%...Temporarilly use arrays to store needed second derivatives of x and y
%      a = xxixi
%      b = xxieta
%      c = xetaeta
%      d = yxixi
%      e = yxieta
%      tmp = yetaeta
%

%...interior points

for j = 2:nj-1
    for i = 2:ni-1
        a(i,j) = x(i+1,j) -2.0*x(i,j) + x(i-1,j);
        b(i,j) = 0.25*(x(i+1,j+1) - x(i-1,j+1)...
            - x(i+1,j-1) + x(i-1,j-1));
        c(i,j) = x(i,j+1) -2.0*x(i,j) + x(i,j-1);
        d(i,j) = y(i+1,j) -2.0*y(i,j) + y(i-1,j);
        e(i,j) = 0.25*(y(i+1,j+1) - y(i-1,j+1)...
            - y(i+1,j-1) + y(i-1,j-1));
        tmp(i,j) = y(i,j+1) -2.0*y(i,j) + y(i,j-1);
    end
end

%...boundary points

%.....airfoil surface

for i = 2:ni-1
    a(i,1) = x(i+1,1) -2.0*x(i,1) + x(i-1,1);
    b(i,1) = - 0.75*(x(i+1,1) - x(i-1,1))...
        + (x(i+1,2) - x(i-1,2))...
        - 0.25*(x(i+1,3) - x(i-1,3));
    c(i,1) = x(i,1) -2.0*x(i,2) + x(i,3);
    d(i,1) = y(i+1,1) -2.0*y(i,1) + y(i-1,1);
    e(i,1) = - 0.75*(y(i+1,1) - y(i-1,1))...
        + (y(i+1,2) - y(i-1,2))...
        - 0.25*(y(i+1,3) - y(i-1,3));
    tmp(i,1) = y(i,1) -2.0*y(i,2) + y(i,3);
end

%.....outer boundary

for i = 2:ni-1
    a(i,nj) = x(i+1,nj) -2.0*x(i,nj) + x(i-1,nj);
    b(i,nj) = 0.75*(x(i+1,nj) - x(i-1,nj))...
        - (x(i+1,nj-1) - x(i-1,nj-1))...
        + 0.25*(x(i+1,nj-2) - x(i-1,nj-2));
    c(i,nj) = x(i,nj) -2.0*x(i,nj-1) + x(i,nj-2);
    d(i,nj) = y(i+1,nj) -2.0*y(i,nj) + y(i-1,nj);
    e(i,nj) = 0.75*(y(i+1,nj) - y(i-1,nj))...
        - (y(i+1,nj-1) - y(i-1,nj-1))...
        + 0.25*(y(i+1,nj-2) - y(i-1,nj-2));
    tmp(i,nj) = y(i,nj) -2.0*y(i,nj-1) + y(i,nj-2);
end

%.....periodic boundaries

for j = 2:nj-1
    a(1,j) = x(2,j) -2.0*x(1,j) + x(ni-1,j);
    b(1,j) = 0.25*(x(2,j+1) - x(ni-1,j+1)...
        - x(2,j-1) + x(ni-1,j-1));
    c(1,j) = x(1,j+1) -2.0*x(1,j) + x(1,j-1);
    d(1,j) = y(2,j) -2.0*y(1,j) + y(ni-1,j);
    e(1,j) = 0.25*(y(2,j+1) - y(ni-1,j+1)...
        - y(2,j-1) + y(ni-1,j-1));
    tmp(1,j) = y(1,j+1) -2.0*y(1,j) + y(1,j-1);
end

for j = 2:nj-1
    a(ni,j) = a(1,j);
    b(ni,j) = b(1,j);
    c(ni,j) = c(1,j);
    d(ni,j) = d(1,j);
    e(ni,j) = e(1,j);
    tmp(ni,j) = tmp(1,j);
end

%.....corners

a(1,1) = a(2,1);
b(1,1) = - 0.75*(x(2,1) - x(ni-1,1))...
    + (x(2,2) - x(ni-1,2))...
    - 0.25*(x(2,3) - x(ni-1,3));
c(1,1) = c(1,2);
d(1,1) = d(2,1);
e(1,1) = - 0.75*(y(2,1) - y(ni-1,1))...
    + (y(2,2) - y(ni-1,2))...
    - 0.25*(y(2,3) - y(ni-1,3));
tmp(1,1) = tmp(1,2);

a(1,nj) = a(2,nj);
b(1,nj) = 0.75*(x(2,nj) - x(ni-1,nj))...
    - (x(2,nj-1) - x(ni-1,nj-1))...
    + 0.25*(x(2,nj-2) - x(ni-1,nj-2));
c(1,nj) = c(1,nj-1);
d(1,nj) = d(2,nj);
e(1,nj) = 0.75*(y(2,nj) - y(ni-1,nj))...
    - (y(2,nj-1) - y(ni-1,nj-1))...
    + 0.25*(y(2,nj-2) - y(ni-1,nj-2));
tmp(1,nj) = tmp(1,nj-1);

a(ni,nj) = a(1,nj);
b(ni,nj) = b(1,nj);
c(ni,nj) = c(1,nj);
d(ni,nj) = d(1,nj);
e(ni,nj) = e(1,nj);
tmp(ni,nj) = tmp(1,nj);

a(ni,1) = a(1,1);
b(ni,1) = b(1,1);
c(ni,1) = c(1,1);
d(ni,1) = d(1,1);
e(ni,1) = e(1,1);
tmp(ni,1) = tmp(1,1);

      return
      end

%********************************************************************
%   initial - Set the initial conditions to be uniform flow
%   at an angle alpha
%********************************************************************

      subroutine initial(ni,nj,x,y,psi,alpha)

      real x(ni,nj), y(ni,nj), psi(ni,nj)

      for j = 1:nj
          for i = 1:ni
              psi(i,j) = cos(alpha)*y(i,j) - sin(alpha)*x(i,j)
              %psi(i,j) = 0.0
          end
      end

      return
      end

%********************************************************************
%   restart - start the iteration from a restart file
%********************************************************************

      subroutine restart(ni,nj,itno_rst,psi,alpha)

      real psi(ni,nj)

      open(unit=10,file='restart.dat',status='old',...
        form='unformatted');

      read(10)itno_rst,ni_rst,nj_rst,alpha
      for j = 1:nj
         read(10)(psi(i,j),i=1,ni)
      end

%...Check to make sure we are at the end of the file

      read(unit=10,end=101)dummy
      fprintf(...
     '*** Mismatch between ni,nj and size of restart file ***\n');
      stop

101   continue
      if((ni_rst~=ni)||(nj_rst~=nj))
         fprintf(...
        '*** Error:  Mismatch in ni, nj and ni_rst, nj_rst ***\n');
         stop
      end

      close(10)

      return
      end

%********************************************************************
%   velocity - compute the velocity throughout the flow field
%********************************************************************
      subroutine velocity(ni,nj,x,y,xxi,xeta,yxi,yeta,jinv,...
        psi,alpha,u,v);

      real x(ni,nj), y(ni,nj),...
        xxi(ni,nj), xeta(ni,nj), yxi(ni,nj), yeta(ni,nj), jinv(ni,nj),...
        psi(ni,nj), u(ni,nj), v(ni,nj);


%.....Compute the velocity on the airfoil surface and the outer boundary

u(1,1) = -xeta(1,1)*jinv(1,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + xxi(1,1)*jinv(1,1)...
    *(-1.5*psi(1,1) + 2.0*psi(1,2) - 0.5*psi(1,3));
v(1,1) = -yeta(1,1)*jinv(1,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + yxi(1,1)*jinv(1,1)...
    *(-1.5*psi(1,1) + 2.0*psi(1,2) - 0.5*psi(1,3));
for i = 2:ni-1
    u(i,1) = -xeta(i,1)*jinv(i,1)*0.5*(psi(i+1,1) - psi(i-1,1))...
        + xxi(i,1)*jinv(i,1)...
        *(-1.5*psi(i,1) + 2.0*psi(i,2) - 0.5*psi(i,3));
    v(i,1) = -yeta(i,1)*jinv(i,1)*0.5*(psi(i+1,1) - psi(i-1,1))...
        + yxi(i,1)*jinv(i,1)...
        *(-1.5*psi(i,1) + 2.0*psi(i,2) - 0.5*psi(i,3));
end
u(ni,1) = -xeta(ni,1)*jinv(ni,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + xxi(ni,1)*jinv(ni,1)...
    *(-1.5*psi(ni,1) + 2.0*psi(ni,2) - 0.5*psi(ni,3));
v(ni,1) = -yeta(ni,1)*jinv(ni,1)*0.5*(psi(2,1) - psi(ni-1,1))...
    + yxi(ni,1)*jinv(ni,1)...
    *(-1.5*psi(ni,1) + 2.0*psi(ni,2) - 0.5*psi(ni,3));

for i = 1:ni
    u(i,nj) = cos(alpha);
    v(i,nj) = sin(alpha);
end

%.....interior points

for j = 2:nj-1
    for i = 2:ni-1
        u(i,j) = -xeta(i,j)*jinv(i,j)*0.5*(psi(i+1,j) - psi(i-1,j))...
            + xxi(i,j)*jinv(i,j)*0.5*(psi(i,j+1) - psi(i,j-1));
        v(i,j) = -yeta(i,j)*jinv(i,j)*0.5*(psi(i+1,j) - psi(i-1,j))...
            + yxi(i,j)*jinv(i,j)*0.5*(psi(i,j+1) - psi(i,j-1));
    end
end

%.....periodic boundary

for j = 2:nj-1
    u(1,j) = -xeta(1,j)*jinv(1,j)*0.5*(psi(2,j) - psi(ni-1,j))...
        + xxi(1,j)*jinv(1,j)*0.5*(psi(1,j+1) - psi(1,j-1));
    v(1,j) = -yeta(1,j)*jinv(1,j)*0.5*(psi(2,j) - psi(ni-1,j))...
        + yxi(1,j)*jinv(1,j)*0.5*(psi(1,j+1) - psi(1,j-1));
    u(ni,j) = u(1,j);
    v(ni,j) = v(1,j);
end

      return
      end

%********************************************************************
%   output - write the output files for plotting and the restart file
%********************************************************************
      subroutine output(ni,nj,itno_rst,alpha,psi,u,v,cp)

      real psi(ni,nj), u(ni,nj), v(ni,nj), cp(ni,nj)
      character*80 myformat

%...Write an unformatted restart file

      open(unit=10,file='restart.dat',status='unknown',...
        form='unformatted')

      write(10)itno_rst,ni,nj,alpha
      for j = 1:nj
         write(10)(psi(i,j),i=1,ni)
      end

      close(10)

%...Write output files for plotting

      open(unit=11,file='psi.dat')
      open(unit=12,file='u.dat')
      open(unit=13,file='v.dat')
      open(unit=14,file='cp.dat')

      if((ni<10)&&(ni>2))
         write(myformat,101)ni
101      format('(',i1,'(3x,e15.7))')
      elseif(ni<100)
         write(myformat,102)ni
102      format('(',i2,'(3x,e15.7))')
      elseif(ni<1000)
         write(myformat,103)ni
103      format('(',i3,'(3x,e15.7))')
      elseif(ni<10000)
         write(myformat,104)ni
104      format('(',i4,'(3x,e15.7))')
      else
         fprintf('*** Error:  ni too large or too small ***\n');
      end

      for j = 1:nj
         write(11,myformat)(psi(i,j),i=1,ni)
         write(12,myformat)(u(i,j),i=1,ni)
         write(13,myformat)(v(i,j),i=1,ni)
         write(14,myformat)(cp(i,j),i=1,ni)
      end

      return
      end

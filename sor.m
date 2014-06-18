% sor.mor -- Program to solve Laplace's equation in
% generalized coordinates using SOR.  Solve for
% potential flow over an airfoil.  Use an O-grid.
%    ni = number of points around the airfoil
%    nj = number of radial points
%
ni=100;
nj=75;



%********************************************************************
%   read_grid
%********************************************************************
function []=read_grid(ni,nj,x,y)

      real x(ni,nj), y(ni,nj)

      open(unit=7,file='100_75_x.dat')
      open(unit=8,file='100_75_y.dat')

      for j = 1:nj
         read(7,*)(x(i,j),i=1,ni)
         read(8,*)(y(i,j),i=1,ni)
      end
end
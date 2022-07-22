function dy=double_well1D(t,x,beta)
%global N r1 r2 r3 D alpha
global c r1 r2 r3 D
% r1=1;
% r2=3;
% r3=5;
% D=0.01;
alpha=1;

dy=-(x-r1).*(x-r2).*(x-r3)+D*beta*x.^alpha +c;
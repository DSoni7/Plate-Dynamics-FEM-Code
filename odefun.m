function yd = odefun(t,y,tspan)
dof=1224;
k=readmatrix('K.txt');
m=readmatrix('M.txt');
f1=readmatrix('F.txt');
c=zeros(dof);

%yd=(K-M*omega^2);
%yd=[y(2); -(c*m)*y(2)-(k*m)*y(1)];
z=zeros(dof);
e=eye(dof);
A=[z, e; -m*k, -m*c];
F=[zeros(dof,1);(m)*f1];
Svec=y(1:dof*2);
yd=A*Svec+F;
end
%% Eigenfrequency Analysis
clc
clear all
%Finite Element Formulation
nx=24;
ny=16;
nn=(nx+1)*(ny+1); %number of nodes
nd=3*nn;          %number of dofs 

%Geometry properties
L=0.4;
W=0.6;
a=L/(2*nx);    %length of each element
b=W/(2*ny);    %width of each element
h=0.00625;

%Structure Material Properties
Y=70e9;
v=0.3;
rho=2700;
Edash=Y/(1-v*v);
G=Y/(2*(1+v));
Dmat=[Edash,   Edash*v,    0;
      Edash*v, Edash,      0;
      0,       0,          G];

%FLUID MATERIAL PROPERTIES
rhof=1000;
C=-1;
muf=pi*(((L^-2)+(W^-2))^0.5);
h1=1;  %hieght of fluid above plate
h2=1;  %hieght of fluid below plate
Zf1=(-rhof/muf)*(1+C*exp(2*muf*h1))/(1-C*exp(2*muf*h1));
Zf2=(-rhof/muf)*(exp(muf*(-2*h2+h))+1)/(exp(muf*(-2*h2+h))-1);
Zf3=Zf1+Zf2;

%LOCAL MATRICES
x=sym('x');
y=sym('y');
aa=1;bb=1;
n1=@(xx,yy)(1/8*((1+x*xx)*(1+y*yy)*(2+x*xx+y*yy-x^2-y^2)));
n2=@(xx,yy)(b/8*(1+x*xx)*(y+yy)*(y^2-1));
n3=@(xx,yy)((-a/8)*(x+xx)*(x^2-1)*(1+y*yy));
N=[];nodes=[-1,-1; 1 -1; 1 1; -1 1];
for i=1:4
    N=[N;n1(nodes(i,1),nodes(i,2));n2(nodes(i,1),nodes(i,2));n3(nodes(i,1),nodes(i,2))];
end

P1=(1/(a^2)).*diff(N,x,2);
P2=(1/(b^2)).*diff(N,y,2); mx3=diff(N,x);
P3=(2/(a*b)).*diff(mx3,y);
P=[P1'; P2'; P3'];
B=eval(int(int(P'*Dmat*P,x,-1,1),y,-1,1));

%ELEMENT MATRICES
mstruct=(rho*h*a*b).*(eval(int(int(N*N',x,-1,1),y,-1,1)));
klocal=((h^3)*a*b/12).*B;
p=10;  %pressure
fe=(a*b*p).*(eval(int(int(N',x,-1,1),y,-1,1)))';
mfluid=(Zf3*a*b).*(eval(int(int(N*N',x,-1,1),y,-1,1)));
mlocal=mstruct+mfluid;
%Connectivity Matrix
for i=1:(nx)*(ny)
        R=fix((i-1)/nx+1);
        C=mod((i-1),nx)+1;
        LG(i,1)=(R-1)*(nx+1)+C;
        LG(i,2)=(R-1)*(nx+1)+C+1;
        LG(i,3)=R*(nx+1)+C+1;
        LG(i,4)=R*(nx+1)+C;
end
j=1;dof_node=[];
for i=1:nn
    dof_node=[dof_node;j j+1 j+2];
    j=j+3;
end

%GLOBAL MATRIX FORMULATION
KWW=zeros(nn,nn);  KWX=zeros(nn,nn);  KWY=zeros(nn,nn);  FX=zeros(nn,1);
KXW=zeros(nn,nn);  KXX=zeros(nn,nn);  KXY=zeros(nn,nn);  FY=zeros(nn,1);
KYW=zeros(nn,nn);  KYX=zeros(nn,nn);  KYY=zeros(nn,nn);  FZ=zeros(nn,1);
MWW=zeros(nn,nn);  MWX=zeros(nn,nn);  MWY=zeros(nn,nn);
MXW=zeros(nn,nn);  MXX=zeros(nn,nn);  MXY=zeros(nn,nn);
MYW=zeros(nn,nn);  MYX=zeros(nn,nn);  MYY=zeros(nn,nn);
for i=1:4
    for j=1:4
        kww(i,j)=klocal(3*i-2,3*j-2);mww(i,j)=mlocal(3*i-2,3*j-2); fz(i,1)=fe(3*i-2,1);
        kwx(i,j)=klocal(3*i-2,3*j-1);mwx(i,j)=mlocal(3*i-2,3*j-1); fx(i,1)=fe(3*i-1,1);
        kwy(i,j)=klocal(3*i-2,3*j);  mwy(i,j)=mlocal(3*i-2,3*j);   fy(i,1)=fe(3*i,1);
        kxw(i,j)=klocal(3*i-1,3*j-2);mxw(i,j)=mlocal(3*i-1,3*j-2);
        kxx(i,j)=klocal(3*i-1,3*j-1);mxx(i,j)=mlocal(3*i-1,3*j-1);
        kxy(i,j)=klocal(3*i-1,3*j);  mxy(i,j)=mlocal(3*i-1,3*j);
        kyw(i,j)=klocal(3*i,3*j-2);  myw(i,j)=mlocal(3*i,3*j-2);
        kyx(i,j)=klocal(3*i,3*j-1);  myx(i,j)=mlocal(3*i,3*j-1);
        kyy(i,j)=klocal(3*i,3*j);    myy(i,j)=mlocal(3*i,3*j);
    end
end
for i=1:4
    for j=1:4
        for e=1:(nx*ny)
            I=LG(e,i);
            J=LG(e,j);
            KWW(I,J)=KWW(I,J)+kww(i,j);  MWW(I,J)=MWW(I,J)+mww(i,j);  
            KWX(I,J)=KWX(I,J)+kwx(i,j);  MWX(I,J)=MWX(I,J)+mwx(i,j);      
            KWY(I,J)=KWY(I,J)+kwy(i,j);  MWY(I,J)=MWY(I,J)+mwy(i,j);  
            KXW(I,J)=KXW(I,J)+kxw(i,j);  MXW(I,J)=MXW(I,J)+mxw(i,j);
            KXX(I,J)=KXX(I,J)+kxx(i,j);  MXX(I,J)=MXX(I,J)+mxx(i,j);       
            KXY(I,J)=KXY(I,J)+kxy(i,j);  MXY(I,J)=MXY(I,J)+mxy(i,j);
            KYW(I,J)=KYW(I,J)+kyw(i,j);  MYW(I,J)=MYW(I,J)+myw(i,j);
            KYX(I,J)=KYX(I,J)+kyx(i,j);  MYX(I,J)=MYX(I,J)+myx(i,j);       
            KYY(I,J)=KYY(I,J)+kyy(i,j);  MYY(I,J)=MYY(I,J)+myy(i,j);
        end
    end
end
for i=1:4
    for e=1:(nx*ny)
        I=LG(e,i);
        FX(I)=FX(I)+fx(i);
        FY(I)=FY(I)+fy(i); 
        FZ(I)=FZ(I)+fz(i); 
    end
end
Kglobal=zeros(nd,nd); Mglobal=zeros(nd,nd); Fglobal=zeros(nd,1);

for i=1:nn
    for j=1:nn
        Kglobal(3*i-2,3*j-2)=KWW(i,j); Mglobal(3*i-2,3*j-2)=MWW(i,j);  
        Kglobal(3*i-2,3*j-1)=KWX(i,j); Mglobal(3*i-2,3*j-1)=MWX(i,j);  
        Kglobal(3*i-2,3*j)=KWY(i,j);   Mglobal(3*i-2,3*j)=MWY(i,j);    
        Kglobal(3*i-1,3*j-2)=KXW(i,j); Mglobal(3*i-1,3*j-2)=MXW(i,j);
        Kglobal(3*i-1,3*j-1)=KXX(i,j); Mglobal(3*i-1,3*j-1)=MXX(i,j);
        Kglobal(3*i-1,3*j)=KXY(i,j);   Mglobal(3*i-1,3*j)=MXY(i,j);
        Kglobal(3*i,3*j-2)=KYW(i,j);   Mglobal(3*i,3*j-2)=MYW(i,j);
        Kglobal(3*i,3*j-1)=KYX(i,j);   Mglobal(3*i,3*j-1)=MYX(i,j);
        Kglobal(3*i,3*j)=KYY(i,j);     Mglobal(3*i,3*j)=MYY(i,j);
    end
end
for i=1:nn
    Fglobal(3*i-2)=FZ(i);
    Fglobal(3*i-1)=FX(i);
    Fglobal(3*i)=FY(i);
end
K=Kglobal;
M=Mglobal;
Force=Fglobal;
%NODES MATRIX ALOND THE BOUNDARY
for i=1:nx+1
    LX(1,i)=i; %y=0 axis
    LX(2,i)=i+ny*(nx+1); %y=b axis
end
for j=1:ny+1
    LY(1,j)=(j-1)*(nx+1)+1; %x=0 axis
    LY(2,j)=j*(nx+1);    %x=a axis
end
%BOUNDARY CONDITION for clamped/simply supported
%rows2remove=[dof_node(LX(1,:),1:3)' dof_node(LX(2,:),1:3)' dof_node(LY(1,:),1:3)' dof_node(LY(2,:),1:3)'];
%cols2remove=[dof_node(LX(1,:),1:3)' dof_node(LX(2,:),1:3)' dof_node(LY(1,:),1:3)' dof_node(LY(2,:),1:3)'];
%BC for Cantilever
rows2remove=[dof_node(LX(1,:),1)' dof_node(LX(1,:),2)' dof_node(LX(1,:),3)'];
cols2remove=[dof_node(LX(1,:),1)' dof_node(LX(1,:),2)' dof_node(LX(1,:),3)'];
K(rows2remove,:)=[]; K(:,cols2remove)=[];
M(rows2remove,:)=[]; M(:,cols2remove)=[];

%Eigenfrequency
[V,D]=eig(K,M); 
omega=sqrt(D);     
freq=(omega/(2*pi));  
Freq=sort(diag(freq));
%PLATE
%var=(nx/2)*(ny-1); %for plate center
var=(nx/2)*(2*ny-1); %for plate end center
%disp=dof_node(LG(var,3),1)-nx-ny; %Simply supported all sides
%disp=dof_node(LG(var,3),1)-3*(nx+1)-6*((ny/2)-1)-3; %clamped all sides
disp=dof_node(LG(var,3),1)-3*(nx+1); %CANTILEVER
dof=size(K,1);
%% MODE SHAPES
% Vibration modes (first, second and third modes)
clear X Z
clc
for i=1:(size(V,1)/3)
    X(i,1)=V(3*i-2,6);
end
R=reshape(X,[nx+1,ny]);
Z=zeros(nx+1,ny+1);
Z(:,2:ny+1)=R;
surf(Z)
%% DYNAMIC RESPONSE
clear z3
%PVC
%k1=2.983e7; d1=79227; k2=1e5*(5.68e7);
%y1=k1; y2=d1+(k1*d1/k2); y4=d1/k2;
%Aluminium
k1=7.13e10; d1=4.79e8; k2=2.3162e10;
y1=k1; y2=d1+(k1*d1/k2); y4=d1/k2;
%a0=7.13e10; b0=1; a1=4.79e5; b1=2.068e-9;
%y1=a0; y2=a1; y4=b1;
i=inv(y4.*M);
j=M;
k=y2.*K;
l=y1.*K;
%ODE METHOD
idx=0;
dt=0.01; 
T=1;
s=[0:dt:T];
nt=(T/dt); %number of time step
%Force on plate
F=0*Force*sin(4*s); Fdot=0*y4*4*Force*cos(4*s);
F(rows2remove,:)=[]; Fdot(rows2remove,:)=[];
y0=ones(dof,1); yd0=zeros(dof,1); ydd0=zeros(dof,1);
%Order 2
%c=(d1).*K; k=(k1).*K; m=inv(M); y_int=[y0;yd0];
%[t1,z2]=ode15s(@(t,z)odefun_order2(t,z,m,k,F,idx,c),s,y_int);
%Order 3
y_int=[y0;yd0;ydd0];
[t1,z3]=ode45(@(t,z)odefun_order3(t,z,F,idx,i,j,k,l,Fdot),s,y_int);
plot(z3(:,disp),'k')
clear F Fdot
%% Newmark method
clear X
%PVC
%k1=2.983e7; d1=0*k1;
%Aluminium
k1=7.13e10; d1=0e-2*k1;%4.79e8;
dt=0.01; 
T=1;
s=[0:dt:T];
nt=(T/dt); %number of time step
%Force on plate
for i=1:nt+1
    F(:,i)=Force;
end
F(rows2remove,:)=[];
%with damping
c=(d1).*K;
k=(k1).*K; m=inv(M);
B=1/4;  %beta=1/4
G=1/2;  %gamma=1/2
X(:,1)=zeros(dof,1);
xd(:,1)=zeros(dof,1);
xdd(:,1)=m*(F(:,1)-c*xd(:,1)-k*X(:,1));
c1=(1/(B*dt*dt)).*M+(G/(B*dt)).*c +k;
ic1=inv(c1);
for i=1:nt
    c2=M*((X(:,i)./(B*dt*dt))+(xd(:,i)./(B*dt))+((1/(2*B))-1).*xdd(:,i));
    c3=(G/(B*dt)).*X(:,i)+((G/B)-1).*xd(:,i)+(((G/B)-2)*dt/2).*xdd(:,i);
    X(:,i+1)=ic1*(F(:,i+1)+c2+c*c3);
    xdd(:,i+1)=(1/(B*dt*dt)).*(X(:,i+1)-X(:,i))-(1/(B*dt)).*xd(:,i)-((1/(2*B))-1).*xdd(:,i);
    xd(:,i+1)=xd(:,i)+(1-G)*dt.*xdd(:,i)+G*dt.*xdd(:,i+1);
end
hold on
plot(X(disp,:),'b')
hold on
%legend('Viscoelastic plate','Elastic plate')
clear F 
%% PLOTS
al1000=readmatrix('al1000.txt');
plot(al1000,'r')
hold on
al100=readmatrix('al100.txt');
plot(al100,'k')
hold on
al10=readmatrix('al10.txt');
plot(al10,'b')
legend('1000k2', "100k2",'10k2','2-element model')
function [xHatUpdated r S Ppost]=EKF(t,y,u,xinit,dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% EXTENDED KF ROUTINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Tuning Parameters  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cbn=20000; %current bias noise (Force)
Wdn=.8; %wave disturbance noise (Position)
Pn=.01; %Position noise
Vn=.0001; %Velocity noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Initialize memory  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent P xHat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Initialize filter  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
if isempty (P)  
%%% Wave Parameters%%% 
    P=zeros(15); 
    xHat=xinit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Initialize system  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M= 10^9*[0.007 0 0;
         0 0.0110 -0.0130;
         0 -0.0130 3.1930];

Minv=inv(M);

D=[200000 0 0;
    0 100000 -700000;
    0 -700000 63900000];


    
T1=0.0001;T2=0.0001;T3=0.0001;
T=diag([T1 T2 T3]);

w01=0.9; w02=0.9; w03=0.9;
z1=0.1; z2=0.1; z3=0.1;
sig1=50; sig2=50; sig3=50;


Omega=diag([w01 w02 w03]);
Z=diag([z1 z2 z3]);
Sigma=diag([sig1 sig2 sig3]);

Sigma2=[zeros(3,3); 
        Sigma];
Omega2=[zeros(3,3) eye(3);
        -Omega^2 -2*Z*Omega];
Psi=eye(3);

psi0=xHat(3);
u0=xHat(4);
v0=xHat(5);
bx0=xHat(7);
by0=xHat(8);

a=Minv(1,1);
b=Minv(2,2);
c=Minv(3,3);
d=Minv(2,3);

A11=[0 0 -sin(psi0)*u0-cos(psi0)*v0;
     0 0  cos(psi0)*u0-sin(psi0)*v0;
     0 0              0];
 
A12=[cos(psi0) -sin(psi0)  0;
     sin(psi0)  cos(psi0)  0;
        0          0       1];

A21=[0 0 -a*sin(psi0)*bx0+a*cos(psi0)*by0;
     0 0 -b*cos(psi0)*bx0-b*sin(psi0)*by0;
     0 0 -d*cos(psi0)*bx0-d*sin(psi0)*by0];
 
A23=[ a*cos(psi0) a*sin(psi0) 0;
     -b*sin(psi0) b*cos(psi0) d;
     -d*sin(psi0) d*cos(psi0) c];
 
A=[A11 A12 zeros(3,9);
   A21 -Minv*D A23 zeros(3,6);
   zeros(3,6) -T zeros(3,6);
    zeros(6,9) Omega2]*dt+eye(15);

E=[zeros(6,6); 
    Psi zeros(3,3);
    zeros(6,3) Sigma2];

C=[eye(3) zeros(3,9) eye(3);
    zeros(3,3) eye(3) zeros(3,9)];

R=[0.1^2 0 0 0 0 0;
    0 0.1^2 0 0 0 0;
    0 0 (.01/180*pi)^2 0 0 0;
    0 0 0 0.1^2 0 0;
    0 0 0 0 0.1^2 0;
    0 0 0 0 0 (.01/180*pi)^2];

Qstd2=diag([Cbn,Cbn,Cbn, Wdn,Wdn, Wdn/180*pi]);
Q2=Qstd2.^2;
Qstd1=diag([Pn,Pn,0,Vn,Vn,Vn/180*pi,0,0,0,0,0,0,0,0,0]);
Q1=Qstd1.^2;
Q=(E*Q2*E'+Q1)*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Propagate  estimate %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xdot=SimVesselNLmodel(t,xHat,u,dt,0,0);
xPrio=xdot(1:15)*dt + xHat(1:15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Propagate Apriori covar.  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pprio=A*P*A'+Q;
S = C*Pprio*C'+R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   Propagate meas.   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

if isnan(y)
    xHat(1:15)=xPrio;
    P=Pprio;
    r=NaN*ones(size(y));
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%  Kalman Gain   %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K=Pprio*C'/(S);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%  State Update  %%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = y-C*xPrio(1:15);
    xHat=xPrio+K*(r);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%  Propagate Posteriori covar. %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P=(eye(15)-K*C)*Pprio;
    P=(P+P')/2;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%   output  %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ppost = P;
xHatUpdated=xHat;

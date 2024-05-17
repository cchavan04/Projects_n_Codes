clc;
dx=0.001;  %input dx
dt=0.00001;  %input dt
c=dt/dx;
T=0.16;
NX=1+1/dx;
NT=floor(1+T/dt);
%% Initialise Q1, Q2, Q3 (at t=0)
for it=1:NT
    for j=1:NX
        x=(j-1)*dx;
        if x<0.5
            Q1(it,j)=0.445; 
            Q2(it,j)=0.311; 
            Q3(it,j)=8.928;
        else
            Q1(it,j)=0.5; 
            Q2(it,j)=0; 
            Q3(it,j)=1.4275;
        end
    end
end
%% Godunov with MUSCL scheme
for it=1:NT-1
    t=dt*it;
    % Compute speed of sound
    for j=1:NX
        p(j)=0.4*(Q3(it,j)-0.5*((Q2(it,j))^2)/Q1(it,j));
        a(j)=sqrt(abs(1.4*p(j)/Q1(it,j)));
    end
    % MUSCL slope limiter
    S1=slopelim(Q1,it,NX,dx);
    S2=slopelim(Q2,it,NX,dx);
    S3=slopelim(Q3,it,NX,dx);
    % u(i+1/2)+-
    [QL1,QR1]=uaprox(Q1,S1,it,NX,dx);
    [QL2,QR2]=uaprox(Q2,S2,it,NX,dx);
    [QL3,QR3]=uaprox(Q3,S3,it,NX,dx);
    % Compute fluxes
    for j=1:NX-1
        EL1(it,j)=QL2(it,j);
        u(j)=QL2(it,j)/QL1(it,j);
        p(j)=0.4*(QL3(it,j)-0.5*((QL2(it,j))^2)/QL1(it,j));
        EL2(it,j)=QL2(it,j)*u(j)+p(j);
        EL3(it,j)=QL3(it,j)*u(j)+p(j)*u(j);
    end
    for j=1:NX-1
        ER1(it,j)=QR2(it,j);
        u(j)=QR2(it,j)/QR1(it,j);
        p(j)=0.4*(QR3(it,j)-0.5*((QR2(it,j))^2)/QR1(it,j));
        ER2(it,j)=QR2(it,j)*u(j)+p(j);
        ER3(it,j)=QR3(it,j)*u(j)+p(j)*u(j);
    end
    for n=1:NX-1
        w=[QL1(it,n) (QL1(it,n)+a(n)) (QL1(it,n)-a(n)) QR1(it,n) (QR1(it,n)+a(n)) (QR1(it,n)-a(n))];
        F1(n)=0.5*(EL1(it,n)+ER1(it,n))-0.5*max(w)*(QR1(it,n)-QL1(it,n));
    end
    for n=1:NX-1
        w=[QL2(it,n) (QL2(it,n)+a(n)) (QL2(it,n)-a(n)) QR2(it,n) (QR2(it,n)+a(n)) (QR2(it,n)-a(n))];
        F2(n)=0.5*(EL2(it,n)+ER2(it,n))-0.5*max(w)*(QR2(it,n)-QL2(it,n));
    end
    for n=1:NX-1
        w=[QL3(it,n) (QL3(it,n)+a(n)) (QL3(it,n)-a(n)) QR3(it,n) (QR3(it,n)+a(n)) (QR3(it,n)-a(n))];
        F3(n)=0.5*(EL3(it,n)+ER3(it,n))-0.5*max(w)*(QR3(it,n)-QL3(it,n));
    end
    % Calculate Q1, Q2, Q3 at t(n+1)
    for k=2:NX-1
        Q1(it+1,k)=Q1(it,k)-c*(F1(k)-F1(k-1));
        Q2(it+1,k)=Q2(it,k)-c*(F2(k)-F2(k-1));
        Q3(it+1,k)=Q3(it,k)-c*(F3(k)-F3(k-1));
    end
end
%% Plotting solution
figure;
d=Q1(NT,:);
v=Q2(NT,:)./Q1(NT,:);
P=0.4.*(Q3(NT,:)-0.5.*((Q2(NT,:).^2)./Q1(NT,:)));
U=(Q3(NT,:)./Q1(NT,:))-0.5.*(v.^2);
x=linspace(0,1,NX);
plot(x,d,"black");
ylim([0 1.5]);
title('Density - Godunov w/ MUSCL scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Density (\rho)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,v,"black");
ylim([0 2]);
title('Velocity - Godunov w/ MUSCL scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Velocity (u)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,P,"black");
ylim([0 5]);
title('Pressure - Godunov w/ MUSCL scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Pressure (P)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,U,"black");
ylim([0 25]);
title('Internal energy - Godunov w/ MUSCL scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Internal Energy (e)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
%% Saving solution for comparing
save('GwMUSCL_density.mat','d');
save('GwMUSCL_velocity.mat','v');
save('GwMUSCL_pressure.mat','P');
save('GwMUSCL_energy.mat','U');

% MUSCL slope limiter function
function S=slopelim(u,it,NX,dx)
for i=2:NX-1
    a=(u(it,i+1)-u(it,i))/dx;
    b=(u(it,i)-u(it,i-1))/dx;
    A=[abs(a) abs(b)];
    C=sign(a);
    M=min(A);
    if a*b>0
        S(i)=C*M;
    else
        S(i)=0;
    end
end
S(1)=(u(it,2)-u(it,1))/dx;
S(NX)=(u(it,NX)-u(it,NX-1))/dx;
end

% u(i+1/2)+- Polynomial approx. function
function [uL,uR]=uaprox(u,S,it,NX,dx)
for m=1:NX-1
    uL(it,m)=u(it,m)+S(m)*dx/2;
    uR(it,m)=u(it,m+1)-S(m+1)*dx/2;
end
end


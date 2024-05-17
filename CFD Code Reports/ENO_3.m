clc;
dx=0.01;  %input dx
dt=0.0001;  %input dt
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
%% Godunov with ENO scheme
for it=1:NT-1
    t=dt*it;
    % Compute speed of sound
    for j=1:NX
        p(j)=0.4*(Q3(it,j)-0.5*((Q2(it,j))^2)/Q1(it,j));
        a(j)=sqrt(abs(1.4*p(j)/Q1(it,j)));
    end
    % u(i+1/2)+-
    [QL1,QR1]=ENOaprox(Q1,it,NX);
    [QL2,QR2]=ENOaprox(Q2,it,NX);
    [QL3,QR3]=ENOaprox(Q3,it,NX);
    % Compute fluxes
    for j=3:NX-2
        EL1(it,j)=QL2(it,j);
        u(j)=QL2(it,j)/QL1(it,j);
        p(j)=0.4*(QL3(it,j)-0.5*((QL2(it,j))^2)/QL1(it,j));
        EL2(it,j)=QL2(it,j)*u(j)+p(j);
        EL3(it,j)=QL3(it,j)*u(j)+p(j)*u(j);
    end
    for j=3:NX-2
        ER1(it,j)=QR2(it,j);
        u(j)=QR2(it,j)/QR1(it,j);
        p(j)=0.4*(QR3(it,j)-0.5*((QR2(it,j))^2)/QR1(it,j));
        ER2(it,j)=QR2(it,j)*u(j)+p(j);
        ER3(it,j)=QR3(it,j)*u(j)+p(j)*u(j);
    end
    for n=3:NX-3
        w=[QL1(it,n) (QL1(it,n)+a(n)) (QL1(it,n)-a(n)) QR1(it,n+1) (QR1(it,n+1)+a(n)) (QR1(it,n+1)-a(n))];
        F1(n)=0.5*(EL1(it,n)+ER1(it,n+1))-0.5*max(w)*(QR1(it,n+1)-QL1(it,n));
    end
    for n=3:NX-3
        w=[QL2(it,n) (QL2(it,n)+a(n)) (QL2(it,n)-a(n)) QR2(it,n+1) (QR2(it,n+1)+a(n)) (QR2(it,n+1)-a(n))];
        F2(n)=0.5*(EL2(it,n)+ER2(it,n+1))-0.5*max(w)*(QR2(it,n+1)-QL2(it,n));
    end
    for n=3:NX-3
        w=[QL3(it,n) (QL3(it,n)+a(n)) (QL3(it,n)-a(n)) QR3(it,n+1) (QR3(it,n+1)+a(n)) (QR3(it,n+1)-a(n))];
        F3(n)=0.5*(EL3(it,n)+ER3(it,n+1))-0.5*max(w)*(QR3(it,n+1)-QL3(it,n));
    end
    % Calculate Q1, Q2, Q3 at t(n+1)
    for k=4:NX-3
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
title('Density - Godunov w/ ENO scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Density (\rho)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,v,"black");
ylim([0 2]);
title('Velocity - Godunov w/ ENO scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Velocity (u)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,P,"black");
ylim([0 5]);
title('Pressure - Godunov w/ ENO scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Pressure (P)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,U,"black");
ylim([0 25]);
title('Internal energy - Godunov w/ ENO scheme @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Internal Energy (e)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
%% Saving solution for comparing
save('GwENO_density.mat','d');
save('GwENO_velocity.mat','v');
save('GwENO_pressure.mat','P');
save('GwENO_energy.mat','U');


% ENO approximation function
function [uL,uR]=ENOaprox(u,it,NX)
% Selecting stencil with divided difference
for i=1:NX-1
    D(1,i)=u(it,i+1)-u(it,i);
end
for i=1:NX-2
    D(2,i)=D(1,i+1)-D(1,i);
end
for i=2:NX-2
    is=i;
    for m=1:2
        if abs(D(m,is-1))<abs(D(m,is))
            is=is-1;
        end
    end
    R(i)=i-is;
end
% u(i+1/2)- & u(i-1/2)+
C=[11/6 -7/6 1/3; 1/3 5/6 -1/6; -1/6 5/6 1/3; 1/3 -7/6 11/6];
for m=3:NX-2
    j=R(m)+2;
    k=R(m);
    uL(it,m)=C(j,1)*u(it,m-k)+C(j,2)*u(it,m-k+1)+C(j,3)*u(it,m-k+2);
    uR(it,m)=C(j-1,1)*u(it,m-k)+C(j-1,2)*u(it,m-k+1)+C(j-1,3)*u(it,m-k+2);
end
end
clc;
dx=0.001;
dt=0.00001;
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
%% Upwind - Flux difference splitting method
for it=1:NT-1
    t=dt*it;
    % Compute fluxes
    for j=1:NX
        E1(it,j)=Q2(it,j);
        u(j)=Q2(it,j)/Q1(it,j);
        p(j)=0.4*(Q3(it,j)-0.5*((Q2(it,j))^2)/Q1(it,j));
        a(j)=sqrt(abs(1.4*p(j)/Q1(it,j))); %speed of sound
        E2(it,j)=Q2(it,j)*u(j)+p(j);
        E3(it,j)=Q3(it,j)*u(j)+p(j)*u(j);
    end
    for n=1:NX-1
        w=[Q1(it,n) (Q1(it,n)+a(n)) (Q1(it,n)-a(n)) Q1(it,n+1) (Q1(it,n+1)+a(n)) (Q1(it,n+1)-a(n))];
        F1(n)=0.5*(E1(it,n)+E1(it,n+1))-0.5*max(w)*(Q1(it,n+1)-Q1(it,n));
    end
    for n=1:NX-1
        w=[Q2(it,n) (Q2(it,n)+a(n)) (Q2(it,n)-a(n)) Q2(it,n+1) (Q2(it,n+1)+a(n)) (Q2(it,n+1)-a(n))];
        F2(n)=0.5*(E2(it,n)+E2(it,n+1))-0.5*max(w)*(Q2(it,n+1)-Q2(it,n));
    end
    for n=1:NX-1
        w=[Q3(it,n) (Q3(it,n)+a(n)) (Q3(it,n)-a(n)) Q3(it,n+1) (Q3(it,n+1)+a(n)) (Q3(it,n+1)-a(n))];
        F3(n)=0.5*(E3(it,n)+E3(it,n+1))-0.5*max(w)*(Q3(it,n+1)-Q3(it,n));
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
title('Density - Upwind-Flux difference splitting @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Density (\rho)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,v,"black");
ylim([0 2]);
title('Velocity - Upwind-Flux difference splitting @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Velocity (u)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,P,"black");
ylim([0 5]);
title('Pressure - Upwind-Flux difference splitting @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Pressure (P)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,U,"black");
ylim([0 25]);
title('Internal energy - Upwind-Flux difference splitting @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Internal Energy (e)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
%% Saving solution for comparing
save('UFDS_density.mat','d');
save('UFDS_velocity.mat','v');
save('UFDS_pressure.mat','P');
save('UFDS_energy.mat','U');

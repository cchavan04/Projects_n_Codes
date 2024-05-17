clc;
dx=0.001;  %input dx
dt=0.0001;  %input dt
c=dt/dx;
T=0.16;
NX=1+1/dx;
NT=floor(1+T/dt);
e=0.1;
%% Initialise Q1, Q2, Q3 (at t=0)
for i=1:NT
    for j=1:NX
        x=(j-1)*dx;
        if x<0.5
            Q1(i,j)=0.445; 
            Q2(i,j)=0.311; 
            Q3(i,j)=8.928;
        else
            Q1(i,j)=0.5; 
            Q2(i,j)=0; 
            Q3(i,j)=1.4275;
        end
    end
end
%% R-K (2nd order) w/ 4th order dissipation
for i=1:NT-1
    % R-K loop 1
    for j=1:NX
        E1(i,j)=Q2(i,j);
        u(j)=Q2(i,j)/Q1(i,j);
        p(j)=0.4*(Q3(i,j)-0.5*((Q2(i,j))^2)/Q1(i,j));
        E2(i,j)=Q2(i,j)*u(j)+p(j);
        E3(i,j)=Q3(i,j)*u(j)+p(j)*u(j);
    end
    for j=2:NX-1
        Q1(i+1,j)=Q1(i,j)-c*(E1(i,j+1)-E1(i,j-1))/4;
        Q2(i+1,j)=Q2(i,j)-c*(E2(i,j+1)-E2(i,j-1))/4;
        Q3(i+1,j)=Q3(i,j)-c*(E3(i,j+1)-E3(i,j-1))/4;
    end
    % R-K loop 2
    for j=1:NX
        E1(i,j)=Q2(i+1,j);
        u(j)=Q2(i+1,j)/Q1(i+1,j);
        p(j)=0.4*(Q3(i+1,j)-0.5*((Q2(i+1,j))^2)/Q1(i+1,j));
        E2(i,j)=Q2(i+1,j)*u(j)+p(j);
        E3(i,j)=Q3(i+1,j)*u(j)+p(j)*u(j);
    end
    for j=2:NX-1
        Q1(i+1,j)=Q1(i,j)-(E1(i,j+1)-E1(i,j-1))*c/2;
        Q2(i+1,j)=Q2(i,j)-(E2(i,j+1)-E2(i,j-1))*c/2;
        Q3(i+1,j)=Q3(i,j)-(E3(i,j+1)-E3(i,j-1))*c/2;
    end
    % Artificial dissipation loop
    for j=3:NX-2
        Q1(i+1,j)=Q1(i+1,j)-e*(Q1(i,j-2)-4*Q1(i,j-1)+6*Q1(i,j)-4*Q1(i,j+1)+Q1(i,j+2));
        Q2(i+1,j)=Q2(i+1,j)-e*(Q2(i,j-2)-4*Q2(i,j-1)+6*Q2(i,j)-4*Q2(i,j+1)+Q2(i,j+2));
        Q3(i+1,j)=Q3(i+1,j)-e*(Q3(i,j-2)-4*Q3(i,j-1)+6*Q3(i,j)-4*Q3(i,j+1)+Q3(i,j+2));
    end
end
%% Plot solutions
figure;
d=Q1(NT,:);
v=Q2(NT,:)./Q1(NT,:);
P=0.4.*(Q3(NT,:)-0.5.*((Q2(NT,:).^2)./Q1(NT,:)));
U=(Q3(NT,:)./Q1(NT,:))-0.5.*(u.^2);
x=linspace(0,1,NX);
plot(x,d,"black");
ylim([0 1.5]);
title('Density - R-K w/ artificial dissipation @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Density (\rho)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,v,"black");
ylim([0 2]);
title('Velocity - R-K w/ artificial dissipation @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Velocity (u)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,P,"black");
ylim([0 5]);
title('Pressure - R-K w/ artificial dissipation @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Pressure (P)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
figure;
plot(x,U,"black");
ylim([0 25]);
title('Internal energy - R-K w/ artificial dissipation @T=0.16 (case 1)');
xlabel('Location (x)');
ylabel('Internal Energy (e)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
grid on;
%% Saving solution for comparing
save('RKwAD_density.mat','d');
save('RKwAD_velocity.mat','v');
save('RKwAD_pressure.mat','P');
save('RKwAD_energy.mat','U');

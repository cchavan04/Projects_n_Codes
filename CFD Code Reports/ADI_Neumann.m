clc;
e=1;
x=-pi;
y=-pi/2;
NX=301; % input NX
NY=151; % input NY
NX=NX+2;
NY=NY+2;
dx=2*pi/(NX-3);
dy=pi/(NY-3);
b=dx/dy;
it=0;
%% Initialise
for j=1:NY
    for i=1:NX
        u(i,j)=0;
    end
end
% Neumann BC
for i=2:NX-1
    x=-pi+(i-2)*dx;
    y=-pi/2;
    f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
    u(i,1)=(u(i-1,1)+u(i+1,1)+(b^2)*(2*u(i,2))+(f*(dx^2)))/(2*(1+(b^2)));
    y=pi/2;
    f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
    u(i,NY)=(u(i-1,NY)+u(i+1,NY)+(b^2)*(2*u(i,NY-1))+(f*(dx^2)))/(2*(1+(b^2)));
end
for j=2:NY-1
    y=-(pi/2)+(j-2)*dy;
    x=-pi;
    f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
    u(1,j)=((2*u(2,j))+(b^2)*(u(1,j-1)+u(1,j+1))+(f*(dx^2)))/(2*(1+(b^2)));
    x=pi;
    f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
    u(NX,j)=((2*u(NX-1,j))+(b^2)*(u(NX,j-1)+u(NX,j+1))+(f*(dx^2)))/(2*(1+(b^2)));
end
%% ADI iteration loop
while e>0.00001
    e=0;
    it=it+1;
    e1=0;
    e2=0;
    %% X-sweep
    for j=3:NY-2
        y=-(pi/2)+(j-2)*dy;
        for i=3:NX-2
            x=-pi+(i-2)*dx;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            c(i,j)=u(i,j);
            A(i)=1;
            B(i)=-2*(1+(b^2));
            C(i)=1;
            D(i)=-(b^2)*(u(i,j+1)+u(i,j-1))-(f*(dx^2));
        end
        B(1,2)=-2*(1+(b^2));
        B(1,NX-1)=-2*(1+(b^2));
        D(1,3)=-(b^2)*(u(3,j+1)+u(3,j-1))-(f*(dx^2))-u(2,j);
        D(1,NX-2)=-(b^2)*(u(NX-2,j+1)+u(NX-2,j-1))-(f*(dx^2))-u(NX-1,j);
        for i=3:(NX-2)
            R=A(1,i)/B(1,i-1);
            B(1,i)=B(1,i)-R*C(1,i-1);
            D(1,i)=D(1,i)-R*D(1,i-1);
        end
        D(1,NX-2)=D(1,NX-2)/B(1,NX-2);
        for i=(NX-3):-1:3
            D(1,i)=(D(1,i)-C(1,i)*D(1,i+1))/B(1,i);
        end
        for i=3:NX-2
            u(i,j)=D(1,i);
        end
    end
    %% Y-sweep
    for i=3:NX-2
        x=-pi+(i-2)*dx;
        for j=3:NY-2
            y=-(pi/2)+(j-2)*dy;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            A1(j)=b^2;
            B1(j)=-2*(1+(b^2));
            C1(j)=b^2;
            D1(j)=-u(i+1,j)-u(i-1,j)-(f*(dx^2));
        end
        B1(1,2)=-2*(1+(b^2));
        B1(1,NY-1)=-2*(1+(b^2));
        D1(1,3)=-u(i+1,3)-u(i-1,3)-(f*(dx^2))-(b^2)*u(i,2);
        D1(1,NY-2)=-u(i+1,NY-2)-u(i-1,NY-2)-(f*(dx^2))-(b^2)*u(i,NY-1);
        for j=3:(NY-2)
            R1=A1(1,j)/B1(1,j-1);
            B1(1,j)=B1(1,j)-R1*C1(1,j-1);
            D1(1,j)=D1(1,j)-R1*D1(1,j-1);
        end
        D1(1,NY-2)=D1(1,NY-2)/B1(1,NY-2);
        for j=(NY-3):-1:3
            D1(1,j)=(D1(1,j)-C1(1,j)*D1(1,j+1))/B1(1,j);
        end
        for j=3:NY-2
            u(i,j)=D1(1,j);
            e1=e1+abs(u(i,j)-c(i,j));
            e2=e2+abs(c(i,j));
        end
    end
    % Neumann BC
    for i=2:NX-1
        x=-pi+(i-2)*dx;
        y=-pi/2;
        f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
        u(i,1)=(u(i-1,1)+u(i+1,1)+(b^2)*(2*u(i,2))+(f*(dx^2)))/(2*(1+(b^2)));
        y=pi/2;
        f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
        u(i,NY)=(u(i-1,NY)+u(i+1,NY)+(b^2)*(2*u(i,NY-1))+(f*(dx^2)))/(2*(1+(b^2)));
    end
    for j=2:NY-1
        y=-(pi/2)+(j-2)*dy;
        x=-pi;
        f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
        u(1,j)=((2*u(2,j))+(b^2)*(u(1,j-1)+u(1,j+1))+(f*(dx^2)))/(2*(1+(b^2)));
        x=pi;
        f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
        u(NX,j)=((2*u(NX-1,j))+(b^2)*(u(NX,j-1)+u(NX,j+1))+(f*(dx^2)))/(2*(1+(b^2)));
    end
    e=e1/e2;
end
%% Plot pressure contours
for i=1:NX-2
    for j=1:NY-2
        u1(i,j)=u(i+1,j+1);
    end
end
x=linspace(-pi,pi,NX-2);
y=linspace(-pi/2,pi/2,NY-2);
[X,Y]=meshgrid(x,y);
figure;
contourf(X,Y,u1');
%clim([-0.7 0.1]);
colorbar;
title('Pressure contours using ADI with Neumann BC ');
xlabel('X (m)');
ylabel('Y (m)');
axis([-pi pi -pi/2 pi/2]);
set(gca,'XTick',-pi:pi/4:pi); 
set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
set(gca,'YTick',-pi/2:pi/4:pi/2); 
set(gca,'YTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
a=colorbar;
a.Label.String = 'Pressure (Pa)';
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
it %display iteration count
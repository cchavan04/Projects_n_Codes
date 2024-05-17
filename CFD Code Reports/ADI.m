clc;
e=1;
x=-pi;
y=-pi/2;
NX=301; 
NY=151; 
dx=2*pi/(NX-1);
dy=pi/(NY-1);
b=dx/dy;
it=0;
%% Initialise
for j=1:NY
    for i=1:NX
        u(i,j)=0;
    end
end
for i=1:NX
    x=-pi+(i-1)*dx;
    y=-pi/2;
    u(i,1)=-(cos(2*x)+cos(2*y))/4;
    y=pi/2;
    u(i,NY)=-(cos(2*x)+cos(2*y))/4;
end
for j=1:NY
    y=-(pi/2)+(j-1)*dy;
    x=-pi;
    u(1,j)=-(cos(2*x)+cos(2*y))/4;
    x=pi;
    u(NX,j)=-(cos(2*x)+cos(2*y))/4;
end
%% ADI iteration loop
while e>0.00001
    e=0;
    it=it+1;
    e1=0;
    e2=0;
    %% X-sweep
    for j=2:NY-1
        y=-(pi/2)+(j-1)*dy;
        for i=2:NX-1
            x=-pi+(i-1)*dx;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            c(i,j)=u(i,j);
            A(i)=1;
            B(i)=-2*(1+(b^2));
            C(i)=1;
            D(i)=-(b^2)*(u(i,j+1)+u(i,j-1))-(f*(dx^2));
        end
        B(1,1)=-2*(1+(b^2));
        B(1,NX)=-2*(1+(b^2));
        D(1,2)=-(b^2)*(u(2,j+1)+u(2,j-1))-(f*(dx^2))-u(1,j);
        D(1,NX-1)=-(b^2)*(u(NX-1,j+1)+u(NX-1,j-1))-(f*(dx^2))-u(NX,j);
        for i=2:(NX-1)
            R=A(1,i)/B(1,i-1);
            B(1,i)=B(1,i)-R*C(1,i-1);
            D(1,i)=D(1,i)-R*D(1,i-1);
        end
        D(1,NX-1)=D(1,NX-1)/B(1,NX-1);
        for i=(NX-2):-1:2
            D(1,i)=(D(1,i)-C(1,i)*D(1,i+1))/B(1,i);
        end
        for i=2:NX-1
            u(i,j)=D(1,i);
        end
    end
    %% Y-sweep
    for i=2:NX-1
        x=-pi+(i-1)*dx;
        for j=2:NY-1
            y=-(pi/2)+(j-1)*dy;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            A1(j)=b^2;
            B1(j)=-2*(1+(b^2));
            C1(j)=b^2;
            D1(j)=-u(i+1,j)-u(i-1,j)-(f*(dx^2));
        end
        B1(1,1)=-2*(1+(b^2));
        B1(1,NY)=-2*(1+(b^2));
        D1(1,2)=-u(i+1,2)-u(i-1,2)-(f*(dx^2))-(b^2)*u(i,1);
        D1(1,NY-1)=-u(i+1,NY-1)-u(i-1,NY-1)-(f*(dx^2))-(b^2)*u(i,NY);
        for j=2:(NY-1)
            R1=A1(1,j)/B1(1,j-1);
            B1(1,j)=B1(1,j)-R1*C1(1,j-1);
            D1(1,j)=D1(1,j)-R1*D1(1,j-1);
        end
        D1(1,NY-1)=D1(1,NY-1)/B1(1,NY-1);
        for j=(NY-2):-1:2
            D1(1,j)=(D1(1,j)-C1(1,j)*D1(1,j+1))/B1(1,j);
        end
        for j=2:NY-1
            u(i,j)=D1(1,j);
            e1=e1+abs(u(i,j)-c(i,j));
            e2=e2+abs(c(i,j));
        end
    end
    % Update BC
    for i=1:NX
        x=-pi+(i-1)*dx;
        y=-pi/2;
        u(i,1)=-(cos(2*x)+cos(2*y))/4;
        y=pi/2;
        u(i,NY)=-(cos(2*x)+cos(2*y))/4;
    end
    for j=1:NY
        y=-(pi/2)+(j-1)*dy;
        x=-pi;
        u(1,j)=-(cos(2*x)+cos(2*y))/4;
        x=pi;
        u(NX,j)=-(cos(2*x)+cos(2*y))/4;
    end
    e=e1/e2;
end
%% Plot pressure contours
x=linspace(-pi,pi,NX);
y=linspace(-pi/2,pi/2,NY);
[X,Y]=meshgrid(x,y);
figure;
contourf(X,Y,u');
%clim([-0.5 0.4]);
colorbar;
title('Pressure contours using ADI method');
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
%% Error value
e=0;
for j=2:NY-1
    y=-(pi/2)+(j-1)*dy;
    for i=2:NX-1
        x=-pi+(i-1)*dx;
        P(i,j)=-(cos(x.*2)+cos(y.*2))/4;
        e=e+((u(i,j)-P(i,j))^2);
    end
end
Er=(sqrt(e))/(NX*NY) %display error
it %display iteration count
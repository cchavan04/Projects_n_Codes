clc;
e=1;
x=-pi;
y=-pi/2;
NX=101; % input NX
NY=101; % input NY
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
%% Gauss-Seidel iteration loop
while e>0.00001
    e=0;
    it=it+1;
    e1=0;
    e2=0;
    for j=3:NY-2
        y=-(pi/2)+(j-2)*dy;
        for i=3:NX-2
            x=-pi+(i-2)*dx;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            c=u(i,j);
            u(i,j)=(u(i-1,j)+u(i+1,j)+(b^2)*(u(i,j-1)+u(i,j+1))+f*(dx^2))/(2*(1+(b^2)));
            e1=e1+abs(u(i,j)-c);
            e2=e2+abs(c);
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
    end
    e=e1/e2;
end
%% Plot pressure contours
% for i=1:NX-2
%     for j=1:NY-2
%         u1(i,j)=u(i+1,j+1);
%     end
% end
% x=linspace(-pi,pi,NX-2);
% y=linspace(-pi/2,pi/2,NY-2);
% [X,Y]=meshgrid(x,y);
% figure;
% contourf(X,Y,u1');
% %clim([-0.7 0.1]);
% colorbar;
% title('Pressure contours using G-S with Neumann BC');
% xlabel('X (m)');
% ylabel('Y (m)');
% axis([-pi pi -pi/2 pi/2]);
% set(gca,'XTick',-pi:pi/4:pi); 
% set(gca,'XTickLabel',{'-\pi','-3\pi/4','-\pi/2','-\pi/4','0','\pi/4','\pi/2','3\pi/4','\pi'});
% set(gca,'YTick',-pi/2:pi/4:pi/2); 
% set(gca,'YTickLabel',{'-\pi/2','-\pi/4','0','\pi/4','\pi/2'});
% a=colorbar;
% a.Label.String = 'Pressure (Pa)';
% set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.02,0.005]);
it %display iteration count
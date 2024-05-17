clc;
e=1;
x=-pi;
y=-pi/2;
it=0;
% NX=201;
for NX = [11 16 21 26]
it=it+1;
NY=101;
dx=2*pi/(NX-1);
dy=pi/(NY-1);
b=dx/dy;
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
%% Gauss-Seidel iteration loop
while e>0.00001
    e=0;
    e1=0;
    e2=0;
    for j=2:NY-1
        y=-(pi/2)+(j-1)*dy;
        for i=2:NX-1
            x=-pi+(i-1)*dx;
            f=-(cos(x).^2)+(sin(x).^2)-(cos(y).^2)+(sin(y).^2);
            c=u(i,j);
            u(i,j)=(u(i-1,j)+u(i+1,j)+(b^2)*(u(i,j-1)+u(i,j+1))+f*(dx^2))/(2*(1+(b^2)));
            e1=e1+abs(u(i,j)-c);
            e2=e2+abs(c);
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
    end
    e=e1/e2;
end
%% Plot pressure contours
% x=linspace(-pi,pi,201);
% y=linspace(-pi/2,pi/2,101);
% [X,Y]=meshgrid(x,y);
% contourf(X',Y',u);
% % clim([-0.6 0.3]);
% colorbar;
% title('Pressure contours using Gauss-Seidel method');
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
%% Error plot
e=0;
for j=2:NY-1
    y=-(pi/2)+(j-1)*dy;
    for i=2:NX-1
        x=-pi+(i-1)*dx;
        P(i,j)=-(cos(x.*2)+cos(y.*2))/4;
        e=e+((u(i,j)-P(i,j))^2);
    end
end
Er(1,it)=(sqrt(e))/(NX*NY);
DX(1,it)=dx;
end
figure;
lg=loglog(DX,Er,"Marker","diamond");
acry=0;
for k=2:it
    acry = acry+((log(Er(1,k))-log(Er(1,1)))/(log(DX(1,k))-log(DX(1,1))));
end
acry = acry/(it-1);
txt = "Order of accuracy = " + acry;
text(0.3,0.00035,txt);
title("Gauss-Seidel - Error vs \Deltax with \Deltay=" + dy + " at t=0");
xlabel('\Deltax (m)');
ylabel('Error in pressure (Pa)');
set(gca,'XMinorTick','on','YMinorTick','on','TickLength',[0.03, 0.005]);
xlim([0.25 0.65]);
ylim([0.00001 0.001]);
grid on;

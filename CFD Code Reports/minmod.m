clc;
% a=-3;
% b=-7;
% A=[abs(a) abs(b)];
% c=sign(a);
% M=min(A);
% if a*b>0
%     S=c*M;
% else
%     S=0;
% end

D0=2;
D1=0;
D2=1;
if D0<D1
    if D0<D2
        R=0;
    else
        R=2;
    end
else
    if D1<D2
        R=1;
    else
        R=2;
    end
end
R
if R==0
    u=1;
elseif R==1
    u=2;
else
    u=3;
end
u


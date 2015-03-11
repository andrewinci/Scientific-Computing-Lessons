% Matrix-vector product operator for the 2D BTCS scheme using a matrix-free
% CG methos
% Input: T = temperature
% Output: AT = matrix-vector product = left hand side of the scheme

function AT=matop2D(T)
global dt dx dy b IMAX JMAX     % define some useful global variables

AT=T;
% x-fluxes
for i=1:IMAX
    for j=1:JMAX
        if(i==1)                 % left boundary
            bp=0.5*(b(i,j)+b(i+1,j));   % average b on the right
            fp=bp*(T(i+1,j)-T(i,j))/dx; % right flux
            fm=0;                       % left flux
        elseif(i==IMAX)          % right boundary
            fp=0;                       % right flux
            bm=0.5*(b(i,j)+b(i-1,j));   % average b on the left
            fm=bm*(T(i,j)-T(i-1,j))/dx; % left flux
        else
            bp=0.5*(b(i,j)+b(i+1,j));   % average b on the right
            fp=bp*(T(i+1,j)-T(i,j))/dx; % right flux
            bm=0.5*(b(i,j)+b(i-1,j));   % average b on the left
            fm=bm*(T(i,j)-T(i-1,j))/dx; % left flux
        end
        AT(i,j)=AT(i,j)-dt/dx*(fp-fm);
    end
end

% y-fluxes
for i=1:IMAX
    for j=1:JMAX
        if(j==1)                 % left boundary
            bp=0.5*(b(i,j)+b(i,j+1));   % average b on the right
            gp=bp*(T(i,j+1)-T(i,j))/dy; % right flux
            gm=0;                       % left flux
        elseif(j==JMAX)          % right boundary
            gp=0;                       % right flux
            bm=0.5*(b(i,j)+b(i,j-1));   % average b on the left
            gm=bm*(T(i,j)-T(i,j-1))/dy; % left flux
        else
            bp=0.5*(b(i,j)+b(i,j+1));   % average b on the right
            gp=bp*(T(i,j+1)-T(i,j))/dy; % right flux
            bm=0.5*(b(i,j)+b(i,j-1));   % average b on the left
            gm=bm*(T(i,j)-T(i,j-1))/dy; % left flux
        end
        AT(i,j)=AT(i,j)-dt/dy*(gp-gm);
    end
end
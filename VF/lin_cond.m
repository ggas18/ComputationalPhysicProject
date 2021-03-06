clear all
close all
clc

tc=1; %[0,tc]
Nt=1000;
dt=tc/Nt;


%==========================================================================
% Espace                                                                 %* 
%==========================================================================
Nx=20;                                                                   %*
Ny=20;                                                                   %*
Ly=1;                                                                    %*
Lx=1;                                                                    %*
%Pas d'espace                                                            %*
dx=2*Lx/Nx;                                                              %*
dy=2*Ly/Ny;                                                              %*
                                                                         %*
x=-Lx+dx/2:dx:Lx-dx/2;                                                   %*
y=-Ly+dy/2:dy:Ly-dy/2;                                                   %*
%==========================================================================
% Cercle                                                                 %*
Rc=Lx/4;                                                                 %*
Xc=0;                                                                    %*
Yc=0;                                                                    %*
%==========================================================================

%==========================================================================
%Vitesse
%==========================================================================
u=-a./(x.^2+y.^2);
v=a./(x.^2+y.^2);
%==========================================================================

%==========================================================================
%==========================================================================
%Matrice de rigidité et solution initiale
%==========================================================================
N=Nx*Ny;
A=zeros(N,N);
U0=zeros(N,1);
for j=1:Ny
    for i=1:Nx
        k1=k(Ny,i,j);
        k2=k(Ny,i-1,j);
        k3=k(Ny,i+1,j);
        k4=k(Ny,i,j-1);
        k5=k(Ny,i,j+1);
        if( i>1 && i<Nx && j>1 && j<Ny)
            A(k1,k1)=1-2*dt/dx^2-2*dt/dy^2;
            A(k1,k2)=dt/dx^2;
            A(k1,k3)=dt/dx^2;
            A(k1,k4)=dt/dy^2;
            A(k1,k5)=dt/dy^2;
        end
        if(i==1 || i==Nx ||j==1 || j==Ny )
            A(k1,k1)=1;
            
        end
        
        if ((x(i)-Xc)^2+(y(j)-Yc)^2)<=Rc^2
            U0(k1,1)=1;
        end
    end
end


%==========================================================================
% Evolution de la solution en fonction du temps
%==========================================================================
U=U0;
U1=U0;
for t=1:Nt
    Uxy=reshape(U1,Nx,Ny);
    surf(x,y,Uxy)
    axis([-Lx Lx -Ly Ly 0 1])
    drawnow
    pause(1)
    U2=A*U1;
    U1=U2;
    U=[U U2];
end


clear all
close all
clc

tc=1; 
Nt=1000;
dt=tc/Nt;
% amplitude de la vitesse a=1 5 10
a=10;

%==========================================================================
% Espace                                                                 %* 
%==========================================================================
Nx=19;                                                                   
Ny=19;                                                                   
Ly=1;                                                                    
Lx=1;                                                                    
%Pas d'espace                                                            
dx=2*Lx/Nx;                                                              
dy=2*Ly/Ny;                                                              
                                                       
x=-Lx+dx/2:dx:Lx-dx/2;                                                
y=-Ly+dy/2:dy:Ly-dy/2;                                                 
%==========================================================================
% Cercle                                                                 
Rc=Lx/4;                                                               
Xc=Lx/2;                                                                    
Yc=Lx/2;                                                                    
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
            %Conduction
            A(k1,k1)=1-2*dt/dx^2-2*dt/dy^2;
            A(k1,k2)=dt/dx^2;
            A(k1,k3)=dt/dx^2;
            A(k1,k4)=dt/dy^2;
            A(k1,k5)=dt/dy^2;
            
            %Vitesse de convection
            ue=u(a,x(i)+dx/2,y(j));
            uw=u(a,x(i)-dx/2,y(j));
            vn=v(a,x(i),y(j)+dy/2);
            vs=v(a,x(i),y(j)-dy/2);
            %Conduction-Convection
            A(k1,k1)=A(k1,k1)-dt/(2*dx)*(ue-uw)-dt/(2*dy)*(vn-vs);
            A(k1,k2)=A(k1,k2)+dt/(2*dx)*uw;
            A(k1,k3)=A(k1,k3)-dt/(2*dx)*ue;
            A(k1,k4)=A(k1,k4)+dt/(2*dy)*vs;
            A(k1,k5)=A(k1,k5)-dt/(2*dy)*vn;
        end
        %condition au bord de Dirichlet homogène
        if(i==1 || i==Nx ||j==1 || j==Ny )
            A(k1,k1)=1;
            
        end
        % À l'interieur du cercle où il y a le colerant initialement
        if (((x(i)-Xc)^2+(y(j)-Yc)^2)<=Rc^2)
            U0(k1,1)=1;
        end
    end
end


%==========================================================================
% Evolution de la solution en fonction du temps
%==========================================================================

U1=U0;
U= zeros(N,Nt);
U(:,1)=U0;
for t=1:Nt
    U(:,t+1)=A*U(:,t);
    Uxy=reshape(U(:,t),Nx,Ny);
    surf(x,y,Uxy)
    axis([-Lx Lx -Ly Ly 0 1])
    drawnow
    pause(0.1)
   
end

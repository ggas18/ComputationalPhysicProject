clear all
close all
clc

tc=1;
Nt=1000;
dt=tc/Nt;
% amplitude de la vitesse a=1 5 10
a=5;

%==========================================================================
% Espace                                                                 %*
%==========================================================================
Nx=25;
Ny=25;
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
Yc=0;
%==========================================================================


%==========================================================================
%==========================================================================
%Matrice de rigidité et solution initiale
%==========================================================================
N=Nx*Ny;
A=zeros(N,N);
U=zeros(N,Nt+1);
for j=1:Ny
    for i=1:Nx
        k1=k(Nx,i,j);
        k2=k(Nx,i-1,j);
        k3=k(Nx,i+1,j);
        k4=k(Nx,i,j-1);
        k5=k(Nx,i,j+1);
        k22=k(Nx,i-2,j);
        k33=k(Nx,i+2,j);
        k44=k(Nx,i,j-2);
        k55=k(Nx,i,j+2);
        %Calcul du champ de vitesse.
        %Vitesse de convection
        ue=u(a,x(i)+dx/2,y(j));
        uw=u(a,x(i)-dx/2,y(j));
        vn=v(a,x(i),y(j)+dy/2);
        vs=v(a,x(i),y(j)-dy/2);
        if( i>1 && i<Nx && j>1 && j<Ny)
            %Conduction
            A(k1,k1)=1-2*dt/dx^2-2*dt/dy^2;
            A(k1,k2)=dt/dx^2;
            A(k1,k3)=dt/dx^2;
            A(k1,k4)=dt/dy^2;
            A(k1,k5)=dt/dy^2;
            
            %Conduction-Convection__quick
            if( i>2 && i<Nx-1 && j>2 && j<Ny-1)
                ue_plus=(-dt/dx)*(max(ue,0)/8);
                ue_moins=(-dt/dx)*(min(ue,0)/8);
                uw_plus=(-dt/dx)*(max(uw,0)/8);
                uw_moins=(-dt/dx)*(min(uw,0)/8);
                vs_plus=(-dt/dy)*(max(vs,0)/8);
                vs_moins=(-dt/dy)*(min(vs,0)/8);
                vn_plus=(-dt/dy)*(max(vn,0)/8);
                vn_moins=(-dt/dy)*(min(vn,0)/8);
                
                % On a deja fait la multiplication des facteurs -dt/dy*1/8
                % dans le calcul des vitesses.
                A(k1,k1)=A(k1,k1)+6*ue_plus+3*ue_moins-3*uw_plus-6*uw_moins+6*vn_plus+3*vn_moins-3*vs_plus-6*vs_moins;
                A(k1,k2)=A(k1,k2)-ue_plus-6*uw_plus-3*uw_moins;
                A(k1,k3)=A(k1,k3)+3*ue_plus+6*ue_moins+uw_moins;
                A(k1,k4)=A(k1,k4)-vn_plus+6*vs_plus-3*vs_moins;
                A(k1,k5)=A(k1,k5)+3*vn_plus+6*vn_moins+vn_moins;
                A(k1,k22)= A(k1,k22)+uw_plus;
                A(k1,k33)= A(k1,k33)-ue_moins;
                A(k1,k44)= A(k1,k44)+vs_plus;
                A(k1,k55)= A(k1,k55)-vn_moins;
            end
            %si on est presqu'au bord on fait du upwind
            if(i==2 || i==Nx-1 ||j==2 || j==Ny-1)
                A(k1,k1)=A(k1,k1)-dt/(dx)*(max(ue,0)-min(uw,0))-dt/(dy)*(max(vn,0)-min(vs,0));
                A(k1,k2)=A(k1,k2)+dt/(dx)*max(uw,0);
                A(k1,k3)=A(k1,k3)-dt/(dx)*min(ue,0);
                A(k1,k4)=A(k1,k4)+dt/(dy)*max(vs,0);
                A(k1,k5)=A(k1,k5)-dt/(dy)*min(vn,0);
            end
        end
        
        if(i==1 || i==Nx ||j==1 || j==Ny )
            A(k1,k1)=1;
            
        end
        
        if (((x(i)-Xc)^2+(y(j)-Yc)^2)<=Rc^2)
            U(k1,1)=1;
        end
    end
end


%==========================================================================
% Evolution de la solution en fonction du temps
%==========================================================================

for t=1:Nt
    U(:,t+1)=A*U(:,t);
    Uxy=reshape(U(:,t),Nx,Ny);
    surf(x,y,Uxy)
    axis([-Lx Lx -Ly Ly 0 1])
    
    title('Schema quick')
    xlabel('x');
    ylabel('y');
    zlabel('{\theta}')
    drawnow
end

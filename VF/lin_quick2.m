clear all
close all
clc

tc=1;
Nt=10000;
dt=tc/Nt;
% amplitude de la vitesse a=1 5 10
a=5;

%==========================================================================
% Espace                                                                 %*
%==========================================================================
Nx=10;
Ny=10;
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
U0=zeros(N,1);
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
            tx=dt/dx;
            ty=dt/dy;
            %Conduction-Convection__quick
            if( i>2 && i<Nx-1 && j>2 && j<Ny-1)
                
                A( k1, k1 ) = A( k1, k1 ) + tx * ( 6/8 * min( uw, 0 ) + 3/8 * max( uw, 0 ) );
                A( k1, k1 ) = A( k1, k1 ) - tx * ( 6/8 * max( ue, 0 ) + 3/8 * min( ue, 0 ) );
                A( k1, k1 ) = A( k1, k1 ) + ty * ( 6/8 * min( vs, 0 ) + 3/8 * max( vs, 0 ) );
                A( k1, k1 ) = A( k1, k1 ) - ty * ( 6/8 * max( vn, 0 ) + 3/8 * min( vn, 0 ) );
                A( k1, k2 ) = A( k1, k2 ) + tx * ( 3/8 * min( uw, 0 ) + 6/8 * max( uw, 0 ) );
                A( k1, k2 ) = A( k1, k2 ) + tx * ( 1/8 * max( ue, 0 ) );
                A( k1, k3 ) = A( k1, k3 ) - tx * ( 3/8 * max( ue, 0 ) + 6/8 * min( ue, 0 ) );
                A( k1, k3 ) = A( k1, k3 ) - tx * ( 1/8 * min( uw, 0 ) );
                A( k1, k4 ) = A( k1, k4 ) + ty * ( 3/8 * min( vs, 0 ) + 6/8 * max( vs, 0 ) );
                A( k1, k4 ) = A( k1, k4 ) + ty * ( 1/8 * max( vn, 0 ) );
                A( k1, k5 ) = A( k1, k5 ) - ty * ( 3/8 * max( vn, 0 ) + 6/8 * min( vn, 0 ) );
                A( k1, k5 ) = A( k1, k5 ) - ty * ( 1/8 * min( vs, 0 ) );
                A( k1, k22 ) = A( k1, k22 ) - tx * ( 1/8 * max( uw, 0 ) );
                A( k1, k33 ) = A( k1, k33 ) - tx * ( - 1/8 * min( ue, 0 ) );
                A( k1, k44 ) = A( k1, k44 ) - ty * ( 1/8 * max( vs, 0 ) );
                A( k1, k55 ) = A( k1, k55 ) - ty * ( - 1/8 * min( vn, 0 ) );
            end
            %si on est presqu'au bord on fait du upwind
            if(i==2 || i==Nx-1 ||j==2 || j==Ny-1)
                A( k1, k1 ) = A( k1, k1 ) - max( ue, 0 ) * tx + min( uw, 0) * tx;
                A( k1, k1 ) = A( k1, k1 ) - max( vn, 0 ) * ty + min( vs, 0) * ty;
                A( k1, k2 ) = A( k1, k2 ) + max( uw, 0 ) * tx;
                A( k1, k3 ) = A( k1, k3 ) - min( ue, 0 ) * tx;
                A( k1, k4 ) = A( k1, k4 ) + max( vs, 0 ) * ty;
                A( k1, k5 ) = A( k1, k5 ) - min( vn, 0 ) * ty;
            end
        end
       %Condition de Dirichlet au bord
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
    shading interp
    axis([-Lx Lx -Ly Ly 0 1])
    title('Schema quick')
    drawnow
    pause(0.1)
    U2=A*U1;
    U1=U2;
    U=[U U2];
end

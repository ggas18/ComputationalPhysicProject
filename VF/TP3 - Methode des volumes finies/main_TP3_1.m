clear;
close all;
clc;
%% Définition du domaine
Lx = 1;
Ly = 1;
Nx = 10;
Ny = 10;
Nt = 1000;
dx = 2*Lx / (Nx-1);
dy = 2*Ly / (Ny-1);
dt = 1 / Nt;
%
x = -Lx:dx:Lx;
y = -Ly:dy:Ly;
%
Nnoeuds = Nx * Ny;
% Coordonnées du cercle : conditions initiales
xc = 0;
yc = 0;
rc = Lx / 4;
%% Initialisation du problème et mise sous forme matricielle
txx = dt / (dx^2);
tyy = dt / (dy^2);
tx = dt / dx;
ty = dt / dy;
A = zeros( Nnoeuds, Nnoeuds );
U = zeros( Nnoeuds, 1);
for j=1:Ny
    for i=1:Nx
        % Calcul des emplacement dans la matrice globale
        k1=kk(i,j,Nx);
        k2=kk(i-1,j,Nx);
        k3=kk(i+1,j,Nx);
        k4=kk(i,j-1,Nx);
        k5=kk(i,j+1,Nx);
        % A la frontière du domaine
        if ( i == 1 || j == 1 || i == Nx || j == Ny )
            A( k1, k1 ) = 1;
        end
        % A l'interieur du domaine
        if ( i > 1 && i < Nx && j > 1 && j < Ny )
            % Termes de conduction
            A( k1, k1 ) = 1 - 2 * txx - 2 * tyy;
            A( k1, k2 ) = txx;
            A( k1, k3 ) = txx;
            A( k1, k4 ) = tyy;
            A( k1, k5 ) = tyy;
        end
        % A l'interieur du cercle central
        if ( ( x(i) - xc )^2 + ( y(j) - yc )^2 <= rc^2 )
            U( k1, 1 ) = 1;
        end
    end
end
%% Calcul itératif à l'aide d'un schéma explicite
for k=1:Nt
    U( :, k+1 ) = A * U( :, k );
    Us = reshape( U( :, k+1 ), Nx, Ny );
    surf( x, y, Us );
    shading interp
    axis( [ -Lx, Lx, -Ly, Ly, 0, 1 ] );
    drawnow;
end

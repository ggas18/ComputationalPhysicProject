clear;
close all;
clc;
%% D�finition du domaine
Lx = 1;
Ly = 1;
Nx = 15;
Ny = 15;
Nt = 1000;
dx = 2*Lx / Nx;
dy = 2*Ly / Ny;
dt = 1 / Nt;
%
x = -Lx+dx/2:dx:Lx-dx/2;
y = -Ly+dy/2:dy:Ly-dy/2;
%
Nnoeuds = Nx * Ny;
%
method='quick';
% Coordonn�es du cercle : conditions initiales
xc = Lx / 2;
yc = Ly / 2;
rc = Lx / 4;
a = 10;
%% Initialisation du probl�me et mise sous forme matricielle
txx = dt / (dx^2);
tyy = dt / (dy^2);
%
tx = dt / dx;
ty = dt / dy;
%
A = zeros( Nnoeuds, Nnoeuds );
U = zeros( Nnoeuds, 1);
for j=1:Ny
    for i=1:Nx
        % Calcul des emplacement dans la matrice globale
        k1 = kk( i, j, Nx);
        k2 = kk( i-1, j, Nx);
        k3 = kk( i+1, j, Nx);
        k4 = kk( i, j-1, Nx);
        k5 = kk( i, j+1, Nx);
        % A la fronti�re du domaine
        if ( i == 1 || j == 1 || i == Nx || j == Ny )
            A( k1, k1 ) = 1;
        end
        % A l'interieur du domaine pour les termes conductifs
        if ( i > 1 && i < Nx && j > 1 && j < Ny )
            % Termes de conduction
            A( k1, k1 ) = 1 - 2 * txx - 2 * tyy;
            A( k1, k2 ) = txx;
            A( k1, k3 ) = txx;
            A( k1, k4 ) = tyy;
            A( k1, k5 ) = tyy;
        end
        % Choix du mod�le de convection
        ue = u( a, x(i) + dx/2, y(j) );
        uw = u( a, x(i) - dx/2, y(j) );
        vn = v( a, x(i), y(j) + dy/2 );
        vs = v( a, x(i), y(j) - dy/2 );
        switch( method )
            case 'quick'
                % Sur les bords interieurs du domaine, on applique upwind
                if ( i > 1 && i < Nx && j > 1 && j < Ny )
                    if ( i == 2 || j == 2 || i == Nx-1 || j == Ny-1 )
                        A( k1, k1 ) = A( k1, k1 ) - max( ue, 0 ) * tx + min( uw, 0) * tx;
                        A( k1, k1 ) = A( k1, k1 ) - max( vn, 0 ) * ty + min( vs, 0) * ty;
                        A( k1, k2 ) = A( k1, k2 ) + max( uw, 0 ) * tx;
                        A( k1, k3 ) = A( k1, k3 ) - min( ue, 0 ) * tx;
                        A( k1, k4 ) = A( k1, k4 ) + max( vs, 0 ) * ty;
                        A( k1, k5 ) = A( k1, k5 ) - min( vn, 0 ) * ty;
                    end
                end
                % A l'interieur du domaine
                if ( i > 2 && i < Nx-1 && j > 2 && j < Ny-1 )
                    % calcul des coefficients supplementaires
                    k22 = kk( i-2, j, Nx);
                    k33 = kk( i+2, j, Nx);
                    k44 = kk( i, j-2, Nx);
                    k55 = kk( i, j+2, Nx);
                    % 
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
            otherwise
                % A l'interieur du domaine
                if ( i > 1 && i < Nx && j > 1 && j < Ny )
                    switch( method )
                        case 'linear'
                            % Termes de convection
                            A( k1, k1 ) = A( k1, k1 ) - tx/2 * ue + tx/2 * uw;
                            A( k1, k1 ) = A( k1, k1 ) - ty/2 * vn + ty/2 * vs;
                            A( k1, k2 ) = A( k1, k2 ) + tx/2 * uw;
                            A( k1, k3 ) = A( k1, k3 ) - tx/2 * ue;
                            A( k1, k4 ) = A( k1, k4 ) + ty/2 * vs;
                            A( k1, k5 ) = A( k1, k5 ) - ty/2 * vn;
                        case 'upwind'
                            A( k1, k1 ) = A( k1, k1 ) - max( ue, 0 ) * tx + min( uw, 0) * tx;
                            A( k1, k1 ) = A( k1, k1 ) - max( vn, 0 ) * ty + min( vs, 0) * ty;
                            A( k1, k2 ) = A( k1, k2 ) + max( uw, 0 ) * tx;
                            A( k1, k3 ) = A( k1, k3 ) - min( ue, 0 ) * tx;
                            A( k1, k4 ) = A( k1, k4 ) + max( vs, 0 ) * ty;
                            A( k1, k5 ) = A( k1, k5 ) - min( vn, 0 ) * ty;
                    end
                end
        end
        % A l'interieur du cercle central
        if ( ( x(i) - xc )^2 + ( y(j) - yc )^2 <= rc^2 )
            U( k1, 1 ) = 1;
        end
    end
end
%% Calcul it�ratif � l'aide d'un sch�ma explicite
for k=1:Nt
    U( :, k+1 ) = A * U( :, k );
    Us = reshape( U( :, k ), Nx, Ny );
    surf( x, y, Us );
    %shading interp
    axis( [ -Lx, Lx, -Ly, Ly, 0, 1 ] );
    drawnow;
end

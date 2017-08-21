close all
Uxy=reshape(U(:, 1000),Nx,Ny);
surf(y,x,Uxy)
axis([-Lx Lx -Ly Ly 0 1]);
title('Profil à   pas de temps après');
xlabel('x');
ylabel('y');
zlabel('{\theta}')
close all
Uxy=reshape(U(:, 1000),Nx,Ny);
surf(y,x,Uxy)
axis([-Lx Lx -Ly Ly 0 1]);
title('Profil �  pas de temps apr�s');
xlabel('x');
ylabel('y');
zlabel('{\theta}')
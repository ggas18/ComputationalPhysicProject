%% Maillage
clear all
close all
clc

a=1; %gradient de pression
mu= 0.3125;
mode=0;
echelle=1;
rect='triangle_4.mat';
switch mode
    case {0}
        load(rect);
    case {1}
        nomFichMaillage='rectangle_4_4'; % On va generer un fichier mat
        longueur=echelle;
        largeur=longueur;
        Nx=10;
        Ny=10;
        Mailleur_Simple_Poiseuille_Rect(longueur,largeur,echelle,Nx,Ny,nomFichMaillage);
    case {2}
        nomFichMaillage='triangle_4'; % On va generer un fichier triangle_4.mat
        longueur=echelle;
        largeur=longueur;
        Nx=10;
        Ny=10;
        Mailleur_Simple_Poiseuille_TriEqui(longueur,echelle,Nx,nomFichMaillage);
    case {3}
        nomFichMaillage_freefem='Cacaouette_Avec_Trous.msh';
        nomFichMaillage='Cacaouette'; % On va generer un fichier triangle_4.mat
        longueur=echelle;
        Converti_Mailleur_Freefem(nomFichMaillage_freefem,nomFichMaillage,echelle);
end

%% Construction des matrices de rigidité locales

[tab,Som,Num,NS,N,NT] = matLoc( rect );

%% Matrices de rigidités globales
A=zeros(NS,NS);
Fg=zeros(NS,1);
for i=1:NT
    for j=1:3
        for k=1:3
            A(Num(j,i),Num(k,i))=  A(Num(j,i),Num(k,i))+tab{2,i}(j,k);
            Fg(Num(j,i),1)=Fg(Num(j,i),1)+a*tab{1,i}/3;
        end
    end
end

%% Condition de Dirichlet
for i=N+1:NS
    A(i,:)=zeros(1,NS);
    A(i,i)=1;
    Fg(i,1)=0;
end

%% RESOLUTION
U=A\Fg;
%figure()
for i=1:NT
    %On prend les sommets du triangle
    N1=Num(1,i);
    N2=Num(2,i);
    N3=Num(3,i);
    
    %On prend les coordonnées et les regions des sommets
    x1=Som(1,N1);
    y1=Som(2,N1);
    region1=Som(3,N1);
    
    x2=Som(1,N2);
    y2=Som(2,N2);
    region2=Som(3,N2);
    
    x3=Som(1,N3);
    y3=Som(2,N3);
    
    region3=Som(3,N3);
    U1=U(N1);
    U2=U(N2);
    U3=U(N3);
    
    
    patch([x1 x2 x3 x1],[y1 y2 y3 y1],[U1 U2 U3 U1]);
    
    
    
end
title('Solution par approximation élements finis')
shading interp
colorbar

%% Solution théorique pour le triangle.
figure()
l=1;
for i=1:NT
    %On prend les sommets du triangle
    N1=Num(1,i);
    N2=Num(2,i);
    N3=Num(3,i);
    
    %On prend les coordonnées et les regions des sommets
    x1=Som(1,N1);
    y1=Som(2,N1);
    
    
    x2=Som(1,N2);
    y2=Som(2,N2);
    
    
    x3=Som(1,N3);
    y3=Som(2,N3);
    
    region3=Som(3,N3);
 
    U1=u_th(x1-l/2,y1,mu,a,l);
    U2=u_th(x2-l/2,y2,mu,a,l);
    U3=u_th(x3-l/2,y3,mu,a,l);
    
    
    patch([x1 x2 x3 x1],[y1 y2 y3 y1],[U1 U2 U3 U1]);
    
end
title('solution analytique');
shading interp
colorbar


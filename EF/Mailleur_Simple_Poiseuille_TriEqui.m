function Mailleur_Simple_Poiseuille_TriEqui(Longueur,Echelle,Nx,NomfichMaillage) 
% On crée un maillage formé de triangles équilatéraux pour une section qui est un triangle équilatéral.

% Dans un premier temps on ne distingue pas les sommets du bord dans la création des triangles, ce qui
% simplifie la programmation, et on créé un
% tableau d'entiers nommé "Bord" qui contient les numéros d'ordre qu'auront les points du bord une fois triés.
% Note: le tableau des sommets contient également une troisième ligne (non
% utilisée ici) qui indique par une valeur non nulle à quel bord appartient
% ce sommet ou de quel coin il s'agit (pour la signification des valeurs: voir le programme).


%  On trie ensuite le tableau des sommets de
%  manière à ce que les $N$
%  premiers sommets du tableau trié soient les sommets intérieurs et les $NS-N$ restants ceux du bord.

% $NC$ (ici $NC=3$) désigne le nombre de côtés du bord, on crée les $NC$
% listes ordonnées des sommets de chaque arête avec
% Arete{1}= Bas.
% Arete{2}= Gauche.
% Arete{3}= Droite.

%  On effectue également un tracé du maillage et un tracé des arêtes qui servent aussi de contrôles.

% En fin de programme, les résultats utiles sont sauvegardés dans un fichier (instruction "save").
% Il suffit de charger le fichier sauvegardé dans le programme utilisateur,
% par l'instruction "load", pour disposer des tableaux et variables utiles.

% ENTREES:
% Longueur= Longueur d'un côté.
% Echelle= Dimension catactéristique. les coordonnées sont adimensionnées.
% Nx= Nombre de subdivisions d'un côté.
% NomfichMaillage= Nom du fichier de sauvegarde.

% SORTIES:
% Section= Chaîne de caractères décrivant la nature de la section. 
% Echelle = est sauvegardée pour redimensionner au besoin les calculs.
% NT= Nombre total de triangles.
% NS= Nombre total de sommets.
% N= Nombre de sommets intérieurs.

% Som = Tableau à 3 lignes et NS colonnes dont les deux premières lignes indiquent les coordonnées des sommets (1 colonne=1 sommet),
%       les N premiers sommets (i.e. colonnes) sont les sommets du bord.
%       La troisième ligne est un réel si il est non nul il indique à quel bord
%       appartient le sommet, le cas échéant (0.0 si intérieur). Les coordonnées sont sans
%       dimension.
%
% Num= Tableau à 3 lignes et NT colonnes qui contient les numéros d'ordre des 3 sommets de chaque triangle (1 colonne=1 triangle),
%      l'ordre des sommets est celui du tableau "Som".
%
% Notes: Les sommets du bord dans le tableau des sommets sont
% numérotés de manière croissante dans le sens des aiguilles d'une montre alors qu'ils sont
% classés de manière croissante dans le sens trigonométrique sur chaque arête, ce qui n'a pas
% d'importance. La structure Arete et NC ne sont pas sauvegardée.

if nargin ~= 4
    error('Nombre de paramètres incorrect.')
end

if (Nx<3)
    error('Nx est trop petit.')
end

% On peut prévoir d'autres tests.

% Préparation des données.
dx=Longueur/Echelle/Nx;
dy=dx*sqrt(3.0)/2.0;
NT=Nx*Nx;
NS=(Nx+1)*(Nx+2)/2;
NC=3;

% Création des tableaux de sommets et de triangles non triés.
% Réservations
SomProvisoire=zeros(3,NS);
Bord=zeros(1,NS,'int32'); % Numéro d'ordre des points du bord.
NumProvisoire=zeros(3,NT,'int32');

% liste des sommets
k=0;
for j=1:(Nx+1)
    for i=1:Nx+2-j
        k=k+1;
        SomProvisoire(1,k)=(i-1)*dx+(j-1)*dx/2.0;
        SomProvisoire(2,k)=(j-1)*dy;
    end
end
%liste des triangles
k=0;
for j=1:Nx-1
    for i=1:Nx-j
        k=k+1;
        NumProvisoire(1,k)=(j-1)*(Nx+2)-j*(j-1)/2+i;
        NumProvisoire(2,k)=(j-1)*(Nx+2)-j*(j-1)/2+i+1;
        NumProvisoire(3,k)=j*(Nx+2)-j*(j+1)/2+i;
        k=k+1;
        NumProvisoire(1,k)=(j-1)*(Nx+2)-j*(j-1)/2+i+1;
        NumProvisoire(2,k)=j*(Nx+2)-j*(j+1)/2+i;
        NumProvisoire(3,k)=j*(Nx+2)-j*(j+1)/2+i+1;
    end
    k=k+1;
        NumProvisoire(1,k)=(j-1)*(Nx+2)-j*(j-1)/2+i+1;
        NumProvisoire(2,k)=(j-1)*(Nx+2)-j*(j-1)/2+i+2;
        NumProvisoire(3,k)=j*(Nx+2)-j*(j+1)/2+i+1;

end
% Dernier triangle
j=Nx;
k=k+1;
NumProvisoire(1,k)=(j-1)*(Nx+2)-j*(j-1)/2+1;
NumProvisoire(2,k)=(j-1)*(Nx+2)-j*(j-1)/2+2;
NumProvisoire(3,k)=j*(Nx+2)-j*(j+1)/2+1;
% Repérage des noeuds du bord.
% La troisième ligne de SomProvisoire est affectée mais pas utilisée.
% Côté Bas
for i=2:Nx
    Bord(1,i)=NS+1-i; % Bas
    SomProvisoire(3,i)=1.0; % Bas
end
% Côtés gauche et droit
for j=2:Nx
    Bord(1,(j-1)*(Nx+2)-j*(j-1)/2+1)=NS-3*Nx+j-1; % Gauche
    Bord(1,j*(Nx+2)-j*(j+1)/2)=NS-Nx-j+1; % Droit
    SomProvisoire(3,(j-1)*(Nx+2)-j*(j-1)/2+1)=2.0; % Gauche
    SomProvisoire(3,j*(Nx+2)-j*(j+1)/2)=3.0; % Droit
end

% Sommets
Bord(1,1)=NS; % Bas Gauche
Bord(1,Nx+1)=NS-Nx; %Bas Droit
Bord(1,(Nx+1)*(Nx+2)/2)=NS-2*Nx;  % Haut
SomProvisoire(3,1)=0.5; % Bas Gauche 
SomProvisoire(3,Nx+1)=1.5; %Bas droit 
SomProvisoire(3,(Nx+1)*(Nx+2)/2)=2.5;  % Haut


% Tri des sommets
Som=zeros(3,NS);
Num=zeros(3,NT,'int32');
pos1=0;
for i=1:NS
    if (Bord(1,i)==0)
        pos1=pos1+1;
        Som(1,pos1)=SomProvisoire(1,i);
        Som(2,pos1)=SomProvisoire(2,i);
        for j=1:NT
            if (NumProvisoire(1,j)==i)
                Num(1,j)=pos1;
            end
            if (NumProvisoire(2,j)==i)
                Num(2,j)=pos1;
            end
            if (NumProvisoire(3,j)==i)
                Num(3,j)=pos1;
            end
        end
    else
        pos2=Bord(1,i);
        Som(1,pos2)=SomProvisoire(1,i);
        Som(2,pos2)=SomProvisoire(2,i);
        Som(3,pos2)=SomProvisoire(3,i);
        for j=1:NT
            if (NumProvisoire(1,j)==i)
                Num(1,j)=pos2;
            end
            if (NumProvisoire(2,j)==i)
                Num(2,j)=pos2;
            end
            if (NumProvisoire(3,j)==i)
                Num(3,j)=pos2;
            end
        end
    end
end


% Création des arêtes.
% Réservations et formatages.
Arete=cell(NC);
Arete{1}=zeros(1,Nx+1,'int32');
Arete{2}=zeros(1,Nx+1,'int32');
Arete{3}=zeros(1,Nx+1,'int32');

% Affectations
Arete{1}(1)=NS; %Bas
Arete{1}(Nx+1)=NS-Nx;
Arete{2}(1)=NS-2*Nx; %Gauche
Arete{2}(Nx+1)=NS;% Gauche
Arete{3}(1)=NS-Nx; % Droite
Arete{3}(Nx+1)=NS-2*Nx; % Droite
for j=1:NC
    for i=2:Nx
        Arete{j}(i)=Arete{j}(1)-i+1;
    end
end


% tracés
Lx=Longueur/Echelle;
Ly=Lx;
N=NS-3*Nx;
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
figure(1);
axis equal;
rectangle('Position',[-.1,-.1,Lx+.2,Ly+.2]);
xlabel('x');
ylabel('y');
% title('Tracés de la triangulation en coordonnées réduites.');
title(['Tracés de la triangulation en coordonnées réduites.', 'N= ', num2str(N), ',  NS= ', num2str(NS),...
    ',   NT= ', num2str(NT)]);
X=zeros(1,4);
Y=zeros(1,4);
for i=1:NT
    X(1)=Som(1,Num(1,i));
    X(2)=Som(1,Num(2,i));
    X(3)=Som(1,Num(3,i));
    X(4)=Som(1,Num(1,i));
    Y(1)=Som(2,Num(1,i));
    Y(2)=Som(2,Num(2,i));
    Y(3)=Som(2,Num(3,i));
    Y(4)=Som(2,Num(1,i));
    line(X,Y);
end

figure('Position',[scrsz(3)/4 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);
figure(2);
axis equal;
rectangle('Position',[-.5,-.5,Lx+1,Ly+1]);
xlabel('x');
ylabel('y');
title('Tracés des arêtes en coordonnées réduites.');
A=zeros(1,2);
B=zeros(1,2);
for i=1:length(Arete{1})-1
    A(1)=Som(1,Arete{1}(i));
    A(2)=Som(1,Arete{1}(i+1));
    B(1)=Som(2,Arete{1}(i));
    B(2)=Som(2,Arete{1}(i+1));
    line(A,B,'color',[0 0 0],'Marker','*','MarkerSize',2,'LineStyle','none');
end
for i=1:length(Arete{3})-1
    A(1)=Som(1,Arete{3}(i));
    A(2)=Som(1,Arete{3}(i+1));
    B(1)=Som(2,Arete{3}(i));
    B(2)=Som(2,Arete{3}(i+1));
    line(A,B,'color',[0 0 1],'Marker','*','MarkerSize',2,'LineStyle','none');
end
for j=1:length(Arete{2})-1
    A(1)=Som(1,Arete{2}(j));
    A(2)=Som(1,Arete{2}(j+1));
    B(1)=Som(2,Arete{2}(j));
    B(2)=Som(2,Arete{2}(j+1));
    line(A,B,'color',[1 0 0],'LineStyle','-.');
end
for j=1:length(Arete{4})-1
    A(1)=Som(1,Arete{4}(j));
    A(2)=Som(1,Arete{4}(j+1));
    B(1)=Som(2,Arete{4}(j));
    B(2)=Som(2,Arete{4}(j+1));
    line(A,B,'color',[0 1 0],'LineStyle','-.');
end
% Sauvegarde et sortie.
Section='Triangle';
save(NomfichMaillage,'Section','Echelle','NT','NS','N','Som','Num');
% Effaçage des variables locales de la mémoire.
clear all;
return


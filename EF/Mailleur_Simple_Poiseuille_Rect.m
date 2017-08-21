function Mailleur_Simple_Poiseuille_Rect(Longueur,Largeur,Echelle,Nx,Ny,NomfichMaillage)
% On crée un maillage formé de triangles rectangles pour une plaque rectangulaire.

% Dans un premier temps on ne distingue pas les sommets du bord dans la création des triangles, ce qui
% simplifie la programmation, et on créé un
% tableau d'entiers nommé "Bord" qui contient les numéros d'ordre qu'auront les points du bord une fois triés.
% Note: le tableau des sommets contient également une troisième ligne (non
% utilisée ici) qui indique par une valeur non nulle à quel bord appartient
% ce sommet ou de quel coin il s'agit (pour la signification des valeurs: voir le programme).


%  On trie ensuite le tableau des sommets de
%  manière à ce que les $N$
%  premiers sommets du tableau trié soient les sommets intérieurs et les $NS-N$ restants ceux du bord.

%  $NC$ (ici $NC=4$) désigne le nombre d'arêtes du bord, on crée les $NC$
%  listes ordonnées des sommets de chaque arête avec
% Arete{1}= Bas.
% Arete{2}= Droite.
% Arete{3}= Haut.
% Arete{4}= Gauche.

%  On effectue également un tracé du maillage et un tracé des arêtes qui servent aussi de contrôles.

% En fin de programme, les résultats utiles sont sauvegardés dans un fichier (instruction "save").
% Il suffit de charger le fichier sauvegardé dans le programme utilisateur,
% par l'instruction "load", pour disposer des tableaux et variables utiles.

% ENTREES:
% Longueur= Longueur de la plaque (selon Ox).
% Largeur= Largueur de la plaque (selon Oy).
% Echelle= Dimension catactéristique. les coordonnées sont adimensionnées.
% Nx= Nombre de subdivisions de Lx. Il faut $Nx>=2$.
% Ny= Nombre de subdivisions de Ly. Il faut $Ny>=2$.
% NomfichMaillage= Nom du fichier de sauvegarde.

% SORTIES:
% Section= Chaîne de caractères décrivant la nature de la section. 
% Echelle = est sauvegardée pour redimensionner au besoin les calculs dans le programme appelant.
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

if nargin ~= 6
    error('Nombre de paramètres incorrect.')
end

if ((Nx<2) | (Ny<2))
    error('Nx (ou Ny) est trop petit.')
end

% On peut prévoir d'autres tests.

% Préparation des données.
Lx=Longueur/Echelle;
Ly=Largeur/Echelle;
dx=Lx/Nx;
dy=Ly/Ny;
NT=2*Nx*Ny;
NS=(Nx+1)*(Ny+1);
NC=4;

% Création des tableaux de sommets et de triangles non triés.
% Réservations
SomProvisoire=zeros(3,NS);
Bord=zeros(1,NS,'int32'); % Numéro d'ordre des points du bord.
NumProvisoire=zeros(3,NT,'int32');

% liste des sommets
for j=1:(Ny+1)
    for i=1:(Nx+1)
        k=(j-1)*(Nx+1)+i;
        SomProvisoire(1,k)=(i-1)*dx;
        SomProvisoire(2,k)=(j-1)*dy;
    end
end
%liste des triangles
for j=1:Ny
    for i=1:Nx
        k=2*Nx*(j-1)+2*(i-1)+1;
        NumProvisoire(1,k)=(j-1)*(Nx+1)+i;
        NumProvisoire(2,k)=(j-1)*(Nx+1)+i+1;
        NumProvisoire(3,k)=j*(Nx+1)+i;
        k=k+1;
        NumProvisoire(1,k)=(j-1)*(Nx+1)+i+1;
        NumProvisoire(2,k)=j*(Nx+1)+i+1;
        NumProvisoire(3,k)=j*(Nx+1)+i;
    end
end

% Repérage des noeuds du bord.
% La troisième ligne de SomProvisoire est affectée mais pas utilisée.
% Côtés
for i=2:Nx
    Bord(1,i)=NS+1-i; % Bas
    Bord(1,NS-i+1)=NS+1-Nx-Ny-i; % Haut
    SomProvisoire(3,i)=4.0; % Bas
    SomProvisoire(3,NS-i+1)=2.0; % Haut
end
for j=2:Ny
    Bord(1,(j-1)*(Nx+1)+1)=NS+j-2*Nx-2*Ny-1; % Gauche
    Bord(1,j*(Nx+1))=NS+1-Nx-j; % Droit
    SomProvisoire(3,(j-1)*(Nx+1)+1)=1.0; % Gauche
    SomProvisoire(3,j*(Nx+1))=3.0; % Droit
end
% Sommets
Bord(1,1)=NS; % BasGauche
Bord(1,Nx+1)=NS-Nx; %Bas droit
Bord(1,(Nx+1)*(Ny+1))=NS-Nx-Ny;  % HautDroit
Bord(1,(Nx+1)*Ny+1)=NS-2*Nx-Ny; % Haut Gauche
SomProvisoire(3,1)=0.5; % BasGauche
SomProvisoire(3,Nx+1)=3.5; %Bas droit
SomProvisoire(3,(Nx+1)*(Ny+1))=2.5;  % HautDroit
SomProvisoire(3,(Nx+1)*Ny+1)=1.5; % Haut Gauche


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
Arete{2}=zeros(1,Ny+1,'int32');
Arete{3}=zeros(1,Nx+1,'int32');
Arete{4}=zeros(1,Ny+1,'int32');

% Affectations
Arete{1}(1)=NS;
Arete{1}(Nx+1)=NS-Nx;
Arete{3}(1)=NS-Nx-Ny;
Arete{3}(Nx+1)=NS-2*Nx-Ny;
for i=2:Nx
    Arete{1}(i)=NS+1-i;
    Arete{3}(i)=NS+1-Nx-Ny-i;
end
Arete{2}(1)=NS-Nx;
Arete{2}(Ny+1)=NS-Nx-Ny;
Arete{4}(1)=NS-2*Nx-Ny;
Arete{4}(Ny+1)=NS;
for j=2:Ny
    Arete{2}(j)=NS-Nx+1-j;
    Arete{4}(j)=NS+1-2*Nx-Ny-j;
end

N=(Nx-1)*(Ny-1);
% tracés
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
figure(1);
axis equal;
rectangle('Position',[-.1,-.1,Lx+.2,Ly+.2]);
xlabel('x');
ylabel('y');
% title('Tracés de la triangulation en coordonnées réduites.');
title(['Tracés de la triangulation en coordonnées réduites.','N= ', num2str(N), ',  NS= ', num2str(NS),...
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
Section='Rectangle';
save(NomfichMaillage,'Section','Echelle','NT','NS','N','Som','Num');
% Effaçage des variables locales de la mémoire.
clear all;
return


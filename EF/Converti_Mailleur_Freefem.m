function Converti_Mailleur_Freefem(NomFich,NomfichMaillage,Echelle) 
% NomFich=fichier detype .msh créé par freefem++
% NomfichMaillage= nom du fichier de sauvegarde sans extension
fid=fopen(NomFich,'r');
C = textscan(fid,'%d %d %d', 1);
NS=C{1};
NT=C{2};
NB=C{3};
N=NS-NB;
Somm_Int=zeros(3,N);
Somm_Bord=zeros(3,NB);
Num=zeros(3,NT,'int32');
Nouveau_Num=zeros(NS,1,'int32');
i_int=0;
i_bord=0;
for i=1:NS
 C = textscan(fid,'%f %f %d', 1);  
 if (C{3} == 0)
     i_int=i_int+1;
     Somm_Int(1,i_int)=C{1};
     Somm_Int(2,i_int)=C{2};
     Somm_Int(3,i_int)=C{3};
     Nouveau_Num(i)=i_int;
 else
     i_bord=i_bord+1;
     Somm_Bord(1,i_bord)=C{1};
     Somm_Bord(2,i_bord)=C{2};
     Somm_Bord(3,i_bord)=C{3};
     Nouveau_Num(i)=N+i_bord;
 end
end
for i=1:NT
 C = textscan(fid,'%d %d %d %d', 1);
 Num(1,i)=Nouveau_Num(C{1});
 Num(2,i)=Nouveau_Num(C{2});
 Num(3,i)=Nouveau_Num(C{3});
end
Som=cat(2,Somm_Int,Somm_Bord);
fclose(fid);
Section='Cacouette_Trouée';
save(NomfichMaillage,'Section','Echelle','NT','NS','N','Som','Num');

% Tracés
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/10 scrsz(3)/2 scrsz(4)/2]);
figure(1);
axis equal;
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
% Effaçage des variables locales de la mémoire.
clear all;
return


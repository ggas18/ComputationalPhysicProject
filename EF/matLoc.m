function [tab,Som,Num,NS,N,NT] = matLoc( nomFichierMaillage )
load(nomFichierMaillage);
mu=1;
% tab(1,:): les aires des triables.
% tab(2,:): les matrices de rigidités locales
tab=cell(2,NT);
for i=1:NT
    tab{1,i}=0;
    tab{2,i}=zeros(3,3);
end

Kloc=zeros(3,3);
for i=1:NT
    %On prend les sommets du triangle
    N1=Num(1,i);
    N2=Num(2,i);
    N3=Num(3,i);
    
    %On prend les coordonnées et les regions des sommets
    x1=Som(1,N1);
    y1=Som(2,N1);
    % region1=Som(3,N1);
    
    x2=Som(1,N2);
    y2=Som(2,N2);
    % region2=Som(3,N2);
    
    x3=Som(1,N3);
    y3=Som(2,N3);
    
   % region3=Som(3,N3);
    Ligne1=[(x3-x2)^2+(y3-y2)^2 -(x3-x1)*(x3-x2)-(y3-y1)*(y3-y2) (x2-x1)*(x3-x2)+(y2-y1)*(y3-y2)];
    Ligne2=[Ligne1(2) (x3-x1)^2+(y3-y1)^2 -(x2-x1)*(x3-x1)-(y2-y1)*(y3-y1)];
    Ligne3=[Ligne1(3) Ligne2(3) (x2-x1)^2+(y2-y1)^2];
    tab{1,i}=abs((x3-x2)*(y1-y2)-(x1-x2)*(y3-y2))/2;
    Kloc(1,:)=Ligne1;
    Kloc(2,:)=Ligne2;
    Kloc(3,:)=Ligne3;
    tab{2,i}=mu/(4*tab{1,i})*(Kloc);
end
end
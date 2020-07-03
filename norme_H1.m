% //**************************************************************************
% //  Calcul de la norme H^1_0 d'une fonction discrète définie sur les centres des volumes
% //  m : un maillage
% //  f : un vecteur colonne de taille m.nb_vol
% //
% //  NB : pour DDFV, cette fonction ne calcule qu'une partie de la norme H1 !!
% //**************************************************************************

function [norme]=norme_H1(m,f)

global K;
global L;
global DKL;
global MES_K_DEB;
global MES_K_FIN;
global MES_L_DEB;
global MES_L_FIN;

  norme=0;
  gf=zeros(m.nb_are,1);
  for i=1:m.nb_are
     gf(i) =  f(m.aretes(i,K));
     if (m.aretes(i,L)>0) 
       gf(i) = gf(i) - f(m.aretes(i,L));
     end;
     gf(i)=gf(i)/m.aretes(i,DKL);
  end;
  norme=sqrt(sum(...
      (m.aretes(:,MES_K_DEB)+m.aretes(:,MES_K_FIN)+m.aretes(:,MES_L_DEB)+m.aretes(:,MES_L_FIN))...
      .* (gf(:).^2)));
end

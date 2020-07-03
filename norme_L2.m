% //**************************************************************************
% //  Calcul de la norme L2 d'une fonction discrète définie sur les centres des volumes
% //  m : un maillage
% //  f : un vecteur colonne de taille m.nb_vol
% //**************************************************************************

function [norme]=norme_L2(m,f)

global K;
global L;
global MES_K_DEB;
global MES_K_FIN;
global MES_L_DEB;
global MES_L_FIN;

  norme=0;
  for i=1:m.nb_are
     norme = norme + (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_K_FIN))* ...
	     f(m.aretes(i,K))^2;
     if (m.aretes(i,L)>0) 
        norme = norme + (m.aretes(i,MES_L_DEB)+m.aretes(i,MES_L_FIN))* ...
	     f(m.aretes(i,L))^2;
     end;
  end;
  norme=sqrt(norme);
end
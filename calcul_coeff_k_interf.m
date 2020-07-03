% //**************************************************************************
% //    Routine de calcul du coefficient de diffusion aux aretes
% //        selon la méthode choisie
% //**************************************************************************
function [ck]=calcul_coeff_k_interf(m,donnees)

global X;
global Y;
global DEB;
global FIN;
global K;
global L;
global DKsigma;
global DLsigma;
global DKL;

  
  if (~isfield(donnees,'methode')) 
    donnees.methode='exacte';
  end;

  centres_aretes=(m.sommets(m.aretes(:,DEB),:)+m.sommets(m.aretes(:,FIN),:))/2;
  switch donnees.methode
   case 'exacte'  %// On fait le calcul exact du coefficient aux centres
                  %// des aretes
    
     ck=eval_fonction(centres_aretes,donnees.coeff_k);
     
   case 'arithmetique'   %// On fait la moyenne arithmétique des valeurs du
                         %// coefficient aux centres des mailles voisines
  
     ck=eval_fonction(m.centres(m.aretes(:,K),[X Y]),donnees.coeff_k);
     
     temp=find(m.aretes(:,L)>0);
     ck(temp)=0.5*(ck(temp) ...
		   + eval_fonction(m.centres(m.aretes(temp,L),[X Y]),donnees.coeff_k));
     
   case 'harmonique'   %// On fait la moyenne harmonique pondérée des
                       %// valeurs du coefficient aux centres des mailles voisines
 
     ck=eval_fonction(m.centres(m.aretes(:,K),[X Y]),donnees.coeff_k);
     
     temp=find(m.aretes(:,L)>0);
     ckL=eval_fonction(m.centres(m.aretes(temp,L),[X Y]),donnees.coeff_k);
     
     ck(temp)=  m.aretes(temp,DKL).*ck(temp).*ckL ...
	      ./(m.aretes(temp,DKsigma).*ckL+m.aretes(temp,DLsigma).*ck(temp));
   otherwise
     disp('Vous n''avez pas choisi une methode de calcul du coeff de diffusion valable!');
     
  end;
end



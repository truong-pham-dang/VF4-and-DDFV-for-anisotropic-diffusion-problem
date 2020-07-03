% //**************************************************************************
% //    Routine de calcul du tenseur de diffusion sur un diamant
% //    selon la méthode choisie
% //**************************************************************************
function [ADNN,ADNNS,ADNSNS]=calcul_coeff_mat(m,donnees)

global DEB;
global FIN;
global K;
global L;
global DKsigma;
global DLsigma;
global N_KL_X;
global N_KL_Y;
global NS_DEBFIN_X;
global NS_DEBFIN_Y;

  if (~isfield(donnees,'methode')) 
    donnees.methode='exacte';
  end;
  
  centres_aretes=(m.sommets(m.aretes(:,DEB),:)+m.sommets(m.aretes(:,FIN),:))/2;
  
  ADxx=eval_fonction(centres_aretes,donnees.coeff_mat.Axx);
  ADxy=eval_fonction(centres_aretes,donnees.coeff_mat.Axy);
  ADyy=eval_fonction(centres_aretes,donnees.coeff_mat.Ayy);
  
  ADN_X=ADxx.*m.aretes(:,N_KL_X)+ADxy.*m.aretes(:,N_KL_Y);
  ADN_Y=ADxy.*m.aretes(:,N_KL_X)+ADyy.*m.aretes(:,N_KL_Y);
  
  ADNS_X=ADxx.*m.aretes(:,NS_DEBFIN_X)+ADxy.*m.aretes(:,NS_DEBFIN_Y);
  ADNS_Y=ADxy.*m.aretes(:,NS_DEBFIN_X)+ADyy.*m.aretes(:,NS_DEBFIN_Y);
  
  ADNN=ADN_X.*m.aretes(:,N_KL_X)+ADN_Y.*m.aretes(:,N_KL_Y);
  ADNNS=ADN_X.*m.aretes(:,NS_DEBFIN_X)+ADN_Y.*m.aretes(:,NS_DEBFIN_Y);
  ADNSNS=ADNS_X.*m.aretes(:,NS_DEBFIN_X)+ADNS_Y.*m.aretes(:,NS_DEBFIN_Y);
  
  AKxx=eval_fonction(m.centres,donnees.coeff_mat.Axx);
  AKxy=eval_fonction(m.centres,donnees.coeff_mat.Axy);
  AKyy=eval_fonction(m.centres,donnees.coeff_mat.Ayy);
  
  
  for i=1:m.nb_are
    if (m.aretes(i,L)>0) 
      
      nbK=m.aretes(i,K);
      nbL=m.aretes(i,L);
      
      AKN_X=AKxx(nbK).*m.aretes(i,N_KL_X)+AKxy(nbK).*m.aretes(i,N_KL_Y);
      AKN_Y=AKxy(nbK).*m.aretes(i,N_KL_X)+AKyy(nbK).*m.aretes(i,N_KL_Y);
      
      AKNS_X=AKxx(nbK).*m.aretes(i,NS_DEBFIN_X)+AKxy(nbK).*m.aretes(i,NS_DEBFIN_Y);
      AKNS_Y=AKxy(nbK).*m.aretes(i,NS_DEBFIN_X)+AKyy(nbK).*m.aretes(i,NS_DEBFIN_Y);
      
      AKNN=AKN_X.*m.aretes(i,N_KL_X)+AKN_Y.*m.aretes(i,N_KL_Y);
      AKNNS=AKN_X.*m.aretes(i,NS_DEBFIN_X)+AKN_Y.*m.aretes(i,NS_DEBFIN_Y);
      AKNSNS=AKNS_X.*m.aretes(i,NS_DEBFIN_X)+AKNS_Y.*m.aretes(i,NS_DEBFIN_Y);
      
      ALN_X=AKxx(nbL).*m.aretes(i,N_KL_X)+AKxy(nbL).*m.aretes(i,N_KL_Y);
      ALN_Y=AKxy(nbL).*m.aretes(i,N_KL_X)+AKyy(nbL).*m.aretes(i,N_KL_Y);
      
      ALNS_X=AKxx(nbL).*m.aretes(i,NS_DEBFIN_X)+AKxy(nbL).*m.aretes(i,NS_DEBFIN_Y);
      ALNS_Y=AKxy(nbL).*m.aretes(i,NS_DEBFIN_X)+AKyy(nbL).*m.aretes(i,NS_DEBFIN_Y);
      
      ALNN=ALN_X.*m.aretes(i,N_KL_X)+ALN_Y.*m.aretes(i,N_KL_Y);
      ALNNS=ALN_X.*m.aretes(i,NS_DEBFIN_X)+ALN_Y.*m.aretes(i,NS_DEBFIN_Y);
      ALNSNS=ALNS_X.*m.aretes(i,NS_DEBFIN_X)+ALNS_Y.*m.aretes(i,NS_DEBFIN_Y);
      
      switch donnees.methode
        
      case 'arithmetique'   %// On fait la moyenne arithmétique des valeurs du
                            %// coefficient aux centres des mailles voisines
        
        ADNN(i)=0.5*(AKNN+ALNN);
        
        ADNNS(i)=0.5*(AKNNS+ALNNS);
        
        ADNSNS(i)=0.5*(AKNSNS+ALNSNS);
	
      case 'harmonique'   %// On fait la moyenne harmonique pondérée des
                          %// valeurs du coefficient aux centres des mailles voisines
                          %// selon la méthode m-DDFV

        ADNN(i)=((m.aretes(i,DKsigma)+m.aretes(i,DLsigma))*AKNN*ALNN) ...
                /(m.aretes(i,DLsigma)*AKNN+m.aretes(i,DKsigma)*ALNN);
        
        ADNNS(i)=(m.aretes(i,DLsigma)*ALNNS*AKNN+ ...
                 +m.aretes(i,DKsigma)*AKNNS*ALNN) ...
                /(m.aretes(i,DLsigma)*AKNN+m.aretes(i,DKsigma)*ALNN);
        
        ADNSNS(i)=(m.aretes(i,DLsigma)*ALNSNS+m.aretes(i,DKsigma)*AKNSNS) ...
                      /(m.aretes(i,DKsigma)+m.aretes(i,DLsigma)) ...
                  -m.aretes(i,DLsigma)*m.aretes(i,DKsigma)...
		     /(m.aretes(i,DKsigma)+m.aretes(i,DLsigma)) ...
                     *(AKNNS-ALNNS)^2/(m.aretes(i,DLsigma)*AKNN+m.aretes(i,DKsigma)*ALNN);
        
      end
    end
  end

end

% //**************************************************************************
% //    Calcul de la fonction de courant associée à la solution
% //    du problème elliptique calculé par VF4 sur un maillage
% //    orthogonal admissible
% //    Ne fonctionne pas avec DDFV !!!
% //**************************************************************************

function [courant]=calcul_courant(m,donnees,u)

global X;
global Y;
global DEB;
global FIN;
global K;
global L;
global MES;
global DKL;
global N_KL_X;
global N_KL_Y;
global LABEL;
   
   centres_aretes=(m.sommets(m.aretes(:,DEB),:)+m.sommets(m.aretes(:,FIN),:))/2;
   coeff_diff=eval_fonction(centres_aretes,donnees.coeff_k);
   flux=zeros(m.nb_are,1);

   temp=find(m.aretes(:,LABEL)>=0); %// Aretes interieures
   flux(temp)=(u(m.aretes(temp,L))-u(m.aretes(temp,K)))./m.aretes(temp,DKL).*m.aretes(temp,MES);

   if (isfield(donnees,'beta')) 
      temp=find(m.aretes(:,LABEL)==1); %// Aretes correspondant à
                                        %// d'éventuelles fractures
      flux(temp)=flux(temp)*donnees.beta.*m.aretes(temp,DKL)./(1+ ...
						  donnees.beta*m.aretes(temp,DKL));
   end;

   temp=find(m.aretes(:,LABEL)==-1); %// aretes du bord Dirichlet
   ubord=eval_fonction(m.sommets(m.aretes(temp,DEB),[X Y]),donnees.bordD);
   flux(temp)=(ubord-u(m.aretes(temp,K)))./m.aretes(temp,DKL).*m.aretes(temp,MES);
   


%    // Pour eviter des problemes d'orientation

   signe= (m.sommets(m.aretes(:,DEB),X)-m.sommets(m.aretes(:,FIN),X)).*m.aretes(:,N_KL_Y) ...
         -(m.sommets(m.aretes(:,DEB),Y)-m.sommets(m.aretes(:,FIN),Y)).*m.aretes(:,N_KL_X);
   signe=signe./abs(signe);

   traite=zeros(m.nb_som,1);
   courant=zeros(m.nb_som,1);

   traite(1)=1;
   totaltraite=0;
   
   while (totaltraite<m.nb_som)
     for i=1:m.nb_are
      if ((traite(m.aretes(i,DEB))==1) && (traite(m.aretes(i,FIN))==0)) 
          traite(m.aretes(i,FIN))=1;
          courant(m.aretes(i,FIN))=courant(m.aretes(i,DEB))+signe(i)*flux(i);
      end;

      if ((traite(m.aretes(i,DEB))==0) && (traite(m.aretes(i,FIN))==1)) 
          traite(m.aretes(i,DEB))=1;
          courant(m.aretes(i,DEB))=courant(m.aretes(i,FIN))-signe(i)*flux(i);
      end;    
     end;

    totaltraite=sum(traite);
  end;

end
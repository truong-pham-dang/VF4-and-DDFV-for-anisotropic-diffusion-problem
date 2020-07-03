% //**************************************************************************
% //    Le schéma DDFV pour une diffusion scalaire
% //    programmé à l'aide d'opérations vectorielles
% //    Prise en compte des CL de Dirichlet non-homogène (_LABEL=-1)
% //                    des CL de Neumann homogène  (_LABEL<-1)
% //**************************************************************************


function [A,b]=const_schema_DDFV(m,donnees)

global X;
global Y;
global DEB;
global FIN;
global K;
global L;
global MES;
global DKL;
global MES_K_DEB;
global MES_K_FIN;
global MES_L_DEB;
global MES_L_FIN;
global N_KL_X;
global N_KL_Y;
global NS_DEBFIN_X;
global NS_DEBFIN_Y;
global LABEL;

  fprintf('Assemblage DDFV\n');

  %// On évalue le terme source aux centres et aux sommets
  source_centres=eval_fonction(m.centres,donnees.source);
  source_sommets=eval_fonction(m.sommets,donnees.source);
  
  centres_aretes=(m.sommets(m.aretes(:,DEB),:)+m.sommets(m.aretes(:,FIN),:))/2;
  
  %// On calcule la valeur du coefficient de diffusion au centre de l'arête primale
  coeff_diff=calcul_coeff_k_interf(m,donnees);
  
  %// On calcule la mesure des diamants
  md=m.aretes(:,MES_K_DEB)+m.aretes(:,MES_K_FIN)+m.aretes(:,MES_L_DEB)+m.aretes(:,MES_L_FIN);

  %// On calcule le produit scalaire entre la normale primale et la normale duale
  %// Ce coefficient est nul si le maillage est orthogonal
  N_KLscaNS_KL=m.aretes(:,N_KL_X).*m.aretes(:,NS_DEBFIN_X)+m.aretes(:,N_KL_Y).*m.aretes(:,NS_DEBFIN_Y);
  
  %// On calcule maintenant les coefficients qui serviront dans la construction du schéma
  flux_NN=coeff_diff(:).*0.5.*m.aretes(:,MES).*m.aretes(:,MES)./md;
  flux_NNS=coeff_diff(:).*0.5.*m.aretes(:,MES).*m.aretes(:,DKL)./md.*N_KLscaNS_KL;
  flux_NSN=flux_NNS;
  flux_NSNS=coeff_diff(:).*0.5.*m.aretes(:,DKL).*m.aretes(:,DKL)./md;
  
%   // Ces coefficients interviennent de la façon suivante
%   // Dans le flux numérique
%   // m.aretes(:,_MES)*k_sig*grad(u).N_KL=flux_NN*(u_L-u_K)+flux_NNS*(u_FIN-u_DEB)
%   // m.aretes(:,_DKL)*k_sig*grad(u).NS_KL=flux_NSN*(u_L-u_K)+flux_NSNS*(u_FIN-u_DEB)
%   
%   // On a besoin de repérer les sommets au bord pour assembler le terme source
%   // som_au_bord vaudra 1 a l'interieur et 0 au bord

  som_au_bord=ones(m.nb_som,1);                             
  temp=find(m.aretes(:,L)<=0);
  som_au_bord(m.aretes(temp,DEB))=0;
  som_au_bord(m.aretes(temp,FIN))=0;
  
  %// On évalue la donnée de Dirichlet sur tous les sommets
  
  ub_som=zeros(m.nb_som,1); 
  ub_som=eval_fonction(m.sommets(:,[X Y]),donnees.bordD);     
 
  %// Début de la construction de la matrice et du second membre
  
  indiA=[];
  indjA=[];
  valA=[];
  
  indib=[];
  indjb=[];
  valb=[];
  
  %// On cherche les arêtes du bord Dirichlet non homogène
  temp=find(m.aretes(:,LABEL)==-1);

  %// On évalue la donnée de Dirichlet au centre de ces arêtes
  ubord=eval_fonction(centres_aretes(temp,:),donnees.bordD);
  
  
  indiA=m.aretes(temp,K);
  indjA=indiA;
  valA=flux_NN(temp);                     

  indib=m.aretes(temp,K);
  valb=flux_NN(temp).*ubord - flux_NNS(temp).*ub_som(m.aretes(temp,DEB)) ...
                            + flux_NNS(temp).*ub_som(m.aretes(temp,FIN));
  
  %// Les autres arêtes du bord sont supposées correspondre a une condition
  %// de Neumann homogène qui ne nécessite aucun traitement particulier
  
  
  %// On traite maintenant les arêtes intérieures.
  
  temp=find(m.aretes(:,L)>0);
  
  %// Contributions à l'équation de la maille K
  
  indiA=[indiA;m.aretes(temp,K)];
  indjA=[indjA;m.aretes(temp,K)];
  valA=[valA;flux_NN(temp)];
  
  indiA=[indiA; m.aretes(temp,K)];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; -flux_NN(temp)];
  
  indiA=[indiA; m.aretes(temp,K)];
  indjA=[indjA; m.aretes(temp,DEB)+m.nb_vol];
  valA=[valA; flux_NNS(temp)];
  
  indiA=[indiA; m.aretes(temp,K)];
  indjA=[indjA; m.aretes(temp,FIN)+m.nb_vol];
  valA=[valA; -flux_NNS(temp)];  
  
  indib=[indib;m.aretes(:,K)];
  valb=[valb; (m.aretes(:,MES_K_DEB)+m.aretes(:,MES_K_FIN)).*source_centres(m.aretes(:,K))]; 
  
  %// Contributions à l'équation de la maille L
  
  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; -flux_NN(temp)];  

  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; flux_NN(temp)];
    
  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,DEB)+m.nb_vol];
  valA=[valA; -flux_NNS(temp)];
  
  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,FIN)+m.nb_vol];
  valA=[valA; flux_NNS(temp)];
  
  indib=[indib; m.aretes(temp,L)];
  valb=[valb; (m.aretes(temp,MES_L_DEB)+m.aretes(temp,MES_L_FIN)).*source_centres(m.aretes(temp,L))];

%   // Contributions à l'équation du sommet DEB
%   // En prenant garde aux sommets du bord
  
  indiA=[indiA; m.aretes(:,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(:,DEB)+m.nb_vol];
  valA=[valA; flux_NSNS(:).*som_au_bord(m.aretes(:,DEB))+1-som_au_bord(m.aretes(:,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,FIN)+m.nb_vol];
  valA=[valA; -flux_NSNS(temp).*som_au_bord(m.aretes(temp,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; flux_NSN(temp).*som_au_bord(m.aretes(temp,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; -flux_NSN(temp).*som_au_bord(m.aretes(temp,DEB))];
  
  %// Le second membre est différent si on est au bord ou a l'interieur
  test= (m.aretes(:,MES_K_DEB)+m.aretes(:,MES_L_DEB)).*source_sommets(m.aretes(:,DEB)) ...
           .*som_au_bord(m.aretes(:,DEB)) ...
       +(1-som_au_bord(m.aretes(:,DEB))).*ub_som(m.aretes(:,DEB)) ;
     
    
  indib=[indib; m.aretes(:,DEB)+m.nb_vol];
  valb=[valb;test ];
 
%   // Contributions à l'équation du sommet FIN
%   // En prenant garde aux sommets du bord
  
  indiA=[indiA; m.aretes(:,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(:,FIN)+m.nb_vol];
  valA=[valA; flux_NSNS(:).*som_au_bord(m.aretes(:,FIN))+1-som_au_bord(m.aretes(:,FIN))];
  
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,DEB)+m.nb_vol];
  valA=[valA; -flux_NSNS(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; -flux_NSN(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; flux_NSN(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  
  %// Le second membre est différent si on est au bord ou a l'interieur
  test= (m.aretes(:,MES_K_FIN)+m.aretes(:,MES_L_FIN)).*source_sommets(m.aretes(:,FIN)) ...
         .*som_au_bord(m.aretes(:,FIN)) ...
       +(1-som_au_bord(m.aretes(:,FIN))).*ub_som(m.aretes(:,FIN));
  
  indib=[indib; m.aretes(:,FIN)+m.nb_vol];
  valb=[valb; test];
  
  
  %// Construction proprement dite
  
  indjb=ones(size(indib));
  
  %A=sparse([indiA indjA], valA); % SCILAB style
  %b=sparse([indib indjb], valb);
  
  A=sparse(indiA, indjA, valA); % MATLAB style
  b=sparse(indib, indjb, valb);
  b=full(b);
  
end

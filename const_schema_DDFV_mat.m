% //**************************************************************************
% //    Le schéma DDFV pour une diffusion matricielle (non isotrope)
% //    programmé à l'aide d'opérations vectorielles
% //    Prise en compte des CL de Dirichlet non-homogène (_LABEL=-1)
% //                    des CL de Neumann homogène  (_LABEL<-1)
% //**************************************************************************

function [A,b]=const_schema_DDFV_mat(m,donnees)

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
global LABEL;

  fprintf('Assemblage DDFV Anisotrope\n');

  %// On évalue le terme source aux centres et aux sommets
  source_centres=eval_fonction(m.centres,donnees.source);
  source_sommets=eval_fonction(m.sommets,donnees.source);
  
  centres_aretes=( m.sommets(m.aretes(:,DEB),:)...
		  +m.sommets(m.aretes(:,FIN),:))/2;
  
  
%   // Dans le cas ou la matrice est isotrope (donnee par coeff_k)
%   // On construit explicitement la matrice dans coeff_mat
  
  if (~isfield(donnees,'coeff_mat')) 
% Truong ajoute 06/05/2018
% //**************************************************************************
% SCILAB inline function as struct field, which is not allowed in MATLAB
%     deff('[z]=temp_axx(x,y)','z=donnees.coeff_k(x,y)');
%     deff('[z]=temp_axy(x,y)','z=0');
%     deff('[z]=temp_ayy(x,y)','z=donnees.coeff_k(x,y)');
%     donnees.coeff_mat=struct('Axx',temp_axx,'Axy',temp_axy,'Ayy',temp_ayy);

% MATLAB implementation:
     temp_axx = donnees.coeff_k;
     temp_axy = @(x,y) 0. ;
     temp_ayy = donnees.coeff_k;
     donnees.coeff_mat=struct('Axx',temp_axx,'Axy',temp_axy,'Ayy',temp_ayy);
% //**************************************************************************
% End of Truong ajoute 06/05/2018
  end;
  
%   // On calcule la valeur du coefficient de diffusion
%   // par l'une des méthodes disponibles
  
    [ADNN,ADNNS,ADNSNS]=calcul_coeff_mat(m,donnees); 
    
  
%   // On calcule la mesure des diamants
  md=m.aretes(:,MES_K_DEB)+m.aretes(:,MES_K_FIN)+m.aretes(:,MES_L_DEB)+m.aretes(:,MES_L_FIN);
  
%   // On calcule maintenant les coefficients qui serviront dans la construction du schéma
  flux_NN=0.5.*m.aretes(:,MES).*m.aretes(:,MES).*ADNN./md;
  flux_NNS=0.5.*m.aretes(:,MES).*m.aretes(:,DKL).*ADNNS./md;
  flux_NSN=flux_NNS;
  flux_NSNS=0.5.*m.aretes(:,DKL).*m.aretes(:,DKL).*ADNSNS./md;
  
%   // Ces coefficients interviennent de la façon suivante :
%   // dans la définition du schéma numérique
%   // m.aretes(:,_MES)*AD*grad(u).N_KL=flux_NN*(u_L-u_K)+flux_NNS*(u_FIN-u_DEB)
%   // m.aretes(:,_DKL)*AD*grad(u).NS_KL=flux_NSN*(u_L-u_K)+flux_NSNS*(u_FIN-u_DEB)
  
%   // On a besoin de repérer les sommets au bord pour assembler le terme source
%   // som_au_bord vaudra 1 a l'interieur et 0 au bord

  som_au_bord=ones(m.nb_som,1);                             
  temp=find(m.aretes(:,L)<=0);
  som_au_bord(m.aretes(temp,DEB))=0;
  som_au_bord(m.aretes(temp,FIN))=0;
  
%   // On évalue la donnée de Dirichlet sur tous les sommets
  
  ub_som=zeros(m.nb_som,1); 
  ub_som=eval_fonction(m.sommets(:,[X Y]),donnees.bordD);     
 
%   // Début de la construction de la matrice et du second membre
  
  indiA=[];
  indjA=[];
  valA=[];
  
  indib=[];
  indjb=[];
  valb=[];
  
%   // On cherche les arêtes du bord Dirichlet non homogène
  temp=find(m.aretes(:,LABEL)==-1);

%   // On évalue la donnée de Dirichlet au centre de ces arêtes
  ubord=eval_fonction(centres_aretes(temp,:),donnees.bordD);
  
  
  indiA=m.aretes(temp,K);
  indjA=indiA;
  valA=flux_NN(temp);                     

  indib=m.aretes(temp,K);
  valb=flux_NN(temp).*ubord - flux_NNS(temp).*ub_som(m.aretes(temp,DEB)) ...
                            + flux_NNS(temp).*ub_som(m.aretes(temp,FIN));
  
%   // Les autres arêtes du bord sont supposées correspondre a une condition
%   // de Neumann homogène qui ne nécessite aucun traitement particulier
  
  
%   // On traite maintenant les arêtes intérieures.
  
  temp=find(m.aretes(:,L)>0);
  
%   // Contributions à l'équation de la maille K
  
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
 
%   // Contributions à l'équation de la maille L
  
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
  valA=[valA; flux_NSNS(:).*som_au_bord(m.aretes(:,DEB))+ones(m.nb_are,1)-som_au_bord(m.aretes(:,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,FIN)+m.nb_vol];
  valA=[valA; -flux_NSNS(temp).*som_au_bord(m.aretes(temp,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; flux_NSN(temp).*som_au_bord(m.aretes(temp,DEB))];
  
  indiA=[indiA; m.aretes(temp,DEB)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; -flux_NSN(temp).*som_au_bord(m.aretes(temp,DEB))];
  
%   // Le second membre est différent si on est au bord ou a l'interieur
  test= (m.aretes(:,MES_K_DEB)+m.aretes(:,MES_L_DEB)).*source_sommets(m.aretes(:,DEB)) ...
           .*som_au_bord(m.aretes(:,DEB)) ...
       +(1-som_au_bord(m.aretes(:,DEB))).*ub_som(m.aretes(:,DEB)) ;
     
    
  indib=[indib; m.aretes(:,DEB)+m.nb_vol];
  valb=[valb;test ];
  
%   // Contributions à l'équation du sommet FIN
%   // En prenant garde aux sommets du bord
  
  indiA=[indiA; m.aretes(:,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(:,FIN)+m.nb_vol];
  valA=[valA; flux_NSNS(:).*som_au_bord(m.aretes(:,FIN))+ones(m.nb_are,1)-som_au_bord(m.aretes(:,FIN))];
  
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,DEB)+m.nb_vol];
  valA=[valA; -flux_NSNS(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; -flux_NSN(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  indiA=[indiA; m.aretes(temp,FIN)+m.nb_vol];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; flux_NSN(temp).*som_au_bord(m.aretes(temp,FIN))];
  
  
%   // Le second membre est différent si on est au bord ou a l'interieur
  test= (m.aretes(:,MES_K_FIN)+m.aretes(:,MES_L_FIN)).*source_sommets(m.aretes(:,FIN)) ...
         .*som_au_bord(m.aretes(:,FIN)) ...
       +(1-som_au_bord(m.aretes(:,FIN))).*ub_som(m.aretes(:,FIN));
  
  indib=[indib; m.aretes(:,FIN)+m.nb_vol];
  valb=[valb; test];
  
  
%   // Construction proprement dite
  
  indjb=ones(size(indib));
  
  A=sparse(indiA, indjA, valA); % MATLAB style
  b=sparse(indib, indjb, valb);
  b=full(b);
  
end
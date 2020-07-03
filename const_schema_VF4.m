% //**************************************************************************
% //    Le schéma VF4 
% //    programmé à l'aide d'opérations vectorielles pour des raisons d'efficacité 
% //    Prise en compte des CL de Dirichlet non-homogène (_LABEL=-1)
% //                    des CL de Neumann homogène  (_LABEL<-1)
% //    Prise en compte de fractures à l'intérieur du domaine (_LABEL=1)
% //**************************************************************************

function [A,b]=const_schema_VF4(m,donnees)

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
  
  fprintf('Assemblage VF4\n');
  
  %// On évalue le terme source aux centres des mailles
  
  source=eval_fonction(m.centres,donnees.source);
  
  %// On calcule le coefficient de diffusion (scalaire) sur chaque arête
  
  coeff_diff=calcul_coeff_k_interf(m,donnees);
  
  indiA=[];
  indjA=[];
  valA=[];
  
  indib=[];
  indjb=[];
  valb=[];
  
  %// Calcul de la transmissivite de chaque arête 
  
  tauKL=coeff_diff(:).*m.aretes(:,MES)./m.aretes(:,DKL);

  %// On corrige la transmissivite en cas de fissures dans le maillage
  if (isfield(donnees,'beta')) 
    
       %// Recherche des aretes labellisees !  
             temp=find(m.aretes(:,LABEL)==1);
             tauKL(temp)=donnees.beta*tauKL(temp).*m.aretes(temp,DKL)...
		          ./(1+donnees.beta*m.aretes(temp,DKL));
 
  end;
  
  %/// Les aretes du bord

  %/// Recherche des arêtes pour la condition de Dirichlet non homogène
  temp=find(m.aretes(:,LABEL)==-1);

  %// Evaluation des données au bord
  ubord=eval_fonction(...
      0.5*(m.sommets(m.aretes(temp,DEB),[X Y])+m.sommets(m.aretes(temp,FIN),[X Y]))...
      ,donnees.bordD);

  indiA=m.aretes(temp,K);
  indjA=indiA;
  valA=tauKL(temp);                     

  indib=[indib; m.aretes(temp,K)];
  valb=[valb; tauKL(temp).*ubord];


  %/// Toutes les autres arêtes du bord sont supposées 
  %/// correspondre à des CL de Neumann non homogène
  %/// Donc elles ne contribuent pas au schéma
  
  
  %/// La partie du terme source qui concerne tout le monde
  
  indib=[indib; m.aretes(:,K)];
  valb=[valb; (m.aretes(:,MES_K_DEB)+m.aretes(:,MES_K_FIN)).*source(m.aretes(:,K))];

  %/// On traite maintenant les aretes interieures
  
  temp=find(m.aretes(:,L)>0);

  indiA=[indiA; m.aretes(temp,K)];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; tauKL(temp)];  
  
  indiA=[indiA; m.aretes(temp,K)];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; -tauKL(temp)];
  
  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,K)];
  valA=[valA; -tauKL(temp)];  

  indiA=[indiA; m.aretes(temp,L)];
  indjA=[indjA; m.aretes(temp,L)];
  valA=[valA; tauKL(temp)];  
  
  indib=[indib; m.aretes(temp,L)];
  valb=[valb; (m.aretes(temp,MES_L_DEB)+m.aretes(temp,MES_L_FIN)).*source(m.aretes(temp,L))];


  %/// Construction proprement dite
  
  indjb=ones(size(indib));  
  A=sparse(indiA, indjA, valA); % MATLAB style
  b=sparse(indib, indjb, valb);
  b=full(b);
  
end
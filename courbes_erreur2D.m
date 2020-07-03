%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%           Divers schemas de volumes finis        %%%
%%%               pour pbs elliptiques 2D            %%%
%%%                                                  %%%
%%%                 Programme principal              %%%
%%%           permettant le calcul par VF4/DDFV      %%%
%%%                                                  %%%
%%%                F. Boyer  &  S. Krell             %%%
%%%                     Mars 2010                    %%%
%%%                                                  %%%
%%%            Version MATLAB par DANG Truong        %%%
%%%                     April 2018                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**************************************************************************
% On efface toutes les variables et on charge les fonctions utiles
%**************************************************************************
clear all
close all
clc
run('maillages2D');
run('castest2D');

%**************************************************************************
% Saisie par l'utilisateur des données du calcul
%**************************************************************************

% choix du cas test


fprintf('========================================\n');
fprintf('Choix du cas test :\n');
for i=1:length(cas_test)
  fprintf('%d) %s\n',i,cas_test{i}.nom);
end;

choix=-1;
while (int16(choix)<=0) || (int16(choix)>length(cas_test))
  choix=input('Faites votre choix :');
end;
choix=int16(choix);

% Chargement des donnees

donnees=cas_test{choix};


% On fait un choix parmi les possibilités restantes

fprintf('========================================\n');
fprintf('Choix du maillage :\n');
for i=1:length(liste_maillages)
  fprintf('%d) %s\n',i,liste_maillages(i).nom);
end;

choix=-1;
while (int16(choix)<=0) || (int16(choix)>length(liste_maillages))
  choix=input('Faites votre choix :');
end;

choix=int16(choix);
maillage=liste_maillages(choix);


% choix des différentes tailles de maillage pour le trace des courbes
if (isfield(maillage,'indice_min')) 
  fprintf('========================================\n');
  fprintf('Niveaux de raffinement utilises entre %d et %d :\n',maillage.indice_min,maillage.indice_max);
  nb_dep=-1;
  while (int16(nb_dep)<maillage.indice_min) || (int16(nb_dep)>maillage.indice_max)
    nb_dep=input('Choix du niveau le plus grosier :');
  end;
  nb_fin=-1;
  while (int16(nb_fin)<nb_dep) || (int16(nb_fin)>maillage.indice_max)
    nb_fin=input('Choix du niveau le plus fin :');
  end;
end
% si le coeff de diffusion est variable : choix de la methode

if (isfield(donnees,'methode')) 
  fprintf('========================================\n');
  fprintf('Choix de la méthode de calcul du coeff de diffusion :\n');
  fprintf('1) exacte\n');
  fprintf('2) moyenne arithmetique\n');
  fprintf('3) moyenne harmonique\n');

  choix=-1;
  while (int16(choix)<=0) || (int16(choix)>3)
    choix=input('Faites votre choix :');
  end;
  
  switch choix
   case 1
    donnees.methode='exacte';
   case 2
    donnees.methode='arithmetique';
   case 3
    donnees.methode='harmonique';
  end;
end;
 
NN=nb_fin-nb_dep+1;

errinf = zeros(NN,1);
errL2  = zeros(NN,1);
errH1  = zeros(NN,1);
pas    = zeros(NN,1);

errL2_Pascal = zeros(NN,1);

% choix du schéma : VF4 ou DDFV
if (isfield(donnees,'coeff_mat')) 
  choix=3;
  fprintf('========================================\n');
  fprintf('Schéma DDFV\n');
else
  fprintf('========================================\n');
  fprintf('Choix du schéma à utiliser :\n');
  fprintf('1) VF4\n');
  fprintf('2) DDFV\n');
  fprintf('3) DDFV hétérogène anisotrope\n');
  choix=-1;
  while ((int16(choix)<1) || (int16(choix)>3))
    choix=input('Faites votre choix :');
  end;
end;
fprintf('========================================\n');


for i=nb_dep:nb_fin     % boucle sur les différents maillages
  
  % Chargement du maillage
  nom_maillage=strcat(maillage.nom,'_',num2str(i));
  nom = strcat(rep_maillages, nom_maillage);
  m = lecture_maillage(nom);
  
  fprintf('Calcul sur le maillage %s \n',nom_maillage);
  
  % Calcul de la solution exacte aux centres des volumes de contrôle
  if (isfield(donnees,'uexacte')) 
   solexacte=eval_fonction(m.centres,donnees.uexacte);
   solexacte_sommet=eval_fonction(m.sommets,donnees.uexacte);
  end;
  
  %**************************************************************************
  % Construction du schéma et résolution du système linéaire
  %**************************************************************************

  % Construction de la matrice et du second membre

  switch choix
   case 1
     [A,b]=const_schema_VF4(m,donnees);
   case 2
     [A,b]=const_schema_DDFV(m,donnees);
   case 3
    [A,b]=const_schema_DDFV_mat(m,donnees);
  end;

  fprintf('Resolution\n');
  sol=A\b;

  % Calcul des erreurs et de la taille du pas du maillage
  
  errinf(i+1-nb_dep)=max(abs(sol(1:m.nb_vol)-solexacte));
  errL2(i+1-nb_dep)=norme_L2(m,sol(1:m.nb_vol)-solexacte);
  errH1(i+1-nb_dep)=norme_H1(m,sol(1:m.nb_vol)-solexacte);
  pas(i+1-nb_dep)=max(m.aretes(:,MES));
  
  errL2_Pascal(i+1-nb_dep)=norme_L2_Pascal(m,sol(1 : m.nb_vol)-solexacte,...
      sol(m.nb_vol+1 : m.nb_som+m.nb_vol)-solexacte_sommet,...
      solexacte,solexacte_sommet);
  
end


%**************************************************************************
%     Tracé des courbes d'erreur
%**************************************************************************

figure   %MATLAB

% Erreur en norme infinie

subplot(1,3,1)
c1=regress(log(errinf),[log(pas) ones(size(pas))]);
loglog(pas,errinf,'-+b');
title({ donnees.nom, m.nom, 'Ordre de convergence norme infinie =',num2str(c1(1))});


% Erreur en norme L2

subplot(1,3,2)
c2=regress(log(errL2),[log(pas) ones(size(pas))]);
loglog(pas,errL2,'-+b');
title({ donnees.nom, m.nom,'Ordre de convergence norme L2 = ',num2str(c2(1))});


% Erreur en norme H1

subplot(1,3,3)
c3=regress(log(errH1),[log(pas) ones(size(pas))]);
loglog(pas,errH1,'-+b');
title({ donnees.nom, maillage.nom, 'Ordre de convergence norme H1 = ',num2str(c3(1))});

figure 

% Erreur en norme L2 (Komla Domelevo et Pascal Omnes)
c4 = regress(log(errL2_Pascal), [log(pas) ones(size(pas))]);
loglog(pas,errL2_Pascal,'-+b');
title({ donnees.nom, m.nom,'Ordre de convergence norme L2 (Komla Domelevo et Pascal Omnes) = ',num2str(c4(1))});
fig_nom = strcat(m.nom,' norme L2 Pascal.png');
saveas(gcf,fig_nom);

fprintf('========================================\n');
fprintf('Ordre de convergence en norme infinie : %g\n',c1(1));
fprintf('Ordre de convergence en norme L2 : %g\n',c2(1));
fprintf('Ordre de convergence en norme H1 : %g\n',c3(1));
fprintf('==============  FIN  ===================\n');

%**************************************************************************
%           FIN
%**************************************************************************
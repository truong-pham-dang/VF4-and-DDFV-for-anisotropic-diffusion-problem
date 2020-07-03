% //////////////////////////////////////////////////////////////
% /////           Divers schemas de volumes finis        ///////
% /////               pour pbs elliptiques 2D            ///////
% /////                                                  ///////
% /////                 Programme principal              ///////
% /////           permettant le calcul par VF4/DDFV      ///////
% /////                                                  ///////
% /////                F. Boyer  &  S. Krell             ///////
% /////                     Mars 2010                    ///////
% /////                                                  ///////
% /////            Version MATLAB par DANG Truong        ///////
% /////                     April 2018                   ///////
% //////////////////////////////////////////////////////////////
%
% //**************************************************************************
% // On efface toutes les variables et on charge les fonctions utiles
% //**************************************************************************
clear all
close all
clc
run('maillages2D');
run('castest2D');

% //**************************************************************************
% // Saisie par l'utilisateur des données du calcul
% //**************************************************************************
% 


% // choix du cas test

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

% // Chargement des donnees

donnees=cas_test{choix};

% // On fait un choix parmi les maillages possibles

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

% // choix des différentes tailles de maillage pour le trace des courbes
if (isfield(maillage,'indice_min')) 
  fprintf('========================================\n');
  fprintf('Choix du niveau de raffinement entre %d et %d :\n',maillage.indice_min,maillage.indice_max);
  choix=-1;
  while (int16(choix)<maillage.indice_min) || (int16(choix)>maillage.indice_max)
    choix=input('Faites votre choix :');
  end;
  nom_maillage=strcat(maillage.nom,'_',num2str(choix));
else
  nom_maillage=maillage.nom;
end;

% // si le coeff de diffusion est variable : choix de la methode

if (isfield(donnees,'methode') )
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


% // choix du schéma : VF4 ou DDFV
if (isfield(donnees,'coeff_mat') )
  choix_schema=3;
  fprintf('========================================\n');
  fprintf('Schema DDFV\n');
else  
  fprintf('========================================\n');
  fprintf('Choix du schéma à utiliser :\n');
  fprintf('1) VF4\n');
  fprintf('2) DDFV\n');
  fprintf('3) DDFV hétérogène anisotrope\n');
  choix_schema=-1;
  while ((choix_schema<1) || (choix_schema>3))
    choix_schema=input('Faites votre choix :');
  end;
end;

fprintf('========================================\n');

nom = strcat(rep_maillages, nom_maillage);

m = lecture_maillage(nom);

if (isfield(donnees,'uexacte')) 
   solexacte=eval_fonction(m.centres,donnees.uexacte);
end;

switch choix_schema
  case 1
    [A,b]=const_schema_VF4(m,donnees);
  case 2
    [A,b]=const_schema_DDFV(m,donnees);
  case 3
    [A,b]=const_schema_DDFV_mat(m,donnees);
end;


fprintf('Resolution\n');
sol=A\b;


trace_fonction(m,sol(1:m.nb_vol),0);

if (isfield(donnees,'uexacte'))
  fprintf('========================================\n');
  trace_fonction(m,solexacte,1);
  fprintf('Erreur en norme L2 = %e\n',norme_L2(m,sol(1:m.nb_vol)-solexacte));
  fprintf('Erreur en norme Linf = %e\n',max(abs(sol(1:m.nb_vol)-solexacte)));
end;

% // Si necessaire (et si VF4 !) on calcule et on trace la fonction de courant

if (isfield(donnees,'courant') && donnees.courant>0 && choix_schema==1) 
  fprintf('========================================\n');
  fprintf('Calcul puis trace des lignes de courant\n');
  courant=calcul_courant(m,donnees,sol);

  Nx=30;
  Ny=30;

  phi=interpole(m,courant,Nx,Ny);
  scf(0);
  xset('fpf',' ');
  min_courant=min(courant)+0.01*(max(courant)-min(courant));
  max_courant=max(courant)-0.01*(max(courant)-min(courant));
  contour2d(linspace(0,1,Nx),linspace(0,1,Ny),phi,...
	    linspace(min_courant,max_courant,donnees.courant));
	    %style=ones(donnees.courant,1));
end;


fprintf('==============  FIN  ===================\n');

% //**************************************************************************
% //           FIN
% //**************************************************************************

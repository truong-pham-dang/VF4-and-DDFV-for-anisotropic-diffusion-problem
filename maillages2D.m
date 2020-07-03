% //////////////////////////////////////////////////////////////
% /////           Divers schemas de volumes finis        ///////
% /////               pour pbs elliptiques 2D            ///////
% /////                                                  ///////
% /////        Definition et outils pour les maillages   ///////
% /////                                                  ///////
% /////                F. Boyer  &  S. Krell             ///////
% /////                     Mars 2010                    ///////
% /////                                                  ///////
% /////            Version MATLAB par DANG Truong        ///////
% /////                     April 2018                   ///////
% //////////////////////////////////////////////////////////////

rep_maillages='./maillages/';   %// Le répertoire dans lequel sont
                                %// stockés les maillages

% //**************************************************************************
% //    Le format des maillages
% //**************************************************************************
% //  Chaque maillage est une structure contenant 5 elements :
% //  *  maillage.nom     :  le nom de la famille de maillages
% //  *  indice_min et indice_max : les niveaux de raffinement maximum disponibles
% //               si ces champs ne sont pas définis, il s'agit d'un maillage unique
% //  *  maillage.nb_vol  :  le nombre de volumes de contrôle
% //  *  maillage.nb_are  :  le nombre d'aretes
% //  *  maillage.nb_som  :  le nombre de sommets
% //  *  maillage.centres :  matrice de taille nb_volx2 contenant les coordonnées
% //                         des centres des volumes de contrôle
% //  *  maillage.sommets :  matrice de taille nb_somx2 contenant les coordonnées
% //                         des sommets
% //  *  maillage.aretes  :  matrice de taille nb_somx17 contenant des infos sur
% //                         les aretes / diamants
% //                         On décrit ci-dessous la signification de chaque colonne
% //**************************************************************************
% 
% // Signification des colonnes
% //   Pour la lisibilité du code, il faut utiliser ces constantes prédéfinies
% //   et non pas les valeurs numériques !!!!

global X;
global Y;
global DEB;
global FIN;
global K;
global L;
global MES;
global DKsigma;
global DLsigma;
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

X=1;   %// Coordonnée x         
Y=2;   %// Coordonnée y

DEB=1; %// Numéro du premier sommet de l'arete
FIN=2; %// Numéro du second sommet de l'arete

K=3;   %// Numéro du premier volume voisin
L=4;   %// Numéro du second volume voisin (par convention ceci est <=0
        %// pour les aretes du bord

MES=5; %// Mesure de l'arete

DKsigma=6;  %// Distance d_K,sigma
DLsigma=7;  %// Distance d_L,sigma
DKL=8;      %// Distance d_K,L

MES_K_DEB=9;  %// Mesure du quart de diamant commun à K et DEB
MES_K_FIN=10; %// Mesure du quart de diamant commun à K et FIN
MES_L_DEB=11; %// Mesure du quart de diamant commun à L et DEB
MES_L_FIN=12; %// Mesure du quart de diamant commun à L et FIN
               %// La somme de ces quatre aires donne l'aire totale du diamant

N_KL_X=13;    %// Coordonnée x de la normale unitaire à l'arete dirigée
               %// de K vers L
N_KL_Y=14;    %// Coordoonée y 

NS_DEBFIN_X=15;  %// Coordonnée x de la normale unitaire à l'arete
                  %// duale dirigée de DEB vers FIN
NS_DEBFIN_Y=16;  %// Coordonnée y

LABEL=17;       %// Label associée à l'arete
                 %// -1 : les aretes du bord dirichlet
                 %//  <-1 : le reste du bord
                 %//  1 : les aretes correspondant aux fractures


% Alternative way to list() in SCILAB: array of struct in MATLAB

for i = 1:6
   liste_maillages(i) = struct('nom',[],'indice_min',[],'indice_max',[]); 
end

% /////////////////////////////////////////////////////////
% /////   Les maillages du carre unite avec CL de Dirichlet
% /////////////////////////////////////////////////////////

maillages_carre_Dir = cell(6,1);
for i=1:6
    maillages_carre_Dir{i} = ['Sample Text ' num2str(i)];
end

index = 0;

liste_maillages(index+1)=struct('nom','maillage_rectangle_uniforme',...
			   'indice_min',1,...
			   'indice_max',6);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;
index = index + 1;

liste_maillages(index+1)=struct('nom','maillage_rectangle_bidomaine',...
			   'indice_min',1,...
			   'indice_max',6);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;
index = index + 1;

liste_maillages(index+1)=struct('nom','maillage_local_raffine',...
			   'indice_min',1,...
			   'indice_max',5);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;
index = index + 1;

liste_maillages(index+1)=struct('nom','maillage_triangle',...
			   'indice_min',1,...
			   'indice_max',6);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;
index = index + 1;
			   
liste_maillages(index+1)=struct('nom','maillage_deforme',...
			   'indice_min',1,...
			   'indice_max',5);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;	
index = index + 1;
			   
liste_maillages(index+1)=struct('nom','maillage_triangle_centregravite',...
			   'indice_min',1,...
			   'indice_max',6);
maillages_carre_Dir{index+1}=liste_maillages(index+1).nom;
index = index + 1;



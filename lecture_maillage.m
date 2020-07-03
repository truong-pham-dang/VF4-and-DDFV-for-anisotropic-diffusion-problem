% //**************************************************************************
% //   Fonction de lecture d'un maillage sur le disque dur
% //             Version MATLAB par DANG Truong
% //                     30/04/2018
% //**************************************************************************

function [maillage]=lecture_maillage(nom)
  

   fichier_sommets = importdata(strcat(nom,'_sommets'));
   [sommets]       = fichier_sommets.data;
   [text_som]      = fichier_sommets.textdata(1);
   
   fichier_centres = importdata(strcat(nom,'_centres'));
   [centres]       = fichier_centres.data;
   [text_cen]      = fichier_centres.textdata(1);
   
   fichier_aretes = importdata(strcat(nom,'_aretes'));
   [aretes]       = fichier_aretes.data;
   [text_aretes]  = fichier_aretes.textdata(1);
   
   maillage=struct('nom',text_aretes(1),...
                   'sommets',sommets,...
		  'centres',centres,...
		   'aretes',aretes,...
		   'nb_vol',size(centres,1),...
		   'nb_are',size(aretes,1),...
		   'nb_som',size(sommets,1));
end 
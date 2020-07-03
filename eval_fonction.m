% //**************************************************************************
% //   Fonction d'évaluation d'une fonction de deux variables sur une 
% //     famille de points
% //
% //    points : une matrice  N x 2
% //         f : une fonction (x,y) -> z
% //
% //   NB : Cette fonction n'est pas du tout optimisée !
% //**************************************************************************

function [v]=eval_fonction(points,f)
    v=zeros(size(points,1),1);
    for i=1:size(points,1)
       v(i)=f(points(i,1),points(i,2));
    end;
end
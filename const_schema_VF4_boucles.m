% //**************************************************************************
% //    Le schéma VF4 avec conditions de Dirichlet homogène
% //    programmé à l'aide de boucles
% //    Celui-ci n'est pas utilisé mais est seulement présenté
% //    pour des raisons pédagogiques
% //**************************************************************************

function [A,b]=const_schema_VF4_boucles(m,donnees)
  
  source=eval_fonction(centres,donnees.source);
  
  A=spzeros(m.nb_vol,m.nb_vol);
  b=zeros(m.nb_vol,1);
  
  
  for i=1:m.nb_are
    
    centre_arete(X)=(m.sommets(m.aretes(i,DEB),X)+m.sommets(m.aretes(i,FIN),X))/2;
    centre_arete(Y)=(m.sommets(m.aretes(i,DEB),Y)+m.sommets(m.aretes(i,FIN),Y))/2;
    coeff_diff=donnees.coeff_k(centre_arete(X),centre_arete(Y));
    tauKL = coeff_diff*m.aretes(i,MES)/m.aretes(i,DKL);
    
    A(m.aretes(i,K),m.aretes(i,K)) = A(m.aretes(i,K),m.aretes(i,K)) + tauKL;
    
    b(m.aretes(i,K)) = b(m.aretes(i,K)) ...
	+ (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_K_FIN))...
	   * donnees.source(m.centres(m.aretes(i,K),X) , m.centres(m.aretes(i,K),Y));
    
    
    if (m.aretes(i,L)>0) then
      A(m.aretes(i,K),m.aretes(i,L)) = - tauKL;
      A(m.aretes(i,L),m.aretes(i,L)) = A(m.aretes(i,L),m.aretes(i,L)) + tauKL;
      A(m.aretes(i,L),m.aretes(i,K)) = - tauKL;

      b(m.aretes(i,L)) = b(m.aretes(i,L)) ...
	+ (m.aretes(i,MES_L_DEB)+m.aretes(i,MES_L_FIN))...
	   * donnees.source(m.centres(m.aretes(i,L),X) , m.centres(m.aretes(i,L),Y));      
    end;

  end;

end
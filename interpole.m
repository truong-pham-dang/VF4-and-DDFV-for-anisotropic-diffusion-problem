% //**************************************************************************
% //    Routine d'interpolation de la fonction de courant (définie sur les sommets)
% //    sur une grille cartésienne.
% //     
% //    Ceci est nécessaire pour tracer les lignes de courant dans Scilab
% //
% //    Ne marche que sur les triangles
% //    et si le domaine est contenu dans le carré unité
% //**************************************************************************  

function [phi]=interpole(m,courant,Nx,Ny)

   max_courant=max(courant);
   min_courant=min(courant);
  
%    // On cherche les triangles
   
   tri=zeros(m.nb_vol,7);
   tri(:,1)=1;
   for i=1:m.nb_are
      j=m.aretes(i,K);
      tri(j,tri(j,1)+1)=m.aretes(i,DEB);
      tri(j,tri(j,1)+2)=m.aretes(i,FIN);
      tri(j,1)=tri(j,1)+2;
      
      j=m.aretes(i,L);
      if (j>0) 
	tri(j,tri(j,1)+1)=m.aretes(i,DEB);
	tri(j,tri(j,1)+2)=m.aretes(i,FIN);
	tri(j,1)=tri(j,1)+2;
      end;
   end;


   tri=tri(:,2:7);

   tri=gsort(tri,'c','i');
   tri=tri(:,[1 3 5]);


   x=linspace(0,1,Nx);
   y=linspace(0,1,Ny);
   [XX , YY] = meshgrid(x,y);
   phi=zeros(Nx,Ny);
   for i=1:size(XX,1)
      for j=1:size(XX,2)
	ptx=XX(i,j);
	pty=YY(i,j);
	coord2= ( (ptx-m.sommets(tri(:,1),X)).*(m.sommets(tri(:,3),Y)-m.sommets(tri(:,1),Y)) ...
                 -(pty-m.sommets(tri(:,1),Y)).*(m.sommets(tri(:,3),X)-m.sommets(tri(:,1),X)))...
              ./(  (m.sommets(tri(:,2),X)-m.sommets(tri(:,1),X)).*(m.sommets(tri(:,3),Y)-m.sommets(tri(:,1),Y)) ...
                 -(m.sommets(tri(:,2),Y)-m.sommets(tri(:,1),Y)).*(m.sommets(tri(:,3),X)-m.sommets(tri(:,1),X)));
	coord3= ( (ptx-m.sommets(tri(:,1),X)).*(m.sommets(tri(:,2),Y)-m.sommets(tri(:,1),Y)) ...
                 -(pty-m.sommets(tri(:,1),Y)).*(m.sommets(tri(:,2),X)-m.sommets(tri(:,1),X)))...
              ./(  (m.sommets(tri(:,3),X)-m.sommets(tri(:,1),X)).*(m.sommets(tri(:,2),Y)-m.sommets(tri(:,1),Y)) ...
                 -(m.sommets(tri(:,3),Y)-m.sommets(tri(:,1),Y)).*(m.sommets(tri(:,2),X)-m.sommets(tri(:,1),X)));
	coord1=1-coord2-coord3;

	[temp]=find(abs(coord1)+abs(coord2)+abs(coord3)-1<10^(-10),1);

	if (isempty(temp)) 
	   phi(i,j)=max_courant+100*(max_courant-min_courant);
	else
	phi(i,j)=coord1(temp)*courant(tri(temp,1))+coord2(temp)*courant(tri(temp,2)) ...
                +coord3(temp)*courant(tri(temp,3));
	end;
      end;
   end;
   
   phi=phi';

end
% //**************************************************************************
% //   Fonction de tracé d'une fonction discrète définie sur les centres
% //
% //   m : un maillage
% //   f : un vecteur colonne de taille m.nb_vol
% //   flag : determine exact solution 0, or approximate solution 1
% //**************************************************************************

function trace_fonction(m,f,flag)

X = 1; DEB = 1;
Y = 2; FIN = 2;
K = 3;
L = 4;

  
  fmin=min(f);
  fmax=max(f);
  

  
    s_X=zeros(3*m.nb_are,1);
    s_X(1:3:end)=m.centres(m.aretes(:,K),X);
    s_X(2:3:end)=m.sommets(m.aretes(:,DEB),X);
    s_X(3:3:end)=m.sommets(m.aretes(:,FIN),X);
    
    s_Y=zeros(3*m.nb_are,1);
    s_Y(1:3:end)=m.centres(m.aretes(:,K),Y);
    s_Y(2:3:end)=m.sommets(m.aretes(:,DEB),Y);
    s_Y(3:3:end)=m.sommets(m.aretes(:,FIN),Y);
    
    tri=zeros(m.nb_are,5);
    tri(:,2)=[1:3:3*m.nb_are]';
    tri(:,3)=[2:3:3*m.nb_are]';
    tri(:,4)=[3:3:3*m.nb_are]';
    
    val=zeros(3*m.nb_are,1);
    val(1:3:3*m.nb_are)=f(m.aretes(:,K));
    val(2:3:3*m.nb_are)=val(1:3:3*m.nb_are);
    val(3:3:3*m.nb_are)=val(1:3:3*m.nb_are);
    
    
    temp=find(m.aretes(:,L)>0);
    
    nb=length(temp);
    
    s_X2=zeros(3*nb,1);
    s_X2(1:3:end)=m.centres(m.aretes(temp,L),X);
    s_X2(2:3:end)=m.sommets(m.aretes(temp,DEB),X);
    s_X2(3:3:end)=m.sommets(m.aretes(temp,FIN),X);
    
    s_Y2=zeros(3*nb,1);
    s_Y2(1:3:end)=m.centres(m.aretes(temp,L),Y);
    s_Y2(2:3:end)=m.sommets(m.aretes(temp,DEB),Y);
    s_Y2(3:3:end)=m.sommets(m.aretes(temp,FIN),Y);
    
    tri2=zeros(nb,5);
    tri2(:,2)=[1:3:3*nb]';
    tri2(:,3)=[2:3:3*nb]';
    tri2(:,4)=[3:3:3*nb]';
    
    val2=zeros(3*nb,1);
    val2(1:3:3*nb)=f(m.aretes(temp,L));
    val2(2:3:3*nb)=val2(1:3:3*nb);
    val2(3:3:3*nb)=val2(1:3:3*nb);
    
    
    % First part of the solution
    switch flag
        case 0
            fid = fopen('exact_solution1.vtk','w');
        case 1
            fid = fopen('approximate_solution1.vtk','w'); 
        otherwise
            fprintf('Please enter 0 or 1. \n');
    end
    
    fprintf(fid,'# vtk DataFile Version 2.0 \n');
    fprintf(fid,'VTK format for unstructured mesh \n');
    fprintf(fid,'ASCII \n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
    fprintf(fid,'POINTS %d float \n',3*m.nb_are);
    for i = 1 : 3*m.nb_are
        fprintf(fid,'%f %f %f \n',s_X(i),s_Y(i),0);
    end
    fprintf(fid,'CELLS %d %d \n',m.nb_are,4*m.nb_are);
    for i = 1 : m.nb_are
        fprintf(fid,'%d %d %d %d \n',3,tri(i,2)-1,tri(i,3)-1,tri(i,4)-1);
    end
    fprintf(fid,'CELL_TYPES %d \n',m.nb_are);
    for i = 1 : m.nb_are
        fprintf(fid,'%d \n',5);
    end
    fprintf(fid,'POINT_DATA %d \n',3*m.nb_are);
    switch flag
        case 0
            fprintf(fid,'SCALARS u_exact double \n');
        case 1
            fprintf(fid,'SCALARS u_approximate double \n');
        otherwise
            fprintf('Please enter 0 or 1. \n');
    end
    fprintf(fid,'LOOKUP_TABLE default \n');
    for i = 1: 3*m.nb_are
        fprintf(fid,'%f \n',val(i));
    end
    fclose(fid);
    
    % Second part of the solution
    switch flag
        case 0
            fid = fopen('exact_solution2.vtk','w');
        case 1
            fid = fopen('approximate_solution2.vtk','w'); 
        otherwise
            fprintf('Please enter 0 or 1. \n');
    end
    
    fprintf(fid,'# vtk DataFile Version 2.0 \n');
    fprintf(fid,'VTK format for unstructured mesh \n');
    fprintf(fid,'ASCII \n');
    fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
    fprintf(fid,'POINTS %d float \n',3*nb);
    for i = 1 : 3*nb
        fprintf(fid,'%f %f %f \n',s_X2(i),s_Y2(i),0);
    end
    fprintf(fid,'CELLS %d %d \n',nb,4*nb);
    for i = 1 : nb
        fprintf(fid,'%d %d %d %d \n',3,tri2(i,2)-1,tri2(i,3)-1,tri2(i,4)-1);
    end
    fprintf(fid,'CELL_TYPES %d \n',nb);
    for i = 1 : nb
        fprintf(fid,'%d \n',5);
    end
    fprintf(fid,'POINT_DATA %d \n',3*nb);
    switch flag
        case 0
            fprintf(fid,'SCALARS u_exact double \n');
        case 1
            fprintf(fid,'SCALARS u_approximate double \n');
        otherwise
            fprintf('Please enter 0 or 1. \n');
    end
    fprintf(fid,'LOOKUP_TABLE default \n');
    for i = 1: 3*nb
        fprintf(fid,'%f \n',val2(i));
    end
    fclose(fid);
    

end
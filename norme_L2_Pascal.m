% Author:  DANG Truong
% Date:    20/05/2018
% Purpose: This function evaluates the L2 norm mentioned in 
%          section "8. Numerical results" of the paper:
%           A FINITE VOLUME METHOD FOR THE LAPLACE EQUATION ON ALMOST
%                       ARBITRARY TWO-DIMENSIONAL GRIDS
%                                   -----
%                       Komla Domelevo and Pascal Omnes

function [norme] = norme_L2_Pascal(m,primal_subtraction,dual_subtraction,primal_solexact,dual_solexact)

global K;
global L;
global DEB;
global MES_K_DEB;
global MES_K_FIN;
global MES_L_DEB;
global MES_L_FIN;

    % Initialize norme
    norme = 0;
    
    denominator = 0;
    
    % Calculate the first sums of the numerator and the denominator 
    % (related to cell volumes of the primal mesh)
    
    for i = 1:m.nb_are
        norme = norme + (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_K_FIN))* ...
	    primal_subtraction(m.aretes(i,K))^2;
    
        denominator = denominator + (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_K_FIN))* ...
	    primal_solexact(m.aretes(i,K))^2;
        
        % Not a boundary edge
        if (m.aretes(i,L)>0)
            norme = norme + (m.aretes(i,MES_L_DEB)+m.aretes(i,MES_L_FIN))* ...
            primal_subtraction(m.aretes(i,L))^2;
        
            denominator = denominator + (m.aretes(i,MES_L_DEB)+m.aretes(i,MES_L_FIN))* ...
            primal_solexact(m.aretes(i,L))^2;
        end
    end
    
    % Calculate the second sums in the numerator and the denominator 
    % (related to cell volumes of the dual mesh)
    
    for i = 1:m.nb_are
        % Not a boundary edge
        if (m.aretes(i,L)>0)
            norme = norme + (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_L_DEB))*...
            dual_subtraction(m.aretes(i,DEB))^2;
        
            denominator = denominator + (m.aretes(i,MES_K_DEB)+m.aretes(i,MES_L_DEB))*...
            dual_solexact(m.aretes(i,DEB))^2;
        % A boundary edge
        else 
            norme = norme + m.aretes(i,MES_K_DEB)*dual_subtraction(m.aretes(i,DEB))^2;
            
            denominator = denominator + m.aretes(i,MES_K_DEB)*dual_solexact(m.aretes(i,DEB))^2;
        end
    end
    
    norme       = 0.5 * norme;
    denominator = 0.5 * denominator;
    
    norme       = norme / denominator;
    norme       = sqrt(norme);
    
end 
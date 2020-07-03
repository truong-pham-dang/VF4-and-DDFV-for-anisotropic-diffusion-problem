function [z]=f6(x,y)
 const=1/(mat_xx_dis(0,0)-mat_xy_dis(0,0));
  if (x>0.5) 
    z=1+const^2;
  else
    z=(mat_xx_dis(0,0)-2*mat_xy_dis(0,0)+mat_yy_dis(0,0))*const^2;
  end;
end
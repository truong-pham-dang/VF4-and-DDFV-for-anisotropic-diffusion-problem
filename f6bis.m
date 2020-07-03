function [z]=f6bis(x,y)
 const=1/(coeff_dis(0,0));
  if (x>0.5) 
    z=1+const^2;
  else
    z=(2*coeff_dis(0,0))*const^2;
  end;
end
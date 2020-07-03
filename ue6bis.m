function [z]=ue6bis(x,y)
  const=1/(coeff_dis(0,0));
  if (x>0.5) 
    z=-0.5*(x-const*y+(const-1)*0.5)^2;
  else
    z=-0.5*(const*(x-y))^2;
  end;
end
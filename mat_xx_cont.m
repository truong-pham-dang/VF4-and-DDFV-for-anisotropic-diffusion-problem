function [z]=mat_xx_cont(x,y)
  z=(10^(-3)*x.^2+y.^2)./(x.^2+y.^2);
end
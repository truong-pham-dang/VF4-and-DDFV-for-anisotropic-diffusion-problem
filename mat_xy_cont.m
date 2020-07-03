function [z]=mat_xy_cont(x,y)
  z=(10^(-3)-1)*x.*y./(x.^2+y.^2);
end
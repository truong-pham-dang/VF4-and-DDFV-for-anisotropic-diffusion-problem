function [z]=mat_yy_cont(x,y)
  z=(x.^2+10^(-3)*y.^2)./(x.^2+y.^2);
end
function result=function_f(t,x)
global delta;

result=(delta^2/2 + 6*x.^2 - 1) .* sin(t)+(x.^2 - x.^4) .* cos(t);

end
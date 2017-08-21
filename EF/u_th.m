function [u_th ] = u_th(x,y,mu,a,l)
 
 u_th=-(a*sqrt(3)/(2*mu*l))*y*(x-y/sqrt(3)+l/2)*(x+y/sqrt(3)-l/2);

end


function [cl, cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)
%forces calculate the lift and drag forces on the aerofoil

cl = -2*circ;
theta_te = thetal(end)+thetau(end);
H = (delstarl(end)+delstaru(end))/(theta_te);
theta_inf = theta_te*sqrt(1-cp(end))^((H+5)/2);
cd = 2*theta_inf;

end


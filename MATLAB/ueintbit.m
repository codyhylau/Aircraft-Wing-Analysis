function f = ueintbit(xa,ua,xb,ub)
%ueintbit Calculates the contribution to the integral between xa and xb
% x'/L=xa(where ue/U = ua) and b x'/L=xb(where ue/U = ub)

ubar = (ua+ub)/2;
delu = ub - ua;
delx = xb - xa;


f = (ubar^5 + (5/6)*(ubar^3)*(delu^2) + (1/16)*ubar*(delu^4))*delx; %handout p.13

end


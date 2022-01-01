function [infa,infb] = panelinf(xa,ya,xb,yb,x,y)
%panelinf returns the influence coefficients fa, fb at (x,y) due to our general panel
%   Detailed explanation goes here

T=[xb-xa;yb-ya]; %vector of panel
del = norm(T); %length of panel
t = T/del; %normalised tangential vector of panel
n = [0,-1;1,0]*t; %normal vector

r = [x-xa;y-ya]; %vector between the generalised point and the (xa,ya)
X = dot(r,t); %component of r tangential to the panel
Y = dot(r,n); %component of r normal to the panel
[infa,infb] = refpaninf(del,X,Y); 

end


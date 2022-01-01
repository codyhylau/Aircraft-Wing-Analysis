function dthickdx = thickdash(xmx0,thick)
%thickdash is f(x,y), and thick is [theta, deltaE], xmx0 is the x variable

global Re ue0 duedx

He = thick(2)/thick(1); %He = deltaE/theta
if He >= 1.46 %equation 17a
    H = (11*He+15)/(48*He-59);
else %equation 17b
    H = 2.803; 
end
% ue = ue0 + duedx*xmx0; %linear varying velocity gradient
Retheta = Re*ue0*thick(1); 

cf = 0.091416*(((H-1)*Retheta)^(-0.232))*exp(-1.26*H); %equation 18
cdiss = 0.010011*((H-1)*Retheta)^(-1/6); %equation 19

dthickdx1= (cf/2) - ((H+2)/ue0)*duedx*thick(1); %momentum integral equation
dthickdx2 = cdiss - (3/ue0)*duedx*thick(2); %KE integral equation
dthickdx = [dthickdx1;dthickdx2]; %'y'
end


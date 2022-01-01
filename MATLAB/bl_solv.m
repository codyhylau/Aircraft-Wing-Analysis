function [int, ils, itr, its, delstar, theta] =  bl_solv(x,cp)
%bl_solv Summary of this function goes here

np = length(x);

global ue0 Re duedx


ue0 = 1;

theta = zeros(1,np); %Theta/L
delE = theta; %Energy thickness
delstar = theta; %displacement thickness

f = 0;
laminar = true;
int = 0; %natural transition
ils = 0; %Laminar separation
itr = 0; %Turbulent reattachment
its = 0; %Turbulent separation
H = 0;
He = zeros(1,np);
He(1) = 1.57258; %Blasius BL for x/L = 0 point

ue = sqrt(1-cp); % Array of velocity ue/U
i = 1;
duedx = (ue(i))/(x(i));
f = f + ueintbit(0,0,x(i),ue(i));
theta(i) = sqrt((0.45/Re)*(ue(i)^(-6))*f); %theta(x)/L (Thwaites method)

Retheta = Re*ue(i)*theta(i);
m = -Re*(theta(i)^2)*duedx;
H = thwaites_lookup(m); %moodle functions
He(i) = laminar_He(H);
delE(i) = He(i) * theta(i); %Energy thickness
delstar(i) = H * theta(i); %displacement thickness

if log(Retheta) >= (18.4*He(i) - 21.74) %Check for natural transition to turbulent
    laminar = false;  
    int = i;
elseif m >= 0.09 %Check for laminar separation
    laminar = false;
    ils = i;
    He(i) = 1.51509; %laminar separation value
    delE(i) = He(i) * theta(i);
end 

while laminar && i<np
    i = i+1;
    duedx = (ue(i-1)-ue(i))/(x(i-1)-x(i));
    f = f + ueintbit(x(i-1),ue(i-1),x(i),ue(i));
    theta(i) = sqrt((0.45/Re)*(ue(i)^(-6))*f); %theta(x)/L (Thwaites method)
    
    Retheta = Re*ue(i)*theta(i);
    m = -Re*(theta(i)^2)*duedx;
    H = thwaites_lookup(m); %moodle functions
    He(i) = laminar_He(H);
    delE(i) = He(i) * theta(i); %Energy thickness
    delstar(i) = H * theta(i); %displacement thickness
    
    if log(Retheta) >= (18.4*He(i) - 21.74) %Check for natural transition to turbulent
        laminar = false;  
        int = i;
    elseif m >= 0.09 %Check for laminar separation
        laminar = false;
        ils = i;
        He(i) = 1.51509; %laminar separation value
        delE(i) = He(i) * theta(i);
    end  
end

while its == 0 && i < np
    i = i+1;
    thick0 = [theta(i-1);delE(i-1)]; %update 'initial' conditions for each panel
    ue0 = ue(i-1); %update initial free stream velocity for each panel
    duedx = (ue(i-1)-ue(i))/(x(i-1)-x(i)); %update initial velocity gradient for each panel
    [delx, thickhist] = ode45(@thickdash,[0,x(i)-x(i-1)],thick0); %solve ode
    theta(i) = thickhist(end,1); %values at the end of each panel
    delE(i) = thickhist(end,2);
    He(i) = delE(i)/theta(i); 
    H = (11*He(i)+15)/(48*He(i)-59);
        
    if ils ~= 0 && itr == 0 && He(i)>1.58 %Check for reattachment
        itr = i;
    end
    if He(i)<1.46 %Check for turbulent separation
        its = i;
        H = 2.803;
    end
    delstar(i) = H*theta(i);
end

if its ~= 0 %If turbulent separation has occured obtain theta on the separated region

    H = 2.803; %Equation 17b from handout, H unchanged after separation
    for j = its+1:np
        theta(j) = theta(j-1)*(ue(j-1)/ue(j))^(H+2);
    end

    delstar(its+1:np) = H * theta(its+1:np);
end

end


global ue0 Re grad

Re = 1e7; % Re_L, 10^6/10^7 for parts a&b, 10^4/10^5/10^6 for parts c&d
grad = -0; % d(ue/U)/d(x/L), 0 for parts a&b, -0.25 for parts c&d
% critical velocity gradient = -0.51 (turbulent separation at x/L = 1)
ue0 = 1;
disp('Critical velocity gradient for Re=10^5: -0.38 (2.s.f.)')

np = 101;
x = linspace(0,1,np); %x/L
ue = linspace(ue0,grad+ue0,np); %ue/U

theta = zeros(1,np); %Theta/L

f = 0;
laminar = true;
int = 0; %natural transition
ils = 0; %Laminar separation
itr = 0; %Turbulent reattachment
its = 0; %Turbulent separation
He = zeros(1,np);
He(1) = 1.57258; %Blasius BL for x/L = 0 point

i = 1;
while laminar && i<np
    i = i+1;
    f = f + ueintbit(x(i-1),ue(i-1),x(i),ue(i));
    theta(i) = sqrt((0.45/Re)*(ue(i)^(-6))*f); %theta(x)/L (Thwaites method)
    
    Retheta = Re*ue(i)*theta(i);
    m = -Re*(theta(i)^2)*grad;
    H = thwaites_lookup(m); %moodle functions
    He(i) = laminar_He(H);
    
    if log(Retheta) >= (18.4*He(i) - 21.74) %Check for natural transition to turbulent
        laminar = false;  
        int = i;
    elseif m >= 0.09 %Check for laminar separation
        laminar = false;
        ils = i;
        He(i) = 1.51509; %laminar separation value
    end  
end

deltae = He(i)*theta(i); % temporary deltaE as only needed to set first thick0
while its == 0 && i < np
    i = i+1;
    thick0 = [theta(i-1);deltae]; %update 'initial' conditions for each panel
    ue0 = ue(i-1); %update initial free stream velocity for each panel
    [delx, thickhist] = ode45(@thickdash,[0,x(i)-x(i-1)],thick0); %solve ode
    theta(i) = thickhist(end,1); %values at the end of each panel
    deltae = thickhist(end,2);
    He(i) = deltae/theta(i); 
        
    if ils ~= 0 && itr == 0 && He(i)>1.58 %Check for reattachment
        itr = i;
    end
    if He(i)<1.46 %Check for turbulent separation
        its = i;
    end
end

if its ~= 0 %If turbulent separation has occured obtain theta on the separated region
    He(its:np) = 1.46; %Value at separation point
    H = 2.803; %Equation 17b from handout, H unchanged after separation
    R = (ue(1)/ue(2))^(H+2); %constant
    for i = its+1:np
        theta(i) = theta(i-1)*R;
    end
end
    
disp(['[ils,int,itr,its] = ' num2str([ils,int,itr,its])])


% PLOTTING METHOD Re-run script but change Re and/or gradient with previous...
% theta or He saved in Matlab Workspace (e.g. plot He1, He2, He3)

theta2 = theta;
He2 = He;

plot(x,theta1, 'linewidth', 2) % plot with parameters saved in workspace
hold on
plot(x,theta2, 'linewidth', 2)
% hold on
% plot(x,He3, 'linewidth', 2)

xlabel("X")
% ylabel("He") %change to theta/L for theta plots
ylabel("\theta/L")
legend("ReL = 10^6", "ReL = 10^7",'Location','northwest')
axis auto
% set figure position and size:
set(gcf,'position',[160 280 600 460])
% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')
%set fonts and frame:
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
title({'Momentum thickness of a boundary','layer for a zero pressure gradient'}) 
% title({'Energy shape factor of a boundary','layer for a zero pressure gradient'}) 
% ax = gca;
% ax.YAxis.Exponent = 0;
% print -dpng -r300 ex6a.png

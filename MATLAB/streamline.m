%    foil.m
%
%  Script to analyse an aerofoil section using potential flow calculation
%  and boundary layer solver.
%

clear all

global Re

%  Read in the parameter file
% caseref = input('Enter case reference: ','s');
% parfile = ['Parfiles/' caseref '.txt'];
% fprintf(1, '%s\n\n', ['Reading in parameter file: ' parfile])
% [section np Re alpha_array] = par_read(parfile);

section = 'naca63-412-9';
np = 500;
Re = 0.5e6;

%  Read in the section geometry
secfile = ['Geometry/' section '.surf'];
[xk yk] = textread ( secfile, '%f%f' );

%  Generate high-resolution surface description via cubic splines
nphr = 5*np;
[xshr yshr] = splinefit ( xk, yk, nphr );

%  Resize section so that it lies between (0,0) and (1,0)
[xsin ysin] = resyze ( xshr, yshr );

%  Interpolate to required number of panels (uniform size)
[xs ys] = make_upanels ( xsin, ysin, np );

%  Assemble the lhs of the equations for the potential flow calculation
A = build_lhs ( xs, ys );
Am1 = inv(A);

alpha = input('Angle of attack to plot: ');

%    rhs of equations
  alfrad = pi * alpha/180;
  b = build_rhs ( xs, ys, alfrad );

%    solve for surface vortex sheet strength
  gam = Am1 * b;

%    calculate cp distribution and overall circulation
  [cp circ] = potential_op ( xs, ys, gam );

%    locate stagnation point and calculate stagnation panel length
  [ipstag fracstag] = find_stag(gam);
  dsstag = sqrt((xs(ipstag+1)-xs(ipstag))^2 + (ys(ipstag+1)-ys(ipstag))^2);

%    upper surface boundary layer calc

%    first assemble pressure distribution along bl
  clear su cpu
  su(1) = fracstag*dsstag;
  cpu(1) = cp(ipstag);
  for is = ipstag-1:-1:1
    iu = ipstag - is + 1;
    su(iu) = su(iu-1) + sqrt((xs(is+1)-xs(is))^2 + (ys(is+1)-ys(is))^2);
    cpu(iu) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag < 1e-6
    su(1) = 0.01*su(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpu(2));
    cpu(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [iunt iuls iutr iuts delstaru thetau] = bl_solv ( su, cpu );

%    lower surface boundary layer calc

%    first assemble pressure distribution along bl
  clear sl cpl
  sl(1) = (1-fracstag) * dsstag;
  cpl(1) = cp(ipstag+1);
  for is = ipstag+2:np+1
    il = is - ipstag;
    sl(il) = sl(il-1) + sqrt((xs(is-1)-xs(is))^2 + (ys(is-1)-ys(is))^2);
    cpl(il) = cp(is);
  end

%    check for stagnation point at end of stagnation panel
  if fracstag > 0.999999
    sl(1) = 0.01*sl(2);    % go just downstream of stagnation
    uejds = 0.01 * sqrt(1-cpl(2));
    cpl(1) = 1 - uejds^2;
  end

%    boundary layer solver
  [ilnt, ills, iltr, ilts, delstarl, thetal] = bl_solv ( sl, cpl );

%    lift and drag coefficients
  [Cl, Cd] = forces ( circ, cp, delstarl, thetal, delstaru, thetau );


  upperbl = sprintf ( '%s', '  Upper surface boundary layer:' );
  if iunt~=0
    is = ipstag + 1 - iunt;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Natural transition at x = ', xs(is) );
  end
  if iuls~=0
    is = ipstag + 1 - iuls;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Laminar separation at x = ', xs(is) );
    if iutr~=0
      is = ipstag + 1 - iutr;
      upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
    end
  end
  if iuts~=0
    is = ipstag + 1 - iuts;
    upperbl = sprintf ( '%s\n%s%5.3f', upperbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
  end
  upperbl = sprintf ( '%s\n', upperbl );
  disp(upperbl)

  lowerbl = sprintf ( '%s', '  Lower surface boundary layer:' );
  if ilnt~=0
    is = ipstag + ilnt;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Natural transition at x = ', xs(is) );
  end
  if ills~=0
    is = ipstag + ills;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Laminar separation at x = ', xs(is) );
    if iltr~=0
      is = ipstag + iltr;
      lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                          '    Turbulent reattachment at x = ', xs(is) );
    end
  end
  if ilts~=0
    is = ipstag + ilts;
    lowerbl = sprintf ( '%s\n%s%5.3f', lowerbl, ... 
                        '    Turbulent separation at x = ', xs(is) );
  end
  lowerbl = sprintf ( '%s\n', lowerbl );
  disp(lowerbl)


%NOW DO STREAMLINE CODE

nx = 300;
ny = 250; %Number of points in space
xmin = -0.5;
xmax = 1.5; 
ymin = -0.5;
ymax = 0.5; %Region to measure psi at

% np = 100; %number of panels around cylinder
% theta = (0:np)*2*pi/np; % array of angles at each panel
% xs = zeros(1,np); % pre-allocation for panel coordinates
% ys = xs;
% gamma = xs;

% for i = 1:np+1
%     xs(i) = cos(theta(i)); % coordinate for each panel
%     ys(i) = sin(theta(i));
%     gamma(i) = -2*sin(theta(i)); % panel vortex strength 
% end

xm = zeros(nx, ny); % pre-allocation
ym = xm;
psi_fs = xm;
psi_cy = xm;
psi = xm;
fa = xm;
fb = xm;

for i=1:nx
    for j=1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1); %Loop through all the points obtaining xm to get r in our psiv function
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1); %Loop through all the points obtaining ym to get r in our psiv function
        psi_fs(i,j) = ym(i,j)*cos(alfrad)-xm(i,j)*sin(alfrad); % free-stream contribition
        for k = 1:length(gam)-1
            [infa,infb] = panelinf(xs(k),ys(k),xs(k+1),ys(k+1),xm(i,j),ym(i,j)); %Uses the function refpaninf to obtaint the influence coefficients
            psi_cy(i,j) = psi_cy(i,j)+ gam(k)*infa + gam(k+1)*infb;  % summing all cylinder streamfunction contributions using influence coefficients 
%               psi_cy(i,j) = psi_cy(i,j) + A(k,:)*gam;
        end
        psi(i,j) = psi(i,j)+ psi_fs(i,j) + psi_cy(i,j); % summing free-stream and cylinder streamfunctions
    end
end

c = min(min(psi)):0.02:max(max(psi));
figure(2)
contour(xm,ym,psi,c,'black')
% set figure position and size:
set(gcf,'position',[160 100 900 690])
% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')
%set fonts and frame:
set(gca,'Fontn','Times','FontSize',18,'linewidth',1)
xlabel("X")
ylabel("Y")

hold on
plot(xs,ys,'linewidth', 2) % plot cylinder shape
axis([-0.2,1.2,-0.5,0.5])
legend("Streamlines", "aerofoil")
% print -dpng -r300 6212_1.png
hold off

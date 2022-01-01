%Plotting

foils = ["naca63-412-9_low"];
foil_names = ["2nd iteration","3rd Iteration","final iteration"];

% figure(2)
% box on
% % for i = 1:length(foils)
% %     hold on
% %     load(strcat("Data/", foils(i),"/",foils(i),"_7.6.mat"))
% %     plot(xs,-cp)
% % end
% 
% hold on
% load(strcat("Data/", foils(1),"/",foils(1),"_5.mat"))
% plot(xs,-cp,'linewidth',1)
% load(strcat("Data/", foils(2),"/",foils(2),"_7.6.mat"))
% plot(xs,-cp,'linewidth',1)
% load(strcat("Data/", foils(3),"/",foils(3),"_7.2.mat"))
% plot(xs,-cp,'linewidth',1)
% hold off
% xlabel("x/c")
% ylabel("-C_P")
% legend([strcat("2nd Iteration @ \alpha = 5",char(176)), strcat("3rd Iteration @ \alpha = 7.6",char(176)),strcat("final Iteration @ \alpha = 7.2",char(176)) ])
% % x-axis ticks:
% % set(gca,'xtick',[0 1 2])
% % y-axis ticks:
% % set(gca,'ytick',[-1 0 1])
% 
% % set figure position and size:
% set(gcf,'position',[160 280 600 460])
% % keep position and size when printing:
% set(gcf,'PaperPositionMode','auto')
% %set fonts and frame:
% set(gca,'Fontn','Times','FontSize',16,'linewidth',1)

figure(2)
plot(alpha, clswp, 'linewidth',1)
xlabel("\alpha")
ylabel("C_L")
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';

% set figure position and size:
set(gcf,'position',[160 280 600 460])
% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')
%set fonts and frame:
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
xticks(-20:2:20)
yticks(-1.6:0.4:2.0)


% 
% % figure(3)
% % plot(alpha,cdswp)
% % xlabel("\alpha")
% % ylabel("C_D")
% % ax = gca;
% % ax.XAxisLocation = 'origin';
% % ax.YAxisLocation = 'origin';
% 
% % set figure position and size:
% set(gcf,'position',[160 280 600 460])
% % keep position and size when printing:
% set(gcf,'PaperPositionMode','auto')
% %set fonts and frame:
% set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
% 
% 


figure(4)

for i = 1:length(foils)
    hold on
    load(strcat("Data/", foils(i)))
    plot(clswp,cdswp, 'linewidth',1)  
end    
box on
xlabel("C_L")
ylabel("C_D")
legend(["naca6212","1st Iteration"], 'location','southeast')

% set figure position and size:
set(gcf,'position',[160 280 600 460])
% keep position and size when printing:
set(gcf,'PaperPositionMode','auto')
%set fonts and frame:
set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
set(gca,'xtick',-1.6:0.4:2.5)
set(gca,'ytick',0:0.004:0.029)
axis([-0.4,2.4,0,0.028])
hold off




% 
% 
% figure(5)
% plot(xs,ys,'black')
% xlabel("x/c")
% ylabel("y/c")
% axis([-0.1,1.1,-0.1,0.2])
% yticks(-0.1:0.1:0.3)
% 
% % set figure position and size:
% set(gcf,'position',[160 280 600 200])
% % keep position and size when printing:
% set(gcf,'PaperPositionMode','auto')
% %set fonts and frame:
% set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
% print(figure(5),'3rd_section','-dpng', '-r300')


% figure(6)
% % foils = ["naca0012" "naca6112" "naca6212" "naca6312"];
% % foils = ["naca6212","naca6212_1"];
% for i = 1:length(foils)
%     hold on
%     load(strcat("Data/", foils(i)))
%     plot(alpha,lovdswp, 'linewidth',1)  
% end    
% box on
% xlabel("\alpha")
% ylabel("L/D")
% % legend(foil_names, 'location','southeast')
% legend(["Re = 0.5e6","Re = 10e6","Re = 20e6"], 'location','southeast')
% set(gca,'xtick',-6:2:12)
% set(gca,'ytick',0:20:210)
% 
% % set figure position and size:
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% set(gcf,'position',[160 280 600 460])
% % keep position and size when printing:
% set(gcf,'PaperPositionMode','auto')
% %set fonts and frame:
% set(gca,'Fontn','Times','FontSize',16,'linewidth',1)
% axis([-6,12, 0, 220])
% hold off


%Create folder and save graphs
% mkdir(section)
% print(figure(1),[section,'/C_P'],'-dpng', '-r300') %FIX FOR PARTICULAR ANGLE, NOW PLOTS LAST ANGLE
% print(figure(2),[section,'/C_L'],'-dpng', '-r300')
% print(figure(3),[section,'/C_D'],'-dpng', '-r300')
% print(figure(4),['Report','/CLCD_initial'],'-dpng', '-r300')
% print(figure(5),[section,'/section'],'-dpng', '-r300')
% print(figure(6),['Report','/LovrD_initial'],'-dpng', '-r300')
% ADD L/D PLOT
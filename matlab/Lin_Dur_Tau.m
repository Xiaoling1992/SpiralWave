clear all;
data=load('/home/xiaoling/SpiralWave/data/F17_T_Tau_sta.dat');
WF=data(1,:);
duration=data(end,:);

taus=[10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];


figure;
plot(taus( taus>100), duration(taus>100), '-O', 'linewidth',3,'markersize' ,7,'MarkerFacecolor','c', 'Color','c');
% hold on
% plot(phis(phis>0.3 ), anchor(1,phis>0.3 )*pi/sqrt(epsilon/2), '-O', 'linewidth',3,'markersize' ,7, 'MarkerFacecolor','m','Color','m');
% legend(gca, {'dynamical (data)','$$\pi d_{struct}/ \sqrt{\frac{\epsilon}{2} }$$'}, 'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex');
% 
% plot([0.41, 0.41], [0, 120], '--', 'LineWidth', 2, 'Color', 'g');

%legend(gca, {'N=100','N=300','N=1000'},  'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex' )

%xlim([2*10^(-2), 2.2*10^(-1)]);
%ylim([2, 50]);
xlabel('\tau', 'fontsize',20)
ylabel('duration', 'fontsize', 20)
%set(gca, 'xticks', 0.41,'xticklabel','$$p_c$$', 'fontsize', 15, 'FontWeight', 'bold','interpreter','latex' )
%set(gca,'xtick',[0.3, 0.35, 0.41, 0.45, 0.5, 0.55, 0.6],'xticklabel',{0.3, 0.35, 'p_c', 0.45, 0.5, 0.55, 0.6}, 'fontsize', 15);

%text(0.48, 42, 'dynamical (data)', 'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.40, 80, '$$\pi d_{struct}/ \sqrt{\frac{\epsilon}{2} }$$',  'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.35, 60, '$$\downarrow$$',  'Color','r','FontSize',70, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.5, 25, '$$\uparrow$$', 'Color','r', 'FontSize',40, 'FontWeight', 'bold','interpreter', 'latex')

print(gcf,['../PapFig/Lin_Dur_Tau.png'],'-dpng','-r400');

dlmwrite(['../PapFig/Duration_Tau.dat'], 'parameters are the same as WT except tau. v_0=0.02, tau_off=5, epsilon=10, phi=0.43, rho=0', 'delimiter', '');
dlmwrite(['../PapFig/Duration_Tau.dat'], 'tau','-append', 'delimiter', '');
dlmwrite(['../PapFig/Duration_Tau.dat'],taus(taus>=100), '-append', 'delimiter', '\t', 'precision', 4);
dlmwrite(['../PapFig/Duration_Tau.dat'], 'duration','-append', 'delimiter', '');
dlmwrite(['../PapFig/Duration_Tau.dat'],duration(taus>=100), '-append', 'delimiter', '\t', 'precision', 4);

clear all;
Ns=[100, 300, 1000];
phis=[0.430, 0.434, 0.438, 0.443, 0.449, 0.456, 0.464, 0.474, 0.485, 0.497, 0.512, 0.529, 0.549, 0.573, 0.600];
ColMap=[
        255,0,255;
        128,0,128;
        255,165,0]/255;

   figure 
for i= 1:length(Ns)
    N=Ns(i);
    data= load(['/home/xiaoling/SpiralWave/data/Spi7_MClu_Phi_N',num2str(N),'.dat']);
    MeanRadius(:, i)= mean( data( length(phis)+1: 2*length(phis), :), 2);
 
    loglog(phis- 0.407, MeanRadius(:,i), '-', 'linewidth',3,'Color',ColMap(i,:))
    hold on
   % loglog(phis- 0.407, MeanRadius, '-O', 'linewidth',3,'markersize' ,7, 'Color',ColMap(i,:))
  
end
legend(gca, {'\ N=100','\ N=300','\ N=1000'}, 'location', 'SouthWest', 'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex' )
legend boxoff


loglog(phis(6:1:9)- 0.407,( (phis(6:1:9)-0.407).^(-43*2/36) )/4,'-',  'linewidth',3,'Color', 'k' )
xlim([0.01,0.22]);
%ylim([1, 50])

set(gca, 'xtick',[0.01, 0.1],'xticklabels',{'10^{-2}','10^{-1}'},'fontsize',15,'ytick',[10, 100, 1000],'yticklabels',{'10^{1}','10^{2}', '10^{3}'},'fontsize',15);
xlabel('$$\phi-p_c$$', 'fontsize',25,'FontWeight', 'bold','interpreter','latex')
ylabel('$$n_{off}$$', 'fontsize', 25,  'FontWeight', 'bold', 'interpreter','latex')

text(0.055,300, '$$-\gamma$$', 'FontSize',25, 'FontWeight', 'bold','interpreter', 'latex')

%print(gcf,['../PapFig/Loglog_MClu_Phi_Pc_1.png'],'-dpng','-r400');

figure();

plot(1,1, 'ks', 'linewidth',1,'markersize' ,10,'MarkerFacecolor','b','MarkerEdgecolor','b');
hold on
plot(2,2, 'ks', 'linewidth',1,'markersize' ,10,'MarkerFacecolor','g','MarkerEdgecolor','g');
plot(2,2, 'ks', 'linewidth',1,'markersize' ,10,'MarkerFacecolor','r','MarkerEdgecolor','r');
%legend(gca,{'firing cell ($$\tau_{on}$$)','non-firing cell ($$\tau_{off}=5$$)'}, 'FontSize',20, 'interpreter','latex');
legend(gca,{'wave death','radial propagation','spiral'}, 'FontSize',15, 'interpreter','latex');
legend boxoff 
print(gcf,['../PapFig/Cartoon_LegOnOff_3.png'],'-dpng','-r500');

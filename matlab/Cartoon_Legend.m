figure();

plot(1,1, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','w','MarkerEdgecolor',[0.5, 0.5, 0.5]);
hold on
plot(2,2, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','k','MarkerEdgecolor', [0.5, 0.5, 0.5]);

%legend(gca,{'firing cell ($$\tau_{on}$$)','non-firing cell ($$\tau_{off}=5$$)'}, 'FontSize',20, 'interpreter','latex');
legend(gca,{' on',' off'}, 'FontSize',15);
legend boxoff 
print(gcf,['../PapFig/Cartoon_LegOnOff.png'],'-dpng','-r400');

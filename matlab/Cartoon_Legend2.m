figure();

plot(1,1, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','w','MarkerEdgecolor',[0.5,0.5, 0.5]);
hold on
plot(2,2, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','k','MarkerEdgecolor',[0.5, 0.5, 0.5]);
leg1={ [' on' char(10) '($$\tau_{on}$$)'], [' off' char(10) '($$\tau_{off}=5$$)']};
leg2={ ['      ' char(10) '      '], ['      ' char(10) '      ']};

%legend(gca,{'firing cell ($$\tau_{on}$$)','non-firing cell ($$\tau_{off}=5$$)'}, 'FontSize',20, 'interpreter','latex');
%legend(gca,{' on($$\tau_{on}$$)', ' off($$\tau_{off}=5$$)'}, 'FontSize',15, 'interpreter','latex');

legend(gca,leg2, 'FontSize',15, 'interpreter','latex');
legend boxoff 
print(gcf,['../PapFig/Cartoon_LegOnOff_2.png'],'-dpng','-r500');

figure();

ColMap=[
       0 0 1
       0 0.5 0.5
       255/255, 192/255, 203/255  %pink 
       1 0 0
       ];
plot(7,9, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','b','MarkerEdgecolor','b');
hold on
plot(7,7.8, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor', ColMap(2, :),'MarkerEdgecolor',ColMap(2, :));
plot(7,6.6, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','r','MarkerEdgecolor','r');
plot(7,5.4, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor',ColMap(3, :),'MarkerEdgecolor',ColMap(3, :));
xlim([1, 10])
ylim([1,10])

text(7.3,9.1,{'wave' } , 'FontSize',15, 'interpreter','latex');
text(7.3,8.65,{'death'} , 'FontSize',15, 'interpreter','latex');

text(7.3,7.9,{'radial'} , 'FontSize',15, 'interpreter','latex');
text(7.3,7.45,{'prop.'} , 'FontSize',15, 'interpreter','latex');

text(7.3,6.7,{'secondary'} , 'FontSize',15, 'interpreter','latex');
text(7.3,6.25,{'spiral'} , 'FontSize',15, 'interpreter','latex');

text(7.3,5.5,{'break-up'} , 'FontSize',15, 'interpreter','latex');
text(7.3,5.05,{'spiral'} , 'FontSize',15, 'interpreter','latex');

% leg1={ [' on' char(10) '($$\tau_{on}$$)'], [' off' char(10) '($$\tau_{off}=5$$)']};
% leg2={ ['      ' char(10) '      '], ['      ' char(10) '      ']};

%legend(gca,{'firing cell ($$\tau_{on}$$)','non-firing cell ($$\tau_{off}=5$$)'}, 'FontSize',20, 'interpreter','latex');
%legend(gca,{' on($$\tau_{on}$$)', ' off($$\tau_{off}=5$$)'}, 'FontSize',15, 'interpreter','latex');

%legend(gca,leg2, 'FontSize',15, 'interpreter','latex');
%legend boxoff 
print(gcf,['../PapFig/Cartoon_LegOnOff_5.png'],'-dpng','-r500');

figure();

plot(7,9, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','b','MarkerEdgecolor','b');
hold on
plot(7,7.8, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','g','MarkerEdgecolor','g');
plot(7,6.6, 'ks', 'linewidth',1,'markersize' ,15,'MarkerFacecolor','r','MarkerEdgecolor','r');
xlim([1, 10])
ylim([1,10])

text(7.3,8.7,{'wave' , 'death'} , 'FontSize',15, 'interpreter','latex');
text(7.3,7.5,{'radial', 'prop.'} , 'FontSize',15, 'interpreter','latex');
text(7.3,6.6,{'spiral'} , 'FontSize',15, 'interpreter','latex');

% leg1={ [' on' char(10) '($$\tau_{on}$$)'], [' off' char(10) '($$\tau_{off}=5$$)']};
% leg2={ ['      ' char(10) '      '], ['      ' char(10) '      ']};

%legend(gca,{'firing cell ($$\tau_{on}$$)','non-firing cell ($$\tau_{off}=5$$)'}, 'FontSize',20, 'interpreter','latex');
%legend(gca,{' on($$\tau_{on}$$)', ' off($$\tau_{off}=5$$)'}, 'FontSize',15, 'interpreter','latex');

%legend(gca,leg2, 'FontSize',15, 'interpreter','latex');
%legend boxoff 
print(gcf,['../PapFig/Cartoon_LegOnOff_3.png'],'-dpng','-r500');

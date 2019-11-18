
a=[ 55   52    63     73    82    88   79     70   54
;37    28   22     12    15     25  31     23   36
];

figure;
plot(a(1,:),a(2,:), '-O', 'color', 'c', 'LineWidth', 4, 'MarkerSize', 7,'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
xlim([1,100])
ylim([1,100])
set(gca,'YDir','reverse')

fprintf('%d, ', a(1,:)-1);
fprintf('\n');
fprintf('%d, ', a(2,:)-1);
%print(gcf,['../PapFig/Spi7_R',num2str(repeat), 'Rho0t240p420c10_SPIHead.png'],'-dpng','-r500');

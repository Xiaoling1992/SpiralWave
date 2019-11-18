clear all;

height=100;
width=100;
ColMap=[1,1,1; 0 0 0]
OCC=load('/home/xiaoling/SpiralWave/data/Spi11Rho95Phi500t9c10_OCCin.dat');
Y=floor(OCC/width);
X=rem(OCC, width);

Y=Y+1;
X=X+1;

SpiHead=[55   52    63     73    82    88   79     70   54
;37    28   22     12    15     25  31     23   36];
     
figure;
biofilm=zeros(height, width);

for i=1:length(X)
    biofilm(Y(i), X(i))=1;
end

imagesc(biofilm);
hold on;
colormap(ColMap);

plot(SpiHead(1,:), SpiHead(2, :), '-O', 'color', [0,1,0], 'LineWidth', 4, 'MarkerFaceColor', [0,1,0], 'MarkerEdgeColor', [0,1,0]);
title('$$\phi=0.50$$', 'fontsize',25,'interpreter','latex','FontWeight', 'bold')

 set(gca,'xtick',[])
 set(gca, 'ytick',[])
   %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   print(gcf,['../PapFig/Spi11Rho95Phi500t9c10_OCCinWF.png'],'-dpng','-r500');
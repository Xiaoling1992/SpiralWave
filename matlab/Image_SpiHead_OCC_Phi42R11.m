
height=100;
width=100;
repeat=11;
ColMap=[1,1,1; 0 0 0]
OCC=load(['/home/xiaoling/SpiralWave/data/Spi7_Rho0t240p420c10R',num2str(repeat),'_OCCin.dat']);
Y=floor(OCC/width);  %OCC= i* width+ j in C
X=rem(OCC, width);

Y=Y+1;  %convert to index in MatLab
X=X+1;
%Where is my explanation? [X, Y] or [j, i] in the 2D matrix
SpiHead=[62    66    70    69    75    69    60    56    55    59    63;
    49    51    59    66    77    75    77    71    62    54    48];
     
figure;
biofilm=zeros(height, width);

for i=1:length(X)
    biofilm(Y(i), X(i))=1;
end

figure
imagesc(biofilm);
hold on;
colormap(ColMap);

plot(SpiHead(1,:), SpiHead(2, :), '-O', 'color', 'c', 'LineWidth', 4, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c');
%title('$$\phi=0.42$$', 'fontsize',25,'interpreter','latex','FontWeight', 'bold')

 set(gca,'xtick',[])
 set(gca, 'ytick',[])
   %axis equal
%    set(gca,'xtick',[1:length(phis)],'xticklabel',phis([1:length(phis)]),'ytick',1:length(taus),'yticklabel',taus);
   %print(gcf,'-depsc',['../PapFig/Spi423_bRho75u2200s57c10_T',num2str( time(i) ),'.eps']);
   %print(gcf,['../PapFig/Spi7_R',num2str(repeat),'_Rho0t240p420c10_OCC_4.png'],'-dpng','-r500');
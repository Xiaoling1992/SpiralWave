clear all;

ColMap=[
       0 0 1
       0 0.5 0.5
       0.5,0.5,1       
       1 0 0
       ];


taus=[10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data= load(['/home/xiaoling/SpiralWave/data/Spi7_b_Phase_Rho75.dat']);
phase_Pre=zeros(length(taus)*length(phis), 4); % Death; Radial; Break up spiral; Secondary;
phase_Pre2=zeros(length(taus),length(phis),4);

for i=1:length(taus)*length(phis)
    phase_Pre(i,1)= sum( data(i,:)==0 );
    phase_Pre(i,2)= sum( data(i,:)==3 )+sum( data(i,:)==5 )+sum( data(i,:)==9 );
    phase_Pre(i,3)= sum( data(i,:)==2 );
    phase_Pre(i,4)= sum( data(i,:)==8  );
end
titles={'Death', 'Radial', 'BS', 'SS'}
names= ['De';'Ra';'BS';'SS'];
phase_Pre2= reshape( phase_Pre, [ length(phis), length(taus), 4] );
for i = 1:4
figure
imagesc(phase_Pre2( :, : , i)' );
hold on
plot(find(phis==0.46), find(taus==530), 'kX', 'LineWidth', 3, 'markersize' ,15,'MarkerEdgeColor','k' )
 plot(find(phis==0.49), find(taus==19), 'kS', 'LineWidth',3,'markersize' ,15,'MarkerEdgeColor','k' )
xlabel('$$\phi$$','fontsize',20,'FontWeight', 'bold','interpreter','latex');
ylabel('$$\tau_{on}$$','fontsize',20,'FontWeight', 'bold','interpreter','latex');
set(gca,'xtick',[1:3: length(phis)],'xticklabel',phis([1:3: length(phis)]),'ytick',1:3:length(taus),'yticklabel',taus(1:3:length(taus) ));  
set(gca,'ydir','normal')
zlim([0, 20])
colorbar( )

title(titles(i) )

print(gcf,['../PapFig/Image_Fre', names(i, :) ,'_TauPhi_Rho75_BS_Review.png'],'-dpng','-r500');
end

%colormap(ColMap);
%cb =colorbar('Ticks',[1, 2 ,3, 4],'TickLabels',{'D','R' ,'BS','SS'} );
%set(cb,'position',[.9 0.5 0.05 0.4])
%c.Label.String = 'phase';  
%   caxis([0,9])


 
 
%c =colorbar('Ticks',[1, 2 ,3, 4],'TickLabels',{'D','R' ,'BS','SS'} );
%c.Label.String = 'phase';  
%   caxis([0,9])


 
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);








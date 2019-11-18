%Spi10_Phase* Review_4.c : FirSamMean_th=0.5  
%Spi10_Phase* Review_5.c : FirSamMean_th=phi

clear all;

ColMap=[
       0 0 1
       0 0.5 0.5
       0.5,0.5,1       
       1 0 0
       ];
   
myCol= zeros( [21, 3] );
for i = 1:21
    
    myCol(i, 2)= i/21;
end
    


taus=[5, 6, 7, 8, 9, 10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data= load(['/home/xiaoling/SpiralWave/data/Spi10_Phase_TauPhi_Rho0_5.dat']);
phase_Pre=zeros(length(taus)*length(phis), 4); % Death; Radial; Break up spiral; Secondary spiral;
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

philabels=[0.2, 0.4, 0.6, 0.8, 1.0];
[P, S]= polyfit( phis, 1: 1: length(phis), 1 );

taulabels= [10, 100, 1000];
[P1, S1]= polyfit( log10(taus), 1: 1: length(taus), 1 );

figure( )
set(gcf, 'PaperUnits', 'inches')
set(gcf, 'PaperSize', [5, 5])
for i = 1:4
ax= subplot(2, 2, i);
imagesc(phase_Pre2( :, : , i)' );
%colormap( parula)
%zlim manual

max_fre= max( phase_Pre(:, i) );
temCol= myCol(1: 1: max_fre+1, :);

colormap(ax, temCol )     %If I don't add ax as a input parameter of colormap, I will apply this commond to every figure
%zlim(ax, [0, max_fre] ) 
%colorbar
%colormap( myCol )
%zlim(ax, [0, 20])

%plot(find(phis==0.46), find(taus==530), 'kX', 'LineWidth', 3, 'markersize' ,15,'MarkerEdgeColor','k' )
%plot(find(phis==0.49), find(taus==19), 'kS', 'LineWidth',3,'markersize' ,15,'MarkerEdgeColor','k' )

set(gca,'xtick', P(1)* philabels + P(2),'xticklabel', philabels,'ytick', P1(1)* log10(taulabels) + P1(2),'yticklabel', taulabels , 'fontsize',8);  
title(titles(i), 'fontsize',12,'FontWeight', 'normal' )
xlabel('$$\phi$$','fontsize',15,'interpreter','latex');
ylabel('$$\tau_{on}$$','fontsize',15,'interpreter','latex');
set(gca,'ydir','normal')





%colorbar( )




end

print(gcf,['../PapFig/Image_Fre_TauPhi_Rho0_Spi10_5.png'],'-dpng','-r500');









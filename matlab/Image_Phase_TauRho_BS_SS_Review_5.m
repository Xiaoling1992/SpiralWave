clear all;

ColMap=[
       0 0 1
       0 0.5 0.5
       255/255, 192/255, 203/255  %pink 
       1 0 0
       ];


taus=[5,6,7,8,9,10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data= load(['/home/xiaoling/SpiralWave/data/Spi10_Phase_TauPhi_Rho50_5.dat']);
phase_Pre=zeros(length(taus)*length(phis), 4); % Death; Radial; Break up spiral; Secondary;
phase_Pre2=zeros(length(taus)*length(phis),1);

for i=1:length(taus)*length(phis)
    phase_Pre(i,1)= sum( data(i,:)==0 );
    phase_Pre(i,2)= sum( data(i,:)==3 )+sum( data(i,:)==5 )+sum( data(i,:)==9 );
    phase_Pre(i,3)= sum( data(i,:)==2 );
    phase_Pre(i,4)= sum( data(i,:)==8  );
    
    if (phase_Pre(i,3)>0 || phase_Pre(i,4) >0)
         [ tem,phase_Pre2(i)]=max( [phase_Pre(i, 3), phase_Pre(i, 4)] );
         phase_Pre2(i)=phase_Pre2(i) +2;
    else        
        [ tem,phase_Pre2(i)]=max( phase_Pre(i, :) );
    end
end

phase=reshape(phase_Pre2, [length(phis), length(taus)]);
phase=phase';

figure
imagesc(phase);
hold on

colormap(ColMap);
%cb =colorbar('Ticks',[1, 2 ,3, 4],'TickLabels',{'D','R' ,'BS','SS'} );
%set(cb,'position',[.9 0.5 0.05 0.4])
%c.Label.String = 'phase';  
caxis([0.5,4.5])





philabels=[0.2, 0.4, 0.6, 0.8, 1.0];
[P, S]= polyfit( phis, 1: 1: length(phis), 1 );

taulabels= [10, 100, 1000];
[P1, S1]= polyfit( log10(taus), 1: 1: length(taus), 1 );


 
 %plot(find(phis==0.46), find(taus==530), 'kX', 'LineWidth', 3, 'markersize' ,15,'MarkerEdgeColor','k' )
 % plot(find(phis==0.49), find(taus==19), 'kS', 'LineWidth',3,'markersize' ,15,'MarkerEdgeColor','k' )
%c =colorbar('Ticks',[1, 2 ,3, 4],'TickLabels',{'D','R' ,'BS','SS'} );
%c.Label.String = 'phase';  
%   caxis([0,9])

%text( find(phis== 0.37)-.5, find(taus==452), 'E','FontSize',20, 'FontWeight','bold','color','k' )
%text( find(phis== 0.64)-1.5, find(taus==14), ' F', 'FontSize',20,'FontWeight','bold','color','k' )

set(gca,'xtick', P(1)* philabels + P(2),'xticklabel', philabels,'ytick', P1(1)* log10(taulabels) + P1(2),'yticklabel', taulabels , 'fontsize',14);  
set(gca,'ydir','normal')
xlabel('$$\phi$$','fontsize',22,'FontWeight', 'bold','interpreter','latex');
ylabel('$$\tau_{on}$$','fontsize',22,'FontWeight', 'bold','interpreter','latex');
%text( 0.64-0.02, find(taus==14), 'D', 'FontSize',20,'FontWeight','bold','color','k' )

 
    
    
    % print('-clipboard','-dbitmap');
 %print(gcf,'-depsc',['../figure/bioPhi55Rho22a0aL0.1b0.004.eps']);
print(gcf,['../PapFig/Image_Spi10_Phase_TauPhi_Rho50_5_1.png'],'-dpng','-r700');







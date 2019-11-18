clear all;
pi=3.1415926;
epsilon=10;
phis=[0.26, 0.28,0.30,0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.50,0.52,0.54,0.56,0.58];
periodPre=load('/home/xiaoling/SpiralWave/data/Spi7_n5_Period_t240.dat');
anchorPre=load('/home/xiaoling/SpiralWave/data/Spi7_n5_Anchor_t240.dat');


for i=1: length(phis)
    period_i=periodPre(i,:);
    periodEff=period_i(period_i >0 )
    
    anchor_i=anchorPre(i,:)
    anchorEff=anchor_i( period_i >0)
    anchorEff=sqrt(anchorEff); %get the radius
     
    if length(periodEff)>0
            period(1,i)=mean(periodEff);
            period(2,i)=std(periodEff)/sqrt(length(periodEff)) ;
            RatioEff(1, i)=length(periodEff) / length(period_i);

            anchor(1,i)=mean(anchorEff);
            anchor(2,i)=std(anchorEff)/sqrt(length(anchorEff));
            RatioEff(2, i)=length(anchorEff) / length(anchor_i) ; 
    end
    

end

figure;
plot(phis( phis>0.3 ), period(1,phis>0.3 ), '-O', 'linewidth',3,'markersize' ,7,'MarkerFacecolor','c', 'Color','c');
hold on
plot(phis(phis>0.3 ), anchor(1,phis>0.3 )*pi/sqrt(epsilon/2), '-O', 'linewidth',3,'markersize' ,7, 'MarkerFacecolor','m','Color','m');
legend(gca, {'\ dynamical',' $$\ \pi \sqrt{n_{off}}/ \sqrt{\frac{\epsilon}{2} }$$'}, 'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex');
legend boxoff
plot([0.41, 0.41], [0, 120], '--', 'LineWidth', 2, 'Color', 'k');

%legend(gca, {'N=100','N=300','N=1000'},  'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex' )

%xlim([2*10^(-2), 2.2*10^(-1)]);
%ylim([2, 50]);
set(gca,'xtick',[0.3, 0.35, 0.41, 0.45, 0.5, 0.55, 0.6],'xticklabel',{0.3, 0.35, 'p_c', 0.45, 0.5, 0.55, 0.6}, 'fontsize', 15);
xlabel('$$\phi$$', 'fontsize',25,'interpreter','latex','FontWeight', 'bold')
ylabel('period, $$T$$', 'fontsize', 25, 'interpreter','latex', 'FontWeight', 'bold')
%set(gca, 'xticks', 0.41,'xticklabel','$$p_c$$', 'fontsize', 15, 'FontWeight', 'bold','interpreter','latex' )


%text(0.48, 42, 'dynamical (data)', 'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.40, 80, '$$\pi d_{struct}/ \sqrt{\frac{\epsilon}{2} }$$',  'FontSize',15, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.35, 60, '$$\downarrow$$',  'Color','r','FontSize',70, 'FontWeight', 'bold','interpreter', 'latex')
%text(0.5, 25, '$$\uparrow$$', 'Color','r', 'FontSize',40, 'FontWeight', 'bold','interpreter', 'latex')

print(gcf,['../PapFig/Lin_PerStruAnc_Phi_1.png'],'-dpng','-r400');



clear all
time=0:0.5:200;
height=100;
width=200;
WinHei=35;

as=[0.02, 0.05, 0.02, 0.02];  %WT, trkA, sinR, ktrA
taus=[300, 550, 300, 650];  
cs=[10, 10, 10, 10];
phis=[0.43, 0.43, 0.74, 0.51];


name={'WT';'Delta trkA';'Delta sinR';'Delta ktrA'};

row=13;
v_th=0.6;
Tab=load('/home/xiaoling/SpiralWave/data/F17_4Strains_ACell_T.dat');

for i=1:length(taus)
    
    
    phi=phis(i);
    a=as(i);
    tau=taus(i);
    c=cs(i);
    
    
    TBio=Tab([1:WinHei]+(i-1)*WinHei,:);
    wave=load(['/home/xiaoling/SpiralWave/data/F17ACPhi', num2str(phi*100),'Rho0a',...
        num2str(a*100),'t', num2str(tau),'tL5c10_wave.dat']);
    
    %select 100 samples;    
    state=zeros(1, width);   
    for j=1:1:100
        sample(j)=randi([1 width],1,1);
        while state( sample(j) ) ~=0
            sample(j)=randi([1 width],1,1);
        end
        state(sample(j))=1;
    end
    
    TSam=TBio(row-1, sample);  %the periods of the picked up cells. Tab started from the secod row. So the row-1 in Tab is the row in wave
    [B, I] = sort(TSam,'descend'); %Sort the TSam 
    
    SamSor=sample(I);
    trace=wave([0: (length(time)-1)]*height+row,SamSor);
    
    dlmwrite(['../PapFig/F17ACPhi', num2str(phi),'Rho0a',...
        num2str(a),'t', num2str(tau),'tL5c10.dat'], name{i}, 'delimiter', '');
     dlmwrite(['../PapFig/F17ACPhi', num2str(phi),'Rho0a',...
        num2str(a),'t', num2str(tau),'tL5c10.dat'], 'Each column is a trace of one cell. time=0:0.5:200','-append', 'delimiter', '');
    dlmwrite(['../PapFig/F17ACPhi', num2str(phi),'Rho0a',...
        num2str(a),'t', num2str(tau),'tL5c10.dat'],trace, '-append', 'delimiter', '\t', 'precision', 4);
    
    
%     figure()
%     imagesc(1:100, time(time<100), trace(time<100,:));
%     colormap(jet);
%     
%     xlabel('single cell ThT traces (trace ID)', 'fontsize',20)
%     ylabel('time', 'fontsize', 20)
%     
%      %title(name(i),'position',[1 1]);
%      
%     print(gcf,['../PapFig/F17ACPhi', num2str(phi*100),'Rho0a',...
%         num2str(a*100),'t', num2str(tau),'tL5c10.png'],'-dpng','-r400');
    
    
    
end



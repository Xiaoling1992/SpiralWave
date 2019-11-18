tau=22;
phi=0.58;

taus=[10, 12, 14, 16, 19, 22, 26, 30, 36, 42, 49, 57, 67, 79, 92, 108, 127, 149, 174, 204, 240, 281, 329, 386, 452, 530, 621, 728, 853, 1000];
phis=[0.13, 0.16, 0.19, 0.22, 0.25, 0.28, 0.31, 0.34,0.37, 0.40, 0.43, 0.46, 0.49, 0.52, 0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.88, 0.91,0.94, 0.97, 1.00];

data= load(['/home/xiaoling/SpiralWave/data/Spi7_b_Phase_Rho75.dat']);
phase_Pre=zeros(length(taus)*length(phis), 4); % Death; Radial; Break up spiral; Secondary;
phase_Pre2=zeros(length(taus)*length(phis),1);

for i=1:length(taus)*length(phis)
    
    itau=floor( (i-1)/length(taus) )+1;
    iphi=rem(i-1, length(taus)) +1;
    
    if( taus(itau) >=0 )
        phase_Pre(i,1)= sum( data(i,:)==0 );
        phase_Pre(i,2)= sum( data(i,:)==3 )+sum( data(i,:)==5 );
        phase_Pre(i,3)= sum( data(i,:)==2 );
        phase_Pre(i,4)= sum( data(i,:)==8 )+sum( data(i,:)==9 );
    else
        phase_Pre(i,1)= sum( data(i,:)==0 );
        phase_Pre(i,2)= sum( data(i,:)==3 )+sum( data(i,:)==5 )+sum( data(i,:)==8 );
        phase_Pre(i,3)= sum( data(i,:)==2 );
        phase_Pre(i,4)= sum( data(i,:)==9 );
    end
    
    if (phase_Pre(i,3)>0 || phase_Pre(i,4) >0)
         [ tem,phase_Pre2(i)]=max( [phase_Pre(i, 3), phase_Pre(i, 4)] );
         phase_Pre2(i)=phase_Pre2(i) +2;
    else        
        [ tem,phase_Pre2(i)]=max( phase_Pre(i, :) );
    end
end

phase=reshape(phase_Pre2, [length(phis), length(taus)]);
phase=phase';


itau= find(taus== tau);
iphi= find(phis==phi);
order= (itau-1)*length(phis) +iphi;
data(order, :)
phase(itau, iphi)
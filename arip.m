function [ARI, mV] = arip(cA,t,fs);
%P can also be the CBFV-ABP dynamic relationship change in ABP during
%spontaneous fluctuations in ABP
%% Tieck' Model Parameters
dados = 'tiecksparameters.xlsx';
local = '';
para = importdata([local,dados]);
para = para.data;

T = para(:,1); %time in seconds
D = para(:,2); %dumping factor
K = para(:,3); %gain
ARI = para(:,4); %ari index
%vt = t/fs; %vetor velocidade
%% Functions
CCP = 12; %critical closure pressure [mm Hg]
ABP = mean(cA);
% cVmean = mean(CBFv_1); %media velocidade do fluxo cerebral

    for n = 1:length(cA)
    
    dP(n) = (cA(n)-ABP / ABP - CCP); %Normalizacao da mean pressure
    x1(1,1:10) = 0; %condicoes iniciais
    x2(1,1:10) = 0; %condicoes iniciais
    
        for i = (1:10)

            x1(n+1,i) = x1(n,i) + (dP(n)-x2(n,i))/fs*T(i); %ED 2 ordem
            x2(n+1,i) = x2(n,i) + (x1(n,i)-2*D(i)*x2(n,i))/fs*T(i); %ED 2 ordem
            mV(n,i) = 1+dP(n)'-K(i)*x2(n,i); %velocidade media
    
        end
    end
   
%     figure(1)
%     plot(mV)
%     xlabel('time(s)')
%     ylabel('flow velocity')
%     xlim([0 15])
%     legend('ARI 0', 'ARI 1','ARI 2','ARI 3','ARI 4','ARI 5','ARI 6','ARI 7','ARI 8','ARI 9')

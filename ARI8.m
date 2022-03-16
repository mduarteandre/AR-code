function [ARI8] = ARI8(ABP,CBFV,fs);
    mABP = mean(ABP);
    Vm = mean(CBFV);
    vt = [0:length(ABP)-1]/fs;
    CCP = 12;

       for n = 1:length(vt)
        
%   dP(n) = (ABP(n) - mABP  / mABP - CCP); %Normalizacao da mean pressure
     dP(n) = (ABP(n)  / 1 - CCP);
 
    
    x1(1,1:10) = 0; %condicoes iniciais
    x2(1,1:10) = 0; %condicoes iniciais
   
   T = 0.87;
   D = 0.52;
   K = 0.97;
   
   x1(n+1) = x1(n) + (dP(n)-x2(n)/fs*T); %ED 2 ordem
   x2(n+1) = x2(n) + (x1(n)-2*D*x2(n)/fs*T); %ED 2 ordem
%    ARI8(n) = Vm*(1+dP(n)-K*x2(n)); %velocidade medis
   ARI8(n) = 1+dP(n)-K*x2(n); 
       end

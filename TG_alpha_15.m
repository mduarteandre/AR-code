clear all
clc

%% Selecao de dados
arquivo = 'P04BA1A1.DAT';
diretorio = '';

dados = importdata([diretorio,arquivo]);
dados = dados.data; %extrai apenas numeros

%% Extracao de dados
fs = 100; %Hz
t = dados(:,1)/fs; %tempo
CBFv_1_left = dados(:,2); %velocidade do fluxo cerebral esq
CBFv_1_right = dados(:,3);
ABP_raw = dados(:,8); % pressao arterial

fs_calculated=1/mean(diff(t)); % find the sampling frequency from the time-base


%% Renomear

CBFv_l = CBFv_1_left;
CBFv_r = CBFv_1_right;
ABP_v = ABP_raw;

%% Achando os pontos Beat-to-beat 

sinal = CBFv_l;
distMin = 50; % dist min entre picos
tipo = 'negativo'; % positivo ou negativo
[indicePicos_l,alturaPicos_l] = achaPicos(sinal, distMin, tipo);

sinal = CBFv_r;
% distMin = 70; % dist min entre picos
tipo = 'negativo'; % positivo ou negativo
[indicePicos_r,alturaPicos_r] = achaPicos(sinal, distMin, tipo);

sinal = ABP_v;
distMin2= 40; % distMin = 50; % dist min entre picos
tipo = 'negativo'; % positivo ou negativo
[indicePicos_ABP,alturaPicos_ABP] = achaPicos(sinal, distMin2, tipo);

%% Media do Sinal CBFv Left
Pos = 0;
beat = 0;
for ii=1:length(indicePicos_l)-1;   
    beat = beat+1;      
        if beat ==1
        Pos (beat)= 0 + (indicePicos_l(ii+1)-indicePicos_l(ii));
        Mean_beat = mean(CBFv_l(indicePicos_l(ii):indicePicos_l(ii+1))); 
        Matriz_final (1, 1:2) = [Pos(1), Mean_beat];
        else
        Pos (beat)= Pos(beat-1) + (indicePicos_l(ii+1)-indicePicos_l(ii));
        Mean_beat = mean(CBFv_l(indicePicos_l(ii):indicePicos_l(ii+1))); 
        Matriz_final (beat, 1) = Pos (beat); 
        Matriz_final (beat, 2) = Mean_beat; 
        end
end

%keyboard
nova_freq = length(Pos)/Pos(length(Pos)); %HZ

%% Interpolation Spline CBFv Left

x_l = Matriz_final (:,1);%intervalo de valore no eixo X do sinal após a media beat to beat
y_l = Matriz_final (:,2); % valores do sinal em y a serem interpolados 
xx_l = Matriz_final (1,1)/100:0.2:floor(x_l(end)/fs); %novo intervalo de valores no eixo x do sinal em uma taxa de geraçao de novos pontos de 5 Hz
yy_l = spline(x_l/100,y_l,xx_l);%plot em segundos

%% Media do Sinal CBFv Right
Pos_r = 0;
beat_r = 0;
for ii=1:length(indicePicos_r)-1;   
    beat_r = beat_r+1;      
        if beat_r ==1
        Pos_r (beat_r)= 0 + (indicePicos_r(ii+1)-indicePicos_r(ii));
        Mean_beat_r = mean(CBFv_r(indicePicos_r(ii):indicePicos_r(ii+1))); 
        Matriz_final_r (1, 1:2) = [Pos_r(1), Mean_beat_r];
        else
        Pos_r (beat_r)= Pos_r(beat_r-1) + (indicePicos_r(ii+1)-indicePicos_r(ii));
        Mean_beat_r = mean(CBFv_r(indicePicos_r(ii):indicePicos_r(ii+1))); 
        Matriz_final_r (beat_r, 1) = Pos_r (beat_r); 
        Matriz_final_r(beat_r, 2) = Mean_beat_r; 
        end
end
 
 nova_freq_r = length(Pos_r)/Pos_r(length(Pos_r)); %HZ
 
 %% Interpolation Spline CBFv Right

x_r = Matriz_final_r (:,1);%intervalo de valore no eixo X do sinal após a media beat to beat
y_r = Matriz_final_r (:,2); % valores do sinal em y a serem interpolados 
xx_r = Matriz_final_r (1,1)/100:0.2:floor(x_r(end)/fs); %novo intervalo de valores no eixo x do sinal em uma taxa de geraçao de novos pontos de 5 Hz
yy_r = spline(x_r/100,y_r,xx_r);

%% Media do Sinal ABP
Pos_abp = 0;
beat_abp = 0;
for ii=1:length(indicePicos_ABP)-1;   
    beat_abp = beat_abp+1;      
        if beat_abp ==1
        Pos_abp (beat_abp)= 0 + (indicePicos_ABP(ii+1)-indicePicos_ABP(ii));
        Mean_beat_abp = mean(ABP_v(indicePicos_ABP(ii):indicePicos_ABP(ii+1))); 
        Matriz_final_abp (1, 1:2) = [Pos_abp(1), Mean_beat_abp];
        else
        Pos_abp (beat_abp)= Pos_abp(beat_abp-1) + (indicePicos_ABP(ii+1)-indicePicos_ABP(ii));
        Mean_beat_abp = mean(ABP_v(indicePicos_ABP(ii):indicePicos_ABP(ii+1))); 
        Matriz_final_abp (beat_abp, 1) = Pos_abp (beat_abp); 
        Matriz_final_abp(beat_abp, 2) = Mean_beat_abp;
        %Mean_ARI_ABP = Mean_beat_abp;
        end
end

for n= 1:length ( Matriz_final_abp)-1;
    mABP(n) =   Matriz_final_abp(n+1,1) - Matriz_final_abp(n,1);
    
end

nova_freq_abp = length(Pos_abp)/Pos_abp(length(Pos_abp)); %HZ
 
 %% Interpolation Spline ABP

x_abp = Matriz_final_abp (:,1);%intervalo de valore no eixo X do sinal após a media beat to beat
y_abp = Matriz_final_abp (:,2); % valores do sinal em y a serem interpolados 
xx_abp = Matriz_final_abp (1,1)/100:0.2:floor(x_abp(end)/fs); %novo intervalo de valores no eixo x do sinal em uma taxa de geraçao de novos pontos de 5 Hz
yy_abp = spline(x_abp/100,y_abp,xx_abp);


%% Matrizão + Deixando os vetores do mesmo tamanho

aaa = xx_abp (:); %Para fazer o teste no CARNET como eixo de tempo nos dados interpolados
aaa1 = yy_r (:);
aaa2 = yy_l (:);
aaa3 = yy_abp (:);



aa =length (aaa); bb = length (aaa1); cc = length (aaa2); dd = length (aaa3);
Size_test = [aa, bb, cc, dd]; adjust = min (Size_test); 

x_axis = aaa(1:adjust); 
Cbfv_r_Carnet = aaa1(1:adjust); 
Cbfv_l_Carnet = aaa2(1:adjust); 
abp_carnet = aaa3(1:adjust); 

yy_r_2 = Cbfv_r_Carnet; 
yy_l_2 = Cbfv_l_Carnet; 
yy_abp_2 = abp_carnet;
xx_abp_2 = x_axis;



matrizao_P04_AVC = [x_axis abp_carnet Cbfv_l_Carnet Cbfv_r_Carnet];

keyboard
%% Média

ABP=yy_abp-mean(yy_abp); %LINHA 118 DA TFA_CAR, RETIRA A MEDIA E USA ESSE VALOR PARO O WELCH
CBFV_l=yy_l-mean(yy_l);
CBFV_r=yy_r-mean(yy_r);

%keyboard

%% Welch

window = hann(512); %white paaper diz hanning acima de 100s
window_length = length(window);
M = window_length;

fs_2 = 5; % HZ nova frequencia do sinal = amostragem

%Nosso código: 
overlap = 0.5078; %50/100; %porcentagem
%Carnet:
%overlap=(window_length-shift)/window_length*100;

Nfft = M;
shift = round((1-overlap)*M); %passo dos segmentos overlapado
M = M(:);%transforma em vetor coluna
N = length(abp_carnet); %comprimento do sinal


x = ABP(:); %transforma em vetor coluna
y = CBFV_l(:); %transforma em vetor coluna
keyboard
%b = CBFV_l(:); %transforma em vetor coluna

X = fft(x(1:M).*window); %fft do primeiro passo
Y = fft(y(1:M).*window);%fft do primeiro passo
%B = fft(b(1:M).*window); %fft do primeiro passo

L = 1; % contador de numero de passos/ no_windows

if shift>0
    i_start=1+shift;
    while i_start+M-1 <= N
      X=[X,fft(x(i_start:i_start+M-1).*window,Nfft)];
      Y=[Y,fft(y(i_start:i_start+M-1).*window,Nfft)];
      i_start=i_start+shift;
      L = L+1; 
      
    end  
   
end
f=[0:Nfft-1]'/Nfft*fs_2; %ajuste de resolucao para ficar em FREQUENCIA
figLenghth = Nfft/2+1;

no_windows = L;

%% Auto-Spectro and Cross-Spectro 

%NAO CONSIDERA PARA L == 1 COMO NO CARNET

D.Pxx=sum(X.*conj(X),2)/L/sum(window.^2)/fs_2; %D serve para armazenarmos os valores de auto/cross e correlaçao 
D.Pyy=sum(Y.*conj(Y),2)/L/sum(window.^2)/fs_2; 
D.Pxy=sum(conj(X).*Y,2)/L/sum(window.^2)/fs_2;
D.coh=D.Pxy./((abs((D.Pxx.*D.Pyy))).^0.5);

Pxx=D.Pxx;
Pyy=D.Pyy;
Pxy=D.Pxy;

%keyboard

%% Triangular Moving Average

e = [1/4 1/2 1/4];

%keyboard

Pxx1=Pxx;
Pxx1(1)=Pxx(2);

Pyy1=Pyy;
Pyy1(1)=Pyy(2);

Pxy1=Pxy;
Pxy1(1)=Pxy(2);

Pxx1raw = Pxx1;
Pyy1raw = Pyy1;
Pxy1raw = Pxy1;


Pxx1=filtfilt(e,1,Pxx1raw);
Pyy1=filtfilt(e,1,Pyy1raw);
Pxy1=filtfilt(e,1,Pxy1raw);

%keyboard

%% Transfer Function 
FT1 = Pxy1./Pxx1; %FT
FT = abs(FT1);
H = FT1;

t_h = [0:length(FT)-1]/fs_2;


keyboard

%% Coherence 
C=Pxy1./(abs(Pxx1.*Pyy1).^0.5); %C É imaginário

%keyboard

%% Thresholds

params.coherence2_thresholds=[3:15;0.51,0.40,0.34,0.29,0.25,0.22,0.20,0.18,0.17,0.15,0.14,0.13,0.12]';
params.remove_negative_phase=1;
params.remove_negative_phase_f_cutoff=0.1;
params.apply_coherence2_threshold=1;
params.vlf=[0.02,0.07];
params.lf=[0.07,0.2];
params.hf=[0.2,0.5];
params.normalize_ABP=0;
params.normalize_CBFV=0;


%% Parametros

i=find(params.coherence2_thresholds(:,1)==no_windows);
if isempty(i)
disp('Warning:no coherence threshold defined for the number of windows obtained - all frequencies will be included');
coherence2_threshold=0;
else
    coherence2_threshold=params.coherence2_thresholds(i,2);
end

G=H; % save for plotting below
%keyboard
if params.apply_coherence2_threshold
i=find(abs(C).^2 < coherence2_threshold); % exclude low coherence
cA = abs(C).^2;
H(i)=nan;

%keyboard

end

P=angle(H);
% keyboard
if params.remove_negative_phase; % exclude negative phase below cut-off frequency
    n=find(f<params.remove_negative_phase_f_cutoff);
    k=find(P(n)<0);
%    keyboard
    if ~isempty(k);
        P(n(k))=nan;
    end
end;


%keyboard

i=find(f>=params.vlf(1) & f<params.vlf(2));
tfa_out.Gain_vlf=nanmean(abs(H(i)));
tfa_out.Phase_vlf=nanmean(P(i))/(2*pi)*360;
tfa_out.Coh2_vlf=nanmean(abs(C(i)).^2);
tfa_out.P_abp_vlf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_vlf=2*sum(Pyy(i))*f(2);

%keyboard 

i=find(f>=params.lf(1) & f<params.lf(2));
tfa_out.Gain_lf=nanmean(abs(H(i)));
tfa_out.Phase_lf=nanmean(P(i))/(2*pi)*360;
tfa_out.Coh2_lf=nanmean(abs(C(i)).^2);
tfa_out.P_abp_lf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_lf=2*sum(Pyy(i))*f(2);
tfa_out.Mean_cbfv=mean(yy_r);

i=find(f>=params.hf(1) & f<params.hf(2));
tfa_out.Gain_hf=nanmean(abs(H(i)));
dummy=angle(H(i));
tfa_out.Phase_hf=nanmean(P(i))/(2*pi)*360;
tfa_out.Coh2_hf=nanmean(abs(C(i)).^2);
tfa_out.P_abp_hf=2*sum(Pxx(i))*f(2);
tfa_out.P_cbfv_hf=2*sum(Pyy(i))*f(2);

if params.normalize_CBFV
    tfa_out.Gain_vlf_norm=tfa_out.Gain_vlf;
    tfa_out.Gain_lf_norm=tfa_out.Gain_lf;
    tfa_out.Gain_hf_norm=tfa_out.Gain_hf;
    tfa_out.Gain_vlf_not_norm=tfa_out.Gain_vlf*tfa_out.Mean_cbfv/100;
    tfa_out.Gain_lf_not_norm=tfa_out.Gain_lf*tfa_out.Mean_cbfv/100;
    tfa_out.Gain_hf_not_norm=tfa_out.Gain_hf*tfa_out.Mean_cbfv/100;
else
    tfa_out.Gain_vlf_not_norm=tfa_out.Gain_vlf;
    tfa_out.Gain_lf_not_norm=tfa_out.Gain_lf;
    tfa_out.Gain_hf_not_norm=tfa_out.Gain_hf;
    tfa_out.Gain_vlf_norm=tfa_out.Gain_vlf/tfa_out.Mean_cbfv*100;
    tfa_out.Gain_lf_norm=tfa_out.Gain_lf/tfa_out.Mean_cbfv*100;
    tfa_out.Gain_hf_norm=tfa_out.Gain_hf/tfa_out.Mean_cbfv*100;
end
    
%keyboard

%% Plots 

% Sinal sem tratamento

figure (1)
subplot(3,1,1)
plot(t,CBFv_1_left,'r')
xlabel('time (s)');
ylabel('cm/s')
title('CBFv_l - Sinal Não Processado')
axis tight

subplot(3,1,2)
plot(t,CBFv_1_right,'b')
xlabel('time (s)');
ylabel('cm/s')
title('CBFv_r - Sinal Não Processado')
axis tight

subplot(3,1,3)
plot(t,ABP_raw,'k')
xlabel('time (s)');
ylabel('mmHg')
title('ABP - Sinal Não Processado')
axis tight

%keyboard

% Moving-Average Filter (encontra a média de séries de diferentes tamanhos)
 
% figure (2)
% subplot(3,1,1)
% plot(t,CBFv_1_mvl,'r')
% xlabel('time (s)');
% ylabel('cm/s')
% title('CBFv_l - Moving Averege Filter')
% axis tight
% 
% subplot(3,1,2)
% plot(t,CBFv_1_mvr,'b')
% xlabel('time (s)');
% ylabel('cm/s')
% title('CBFv_r - Moving Averege Filter')
% axis tight
% 
% subplot(3,1,3)
% plot(t,ABP_mv,'g')
% xlabel('time (s)');
% ylabel('mmHg')
% title('ABP - Moving Averege Filter')
% axis tight

%keyboard

% Achando os pontos Beat-to-beat 

figure(3)
subplot(3,1,1)
plot(indicePicos_l,alturaPicos_l,'>');
xlabel('amostra');
ylabel('cm/s')
title('CBFv_l - Beat to Beat')
axis tight
hold on 
plot(CBFv_l);

subplot(3,1,2)
plot(indicePicos_r,alturaPicos_r,'>');
xlabel('amostra');
ylabel('cm/s')
title('CBFv_r - Beat to Beat')
axis tight
hold on 
plot(CBFv_r);

subplot(3,1,3)
plot(indicePicos_ABP,alturaPicos_ABP,'>');
xlabel('amostra');
ylabel('mmHg')
title('ABP - Beat to Beat')
axis tight
hold on 
plot(ABP_v);

%keyboard

% Media do Sinal CBFv Left & Interpolation Spline CBFv Left


figure(4)
subplot(2,1,1)
plot (Matriz_final(:,1)/100,Matriz_final (:, 2))
xlabel('tempo (segundos)');
ylabel('cm/s')
title('Media do Sinal CBFv Left')
axis tight

hold on

subplot(2,1,2)
plot(x_l/100,y_l,'+',xx_l,yy_l);
xlabel('tempo (segundos)');
ylabel('cm/s')
title('Media do Sinal CBFv Left - Interpolação SPLINE')
axis tight

%keyboard

% Media do Sinal CBFv Right & Interpolation Spline CBFv RIght

figure(5)
subplot(2,1,1)
plot (Matriz_final_r(:,1)/100, Matriz_final_r (:, 2))
xlabel('tempo (segundos)');
ylabel('cm/s')
title('Media do Sinal CBFv Right')
axis tight

hold on

subplot(2,1,2)
plot(x_r/100,y_r,'+',xx_r,yy_r);
xlabel('tempo (segundos)');
ylabel('cm/s')
title('Media do Sinal CBFv Right - Interpolação SPLINE')
axis tight

%keyboard

% Media do Sinal ABP & Interpolation Spline ABP

figure(6)
subplot(2,1,1)
plot (Matriz_final_abp(:,1)/100, Matriz_final_abp (:, 2))
xlabel('tempo (segundos)');
ylabel('mmHg')
title('Media do Sinal ABP')
axis tight

hold on

subplot(2,1,2)
plot(x_abp/100,y_abp,'+');
hold on;
plot(xx_abp,yy_abp);
xlabel('tempo (segundos)');
ylabel('cm/s')
title('Media ABP - Interpolação SPLINE')
axis tight

%keyboard

% Análise Espectral pós Welch

figure(7)
subplot(2,1,1)
plot([1:figLenghth]*fs_2/(2*figLenghth),abs(X(1:figLenghth)))
xlabel('Frequência(Hz)');
ylabel('MCAv Power (cm/s^2)/Hz')
title('Spectral Analysis')


subplot(2,1,2)
plot([1:figLenghth]*fs_2/(2*figLenghth),abs(Y(1:figLenghth)))
xlabel('Frequência(Hz)');
ylabel('MAP Power (mmHg^2)/Hz')
title('Spectral Analysis')
axis tight
hold off

%keyboard

%% Ganho, Fase, Coerência 

params.plot=1;
params.plot_f_range=[0,0.5];
params.plot_title='';

subplot1=[.15,0.65,0.75,0.22];
subplot2=subplot1+[0,-0.25,0,0];
subplot3=subplot2+[0,-0.25,0,0];

% keyboard

if params.plot
    figure (8)
    % set(fig,'Position',pos3);
    t=[0:length(yy_l_2)-1]/fs;
    plot(t,yy_abp_2,t,yy_l_2,'r:');
    title(params.plot_title);
    xlabel('time (s)');
    legend('ABP','CBFV');
    axis tight
    
    %keyboard
    
    figure(9);
    ax(1)=subplot('position',subplot1);%(3,1,1);
    plot(f,abs(G));
    title(params.plot_title);
    ylabel('Gain');   
    
    
    hold on
    plot(params.vlf,[1,1]*tfa_out.Gain_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Gain_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Gain_hf,':r');
    set(ax(1),'XTickLabel','');
    axis tight
    
    ax(2)=subplot('position',subplot2);%(3,1,2);
    set(ax(2),'XTickLabel','');
    plot(f,angle(G)/(2*pi)*360);
    hold on
    ylabel('Phase (deg)');
    plot(params.vlf,[1,1]*tfa_out.Phase_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Phase_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Phase_hf,':r');
    set(ax(2),'XTickLabel','');
    axis tight    
    
    ax(3)=subplot('position',subplot3);%(3,1,3);
    plot(f,abs(C).^2);
    ylabel('|Coh|^2');
    
    hold on
    
    plot(params.vlf,[1,1]*tfa_out.Coh2_vlf,':r');
    plot(params.lf,[1,1]*tfa_out.Coh2_lf,':r');
    plot(params.hf,[1,1]*tfa_out.Coh2_hf,':r');
    plot(params.plot_f_range,coherence2_threshold*ones(1,2),'--k');
    axis tight
    
    xlabel('frequency(Hz)');
    linkaxes(ax,'x');
    xlim(params.plot_f_range);
    
end

keyboard
%save P05_AVC.mat

%% ARI





% Encontra os pontos de inicio de diastole

function [PicoPos, PicoVal] = beattobeat (sinal)

distMin= 70;
PicoPos = 1;
PicoVal = sinal(1);
aPicoPos = 1;
aPicoVal = sinal(1);
j=0;
for i=2:length(sinal)-1
  
  if i>aPicoPos+distMin
    j=j+1;
    PicoPos(j) = aPicoPos; %local no eixo x dos picos
    PicoVal(j) = aPicoVal;
    pPicoPos = i;
    pPicoVal = sinal(i);
    aPicoPos = i;
    aPicoVal = sinal(i);
  end
  
  if sinal(i)<=sinal(i-1) && sinal(i)<=sinal(i+1)
    pPicoPos = i;
    pPicoVal = sinal(i);
    if pPicoVal < aPicoVal
      aPicoVal= pPicoVal;
      aPicoPos = pPicoPos;
    end
  end    
end

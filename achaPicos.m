% A função encontra os picos 

function [PicoPos,PicoVal] = achaPicos(sinal, distMin, tipo)
    PicoPos = 1;
    PicoVal = sinal(1);
    aPicoPos = 1;
    aPicoVal = sinal(1);
    j=0;
    for i=2:length(sinal)-1

      if i>aPicoPos+distMin
        j=j+1;
        PicoPos(j) = aPicoPos;
        PicoVal(j) = aPicoVal;
        pPicoPos = i;
        pPicoVal = sinal(i);
        aPicoPos = i;
        aPicoVal = sinal(i);
      end
      
      if strcmp(tipo,'positivo')==1
        if sinal(i)>sinal(i-1) && sinal(i)>sinal(i+1)
          pPicoPos = i;
          pPicoVal = sinal(i);
          if pPicoVal > aPicoVal
            aPicoVal = pPicoVal;
            aPicoPos = pPicoPos;
          end
        end
      end
      
      if strcmp(tipo,'negativo')==1
        if sinal(i)<=sinal(i-1) && sinal(i)<=sinal(i+1)
          pPicoPos = i;
          pPicoVal = sinal(i);
          if pPicoVal < aPicoVal
            aPicoVal = pPicoVal;
            aPicoPos = pPicoPos;
          end
        end
      end
      
    end
end
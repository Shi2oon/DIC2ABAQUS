function Pla_reg = isplastic(el,yield)

Pla_reg = (abs(mean(el.gpsmises,2))>yield);
Pla_reg =  find(Pla_reg~=0);


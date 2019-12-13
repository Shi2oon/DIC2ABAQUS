function [datum, unit, offset] = unist4Abaqus(datum,unit)
%% Units + 
switch unit
    case 'mm';          datum = datum.*1e-3;                 % conver to m
    case 'um';          datum = datum.*1e-6;                 % conver to m
end

% and for some reason Abaqus is really bad in handling small dim.
if      mean(abs(datum(:,1)))<5e-5;     offset = 1e-6;       unit = 'um';
elseif  mean(abs(datum(:,1)))<5e-2;     offset = 1e-3;       unit = 'mm'; 
else;   offset = 1;                     unit = 'm';                end

datum = datum./offset;
% switch unit
%     case 'mm';      datum = datum./offset;
%     case 'um';      datum = [datum(:,1:2)./offset datum(:,3:end).*1e-3; ];  
%         fprintf('X and Y dim is in mum but displacement is in mm\n');
% end

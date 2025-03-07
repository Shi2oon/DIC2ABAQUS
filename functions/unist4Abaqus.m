function [datum, unit, offset] = unist4Abaqus(datum,unit)
%% Units +
% if ~isempty(datum)
%     switch unit
%         case 'mm'   % conver to m
%             datum.X = datum.X*1e-3;             datum.Y = datum.Y*1e-3;
%             datum.Ux = datum.Ux*1e-3;             datum.Uy = datum.Uy*1e-3;
%         case 'um'	% conver to m
%             datum.X = datum.X*1e-6;             datum.Y = datum.Y*1e-6;
%             datum.Ux = datum.Ux*1e-6;             datum.Uy = datum.Uy*1e-6;
%     end
%     unit = 'm';
%     
%     % and for some reason Abaqus is really bad in handling small dim.
%     if      mean(abs(datum.X(:)))<6e-5;     offset = 1e-6;       unit = 'um';
%     elseif  mean(abs(datum.X(:)))<6e-2;     offset = 1e-3;       unit = 'mm';
%     else;   offset = 1;                      unit = 'm';                end
%     
%     datum.X = datum.X./offset;      datum.Y = datum.Y./offset;
%     datum.Ux = datum.Ux./offset;    datum.Uy = datum.Uy./offset;
% else
    switch unit
        case 'nm'	% conver to m
            offset = 1e-9;
        case 'um'	% conver to m
            offset = 1e-6;
        case 'mm'   % conver to m
            offset = 1e-3;
        case 'm'
            offset = 1;
    end
% end
% % switch unit
% %     case 'mm';      datum = datum./offset;
% %     case 'um';      datum = [datum(:,1:2)./offset datum(:,3:end).*1e-3; ];
% %         fprintf('X and Y dim is in mum but displacement is in mm\n');
% end

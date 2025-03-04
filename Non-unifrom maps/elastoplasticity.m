function [Yield_offset,Exponent] = elastoplasticity(E,yield,n)
% Define material properties
% E = 69.99e9; % Elastic Modulus in Pa
% yield = 193e6; % Yield Stress in Pa
% epsilon_y = 0.002;
% n = 0.13; % Strain Hardening Exponent (assumed)
epsilon_y =  yield/ E; % Yield Strain

% Define strain range
epsilon = linspace(0, 0.05, 1000); % Strain from 0 to 5%

% Calculate stress
sigma = zeros(size(epsilon)); nn=0;
for i = 1:length(epsilon)
    if epsilon(i) <= epsilon_y
        sigma(i) = E * epsilon(i); % Elastic region
    else
        
%         if nn == 1
%             StartPlasticity = i;
%         end
        sigma(i) = yield * (epsilon(i) / epsilon_y)^n; % Plastic region
        if yield <= sigma(i)
           nn = nn+1; 
        sigmaPl(nn) = sigma(i);
        epsilonPl(nn) = epsilon(i);
        end
    end
end

% Plot stress-strain curve
% figure;
% plot(epsilon, sigma / 1e6, 'LineWidth', 2); % Stress in MPa
% xlabel('Strain');
% ylabel('Stress (MPa)');
% title('Stress-Strain Curve for 5052 Aluminum Alloy');
% grid on;

% RO_yield_offset = 1-(epsilon_y*/yield/E);

[~, ~,Yield_offset,Exponent] = ...
    fitforRamberg_Osgood_relationship(sigmaPl, epsilonPl*E,yield);


end
clc;clear
% Define material properties
E = 70.3e9; % Elastic Modulus in Pa
sigma_y = 193e6; % Yield Stress in Pa
nu = 0.33; % Poisson's Ratio
n = 0.2; % Strain Hardening Exponent (assumed)
epsilon_y = sigma_y / E; % Yield Strain
% Define strain range
epsilon = linspace(0, 0.05, 500); % Strain from 0 to 5%
% Calculate stress
sigma = zeros(size(epsilon)); nn=0;
for i = 1:length(epsilon)
    if epsilon(i) <= epsilon_y
        sigma(i) = E * epsilon(i); % Elastic region
    else
        nn = nn+1;
        if nn == 1
            StartPlasticity = i;
        end
        sigma(i) = sigma_y * (epsilon(i) / epsilon_y)^n; % Plastic region
        MatrixAbaqus(nn,:) = [sigma(i) epsilon(i)-epsilon(StartPlasticity)];
    end
end
% Plot stress-strain curve
figure;
plot(epsilon, sigma / 1e6, 'LineWidth', 2); % Stress in MPa
xlabel('Strain');   ylabel('Stress (MPa)');
title('Stress-Strain Curve for 5052 Aluminum Alloy');   grid on;
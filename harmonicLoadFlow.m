clc
clear all;
% Harmonic Load Flow Analysis for IEEE 33-bus radial system
% System parameters
num_buses = 33;
num_harmonics = 5;  % Including fundamental
V_base = 12.66e3;  % Base voltage (V)
S_base = 100e6;  % Base power (VA)
Z_base = V_base^2 / S_base;  % Base impedance (ohm)

% Load IEEE 33-bus system data
[bus_data, branch_data] = load_ieee33_data();

% Initialize matrices
V = ones(num_buses, 1) * [1, 0.01, 0.005, 0.003, 0.002];
I = zeros(num_buses, num_harmonics);
S = zeros(num_buses, num_harmonics);

% Iterative solution
max_iterations = 100;
tolerance = 1e-6;
for iteration = 1:max_iterations
    V_prev = V;
    
    % Solve for each harmonic
    for h = 1:num_harmonics
        % Forward sweep (calculate currents)
        I(:,h) = calculate_currents(bus_data, V(:,h), h);
        
        % Backward sweep (update voltages)
        V(:,h) = update_voltages(branch_data, I(:,h), h);
        
        % Calculate power
        S(:,h) = V(:,h) .* conj(I(:,h));
    end
    
    % Check convergence
    if max(abs(V - V_prev)) < tolerance
        break;
    end
end

% Calculate and display results
display_results(V, I, S, bus_data);

% Helper functions
function [bus_data, branch_data] = load_ieee33_data()
    % Load IEEE 33-bus system data
    % This is a placeholder. You need to implement this function to load actual data.
    V_base = 12.66e3;  % Base voltage (V)
    S_base = 100e6;  % Base power (VA)
    Z_base = V_base^2 / S_base;
    
    raw_bus_data = load('loaddata33bus.m')  % [bus_number, P_load, Q_load]
    raw_branch_data = load('linedata33bus.m') % [from_bus, to_bus, R, X]
    
    % Fill these with actual IEEE 33-bus system data
    S_base = 100; % Assuming 100 MVA base, adjust if different
    bus_data = raw_bus_data;
    bus_data(:, 2) = bus_data(:, 2) / (1000 * S_base); % Convert kW to p.u.
    bus_data(:, 3) = bus_data(:, 3) / (1000 * S_base); % Convert kVAr to p.u.
    
    branch_data = raw_branch_data;
    branch_data(:, 3) = raw_branch_data(:, 3) / Z_base; % Convert R to p.u.
    branch_data(:, 4) = raw_branch_data(:, 4) / Z_base; % Convert X to p.u.
end

function I = calculate_currents(bus_data, V, h)
    % Calculate currents for each bus
    S = bus_data(:, 2) + 1j * bus_data(:, 3);  % Complex power
    I = conj(S ./ V);
    
    % Add harmonic current injections
   if h > 1
    I = I + (0.01/h) * abs(I) .* exp(1j * rand(size(I)) * 2 * pi);
end
end

function V = update_voltages(branch_data, I, h)
    % Backward sweep to update voltages
    V = ones(size(I));  % Start with all voltages at 1.0 pu
    for b = size(branch_data, 1):-1:1
        from_bus = branch_data(b, 1);
        to_bus = branch_data(b, 2);
        Z = branch_data(b, 3) + 1j * branch_data(b, 4);
        
        % Scale impedance for harmonics (Z is already in p.u.)
        Z = Z * h;
        
        V(from_bus) = V(to_bus) + Z * I(to_bus);
    end
end

function display_results(V, I, S, bus_data)
    % Display voltage profile
    figure;
    plot(abs(V(:,1)));
    hold on;
    for h = 2:size(V,2)
        plot(abs(V(:,h)), '--');
    end
    hold off;
    title('Voltage Magnitude Profile');
    xlabel('Bus Number');
    ylabel('Voltage (p.u.)');
    legend('Fundamental', '2nd', '3rd', '4th', '5th');

    % Display harmonic distortion
    disp('Voltage magnitudes:');
disp(abs(V));
    THD = sqrt(sum(abs(V(:,2:end)).^2, 2)) ./ abs(V(:,1)) * 100;
    figure;
    bar(THD);
    title('Total Harmonic Distortion (THD) of Voltage');
    xlabel('Bus Number');
    ylabel('THD (%)');
end
% postprocessing comsol sweeps
% Copyright (C) 2024 Minke Berghuis
% All rights reserved.
%
% This script/function is proprietary software and is provided "as is" without warranty
% of any kind, express or implied, including but not limited to the warranties of
% merchantability, fitness for a particular purpose, and non-infringement.
%
% Redistribution and use in source and binary forms, with or without modification, are
% permitted only with prior written permission from the author.
%
% Author: Minke Berghuis (m.w.berghuis@utwente.nl)
% Date: 2025-02-14
% Version: 1.0
%
% Description:
% loads Comsol result data from .xlsx files and generates
% result/supplementary figures. 

clc, clear, close all

addpath(pwd); % to find the required .m functions
dataFolder = strcat(pwd, '\axisymsensor_results');

% old data (oil 0.1, water, air)
dataFileName1 = '202501151726_mu_0_1_Fsweep_results.xlsx';
dataFile1 = strcat(dataFolder, '\',dataFileName1);
dataFileName2 = '202501151726_mu_0_001_Fsweep_results.xlsx';
dataFile2 = strcat(dataFolder, '\',dataFileName2);

dataFileName3 = '202501161852_mu_0_1_Fsweep_interior_wall_results.xlsx';
dataFile3 = strcat(dataFolder, '\',dataFileName3);
dataFileName4 = '202501161850_mu_0_001_Fsweep_interior_wall_results.xlsx';
dataFile4 = strcat(dataFolder, '\',dataFileName4);

dataFileName5 = '202501151726_mu_air_Fsweep_results.xlsx';
dataFile5 = strcat(dataFolder, '\',dataFileName5);

dataFiles = {dataFile1 dataFile2 dataFile3 dataFile4 dataFile5};
S_old = importResultExcel(dataFiles);

% x/L different indentation location (same oil 0.01, same load 40k)
dataFileName1 = 'mu_0_01_xLsweep_0_01_results.xlsx';
dataFile1 = strcat(dataFolder, '\',dataFileName1);
dataFileName2 = 'mu_0_01_xLsweep_0_05_results.xlsx';
dataFile2 = strcat(dataFolder, '\',dataFileName2);
dataFileName3 = 'mu_0_01_xLsweep_0_09_results.xlsx';
dataFile3 = strcat(dataFolder, '\',dataFileName3);
dataFiles = {dataFile1 dataFile2 dataFile3};
S_xL = importResultExcel(dataFiles);

% oil sweeps: new loadprofile! 
dataFileName1 = '202502060934_mu_0_1_Fsweep_results.xlsx';
dataFile1 = strcat(dataFolder, '\',dataFileName1);
dataFileName2 = '202502060934_mu_0_2_Fsweep_results.xlsx';
dataFile2 = strcat(dataFolder, '\',dataFileName2);
dataFileName3 = '202502060934_mu_1_0_Fsweep_results.xlsx';
dataFile3 = strcat(dataFolder, '\',dataFileName3);
dataFiles = {dataFile1 dataFile2 dataFile3};
S_oil = importResultExcel(dataFiles);

% interior wall
dataFileName1 = '202501161852_mu_0_1_Fsweep_interior_wall_results.xlsx';
dataFile1 = strcat(dataFolder, '\',dataFileName1);
dataFileName2 = '202501161852_mu_0_2_Fsweep_interior_wall_results.xlsx';
dataFile2 = strcat(dataFolder, '\',dataFileName2);
dataFileName3 = '202501161852_mu_1_0_Fsweep_interior_wall_results.xlsx';
dataFile3 = strcat(dataFolder, '\',dataFileName3);
dataFileName4 = '202502081640_mu_air_Fsweep_interior_wall_results.xlsx';
dataFile4 = strcat(dataFolder, '\',dataFileName4);
% dataFileName4 = '202501161850_mu_0_001_Fsweep_interior_wall_results.xlsx';
% dataFile4 = strcat(dataFolder, '\',dataFileName4);

dataFiles = {dataFile1 dataFile2 dataFile3 dataFile4};
S_wall = importResultExcel(dataFiles);

%% prepare for plotting
[lineStyles, plotColors, markerStyles] = plotSettingsSymsensor();
% oil color shades
N = 3;  % Number of shades
oilColor = plotColors{1};
%oilShades = interp1([0 0.5 1], [0 0 0; oilColor; 1 1 1], linspace(0.3,0.7,N)); % dark to light
oilShades = interp1([0 0.5 1], [1 1 1; oilColor; 0 0 0], linspace(0.3,0.7,N)); % light to dark
% water color shades
N = 5;  % Number of shades
waterColor = plotColors{2};
waterShades = interp1([0 0.5 1], [0 0 0; waterColor; 1 1 1], linspace(0.3,0.7,N));
% air color shades
N = 5;  % Number of shades
airColor = plotColors{3};
airShades = interp1([0 0.5 1], [0 0 0; airColor; 1 1 1], linspace(0.3,0.7,N));

%% load profiles
figure(1), clf(1), hold on
xlabel('$t$ (s)')
ylabel('$F_{z}$ (kPa)')
fields = fieldnames(S_oil); % Get all field names
N = numel(fields);      % Get the number of fields
names = {'10 kPa','15 kPa','20 kPa','25 kPa','30 kPa','35 kPa','40 kPa','45 kPa','50 kPa'};
F_max_idx = zeros(1,N); % store the index for reaching max load

for i = 1:N-1 % neglect empty last field
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(1).(field_name); % Access the field value
    if ~isempty(field_value)
        [~, F_max_idx(i)] = max(round(field_value.loadprofile)); % `max` returns first occurrence of max value, round for num.noise
        x = field_value.Time_s_; % s
        y = 10^(-3).*field_value.loadprofile; % kPa
        % instead of peak load, integrate the gaussian profile to obtain
        % total load per crosssect (kN/m)
            % \int (exp(-a(x+b)^2)) dx = \sqrt(pi/a) --> a = 1/(2*\sigma^2),
        sigma = 0.008*6;    % gaussian width = 6*\sigma = d_intenter = 8mm
        y = y.*sqrt(pi*2*sigma^2);
        plot(x,y,'Color',[1 0 0],'DisplayName',field_name)
        text(x(end)-0.013, y(end)-0.25, names{i}, 'FontSize', 12, 'FontWeight', 'bold');
        % plot the point for reaching max load
        p1 = plot(x(F_max_idx(i)),y(F_max_idx(i)),'*','MarkerEdgeColor',[1 0 0])
        % plot the point for reaching max load
        p2 = plot(x(F_max_idx(i)+30),y(F_max_idx(i)+30),'o','MarkerEdgeColor',[1 0 0])
    end
end
legend([p1 p2], {'max load', 'max load + $\Delta t$'}, 'Location', 'NorthWest')

Title = 'load profiles';

%% Oil delta_p over time for different loadcases
figure(2), clf(2), hold on
xlabel('$t$ (s)')
ylabel('$\Delta p$ (kPa)')
fields = fieldnames(S_oil); % Get all field names
N = numel(fields);      % Get the number of fields
N = N-1; % the last field is rubbish
names = {'100 cSt','200 cSt','1000 cSt'};

for j = 1:3 % 3 oil viscosities
for i = 3 % only plot the 40kPa loadcase
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(j).(field_name); % Access the field value
    if ~isempty(field_value)
        x = field_value.Time_s_; % s
        y = 10^(-3).*(field_value.delta_p_Pa_); % kPa
        plot(x,y,'Color',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load
        plot(x(F_max_idx(i)),y(F_max_idx(i)),'*','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load +30*dt
        plot(x(F_max_idx(i)+30),y(F_max_idx(i)+30),'o','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
    end
end 
p(j)=plot(NaN, NaN, 'Color',oilShades(j,:),'DisplayName',names{j}); % for the legend
end
legend(p, 'Location', 'NorthWest')

Title = 'Oil 0_1 delta_p over time for different loadcases';


%% Oil gap over time for different loadcases
figure(3), clf(3), hold on
xlabel('$t$ (s)')
ylabel('$gap$ (mm)')
fields = fieldnames(S_oil); % Get all field names
N = numel(fields);      % Get the number of fields
N = N-1; % the last field is rubbish
names = {'100 cSt','200 cSt','1000 cSt'};

for j = 1:3 % 3 oil viscosities
for i = 3 % 1:N: all loadcases
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(j).(field_name); % Access the field value
    if ~isempty(field_value)
        x = field_value.Time_s_; % s
        y = 10^(3).*(field_value.gap_m_); % mm
        plot(x,y,'Color',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load
        plot(x(F_max_idx(i)),y(F_max_idx(i)),'*','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load +30*dt
        plot(x(F_max_idx(i)+30),y(F_max_idx(i)+30),'o','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
    end
end 
p(j)=plot(NaN, NaN, 'Color',oilShades(j,:),'DisplayName',names{j}); % for the legend
end
offset = 0.1; % mm
p2 = yline(offset,'-k', 'DisplayName', 'offset');
legend([p p2], 'Location', 'NorthEast')
axis([0 0.1 0 2])
Title = 'Oil gap over time for different loadcases';
Title = 'Oil gap over time for 20kPa';

%% Oil delta_z over time for different loadcases
figure(4), clf(4), hold on
xlabel('$t$ (s)')
ylabel('$\Delta z$ (mm)')
fields = fieldnames(S_oil); % Get all field names
N = numel(fields);      % Get the number of fields
N = N-1; % the last field is rubbish
names = {'100 cSt','200 cSt','1000 cSt'};

for j = 1:3 % 3 oil viscosities
for i = 3 % 1:N all loadcases
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(j).(field_name); % Access the field value
    if ~isempty(field_value)
        x = field_value.Time_s_; % s
        y = 10^(3).*(0.005-field_value.delta_z_m_); % mm, initial outer radius = 5mm
        plot(x,y,'Color',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load
        plot(x(F_max_idx(i)),y(F_max_idx(i)),'*','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load +30*dt
        plot(x(F_max_idx(i)+30),y(F_max_idx(i)+30),'o','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
    end
end 
p(j)=plot(NaN, NaN, 'Color',oilShades(j,:),'DisplayName',names{j}); % for the legend
end
legend([p], 'Location', 'NorthWest')
axis([0 0.1 0 2])
Title = 'Oil delta z over time for different loadcases';
Title = 'Oil delta z over time for 20kPa';

%% Oil, air, interior wall delta_p for different loadcases at fixed t
% t_fmax and t_fmax+30*dt
figure(5), clf(5), hold on
xlabel('$\Delta z$ (mm)')
ylabel('$\Delta p$ (kPa)')

% Oil with \mu = 0.1, 0.2, 1.0 Pas
fields = fieldnames(S_oil); % Get all field names
N = numel(fields)-1;      % Get the number of fields, last is rubbish
names = {'100 cSt','200 cSt','1000 cSt'};
for j = 1:3
% time of max load
dp = nan(1,N);
dz = nan(1,N);
% time of max load + dt
dp1 = nan(1,N);
dz1 = nan(1,N);
for i = 1:N
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(j).(field_name); % Access the field value
    if ~isempty(field_value)
        dp(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i))); % kPa
        %dz(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i))); % gap, mm
        dz(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i))); % delta_z, mm
        dp1(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i)+30)); %kPa
        %dz1(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i)+30)); % gap, mm
        dz1(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i)+30)); % delta_z, mm
    end
end 
p(j)=plot(dz,dp,'-*','Color',oilShades(j,:),'DisplayName',names{j})
% plot(dz1,dp1,'-o','Color',oilShades(j,:),'DisplayName',names{j})
end

% air (old data)
fields = fieldnames(S_old); % Get all field names
N = numel(fields)-1;      % Get the number of fields, last is rubbish
names = {'air'};
for j = 5 % air 
% time of max load
dp = nan(1,N);
dz = nan(1,N);
% time of max load + dt
dp1 = nan(1,N);
dz1 = nan(1,N);
for i = 1:N
    field_name = fields{i};  % Get the field name as a string
    field_value = S_old(j).(field_name); % Access the field value
    if ~isempty(field_value)
        dp(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i))); % kPa
        %dz(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i))); % gap, mm
        dz(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i))); % delta_z, mm
        dp1(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i)+30)); %kPa
        %dz1(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i)+30)); % gap, mm
        dz1(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i)+30)); % delta_z, mm
    end
end 
p1=plot(dz,dp,'-*','Color',airShades(3,:),'DisplayName',names{1})
% plot(dz1,dp1,'-o','Color',oilShades(j,:),'DisplayName',names{j})
end

% interior wall
fields = fieldnames(S_wall); % Get all field names
N = numel(fields)-1;      % Get the number of fields, last is rubbish
names = {'100 cSt, wall','200 cSt, wall','1000 cSt, wall','air, wall'};
Shades = [oilShades; airShades(3,:)];
for j = 1:length(names)
dp = nan(1,N);
dz = nan(1,N);
dp1 = nan(1,N);
dz1 = nan(1,N);
for i = 1:N
    field_name = fields{i};  % Get the field name as a string
    field_value = S_wall(j).(field_name); % Access the field value
    if ~isempty(field_value)
        dp(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i)));
        %dz(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i)));
        dz(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i)));
        dp1(i) = 10^(-3).*(field_value.delta_p_Pa_(F_max_idx(i)+30));
        %dz1(i) = 10^(3).*(0.002-field_value.gap_m_(F_max_idx(i)+30));
        dz1(i) = 10^(3).*(0.005-field_value.delta_z_m_(F_max_idx(i)+30));
    end
end 
p2(j) = plot(dz,dp,'-.*','Color',Shades(j,:),'DisplayName',names{j})
% plot(dz1,dp1,'-.o','Color',Shades(j,:),'DisplayName',names{j})
end

dx = 0.1;
axis([-dx 2 -dx 5])
legend([p p1 p2], 'Location', 'NorthWest')

Title = 'Oil air interior wall delta_p for different loadcases at fixed t';

%% x/L different indentation location (same oil 0.01, same load 40k)
figure(6), clf(6), hold on

fields = fieldnames(S_xL); % Get all field names
N = numel(fields);      % Get the number of fields
names = {'x/L = 0.1','x/L = 0.5','x/L = 0.9'};
t = tiledlayout(1,3); % Creates a layout for 1 row and 3 columns
for j = 1:3
    nexttile;
    axis([0 1 -1 4])
    if j ==1
        ylabel('$\Delta p$ (kPa)')
    end
    hold on
    for i = 1:N
    field_name = fields{i};  % Get the field name as a string
    field_value = S_xL(j).(field_name); % Access the field value
    x = field_value.Time_s_; % s
    y = 10^(-3).*(field_value.p_top_Pa_); % kPa
    y2 = 10^(-3).*(field_value.p_bottom_Pa_); % kPa
    if ~isempty(field_value)
        plot(x,y,'Color',plotColors{4})
        plot(x,y2,'Color',plotColors{5})
    end
    end 
    title(names{j});
end
xlabel(t, '$t$ (s)','interpreter','latex'); % Single xlabel

Title = 'x_L different indentation location -oil 0.01 load 40k)';

%% flux for different mu
figure(7), clf(7), hold on
xlabel('$t$ (s)')
ylabel('$Q/W$ (mm$^2$/s)')
fields = fieldnames(S_oil); % Get all field names
N = numel(fields);      % Get the number of fields
N = N-1; % the last field is rubbish
names = {'100 cSt','200 cSt','1000 cSt'};

% the gap for the different \mu under the same loadcase is different! (see gap plot)
% still, clear that flux scales with inversely with \mu 
for j = 1:3 % 3 oil viscosities
for i = 3 % 1 = load case: 10kPa 3= load case 20kPa
    field_name = fields{i};  % Get the field name as a string
    field_value = S_oil(j).(field_name); % Access the field value
    x = field_value.Time_s_; % s
    y = -10^(6).*(field_value.flux_m_3_s_); % mm2/s, change to positive
    if ~isempty(field_value)
        plot(x,y,'Color',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load
        plot(x(F_max_idx(i)),y(F_max_idx(i)),'*','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
        % plot the piont for reaching max load +30*dt
        plot(x(F_max_idx(i)+30),y(F_max_idx(i)+30),'o','MarkerEdgeColor',oilShades(j,:),'DisplayName',field_name)
    end
end 
p(j)=plot(NaN, NaN, 'Color',oilShades(j,:),'DisplayName',names{j}); % for the legend
end
legend(p, 'Location', 'NorthEast')

Title = 'flux over time for different oil viscosities - load 20kPa';



%% Local functions
function S = importResultExcel(dataFiles)
% input: cell array of N datafile locations, with unique column headers for
% Nv variable names. The first columns contains the Np unique parameter
% values of the parameter sweep
% output: struct S containing the N result structs with Np parameter fields with Nv variable fields (automatic naming)
S = struct();  % Initialize an empty struct

for ii = 1:length(dataFiles)
opts = detectImportOptions(dataFiles{ii});
opts = setvartype(opts, 'double');  % Force all columns to double
data = readtable(dataFiles{ii}, opts);
numericData = table2array(data);

varNames = data.Properties.VariableNames;

% Extract unique parameter values from the first column
uniqueParams = unique(data{:,1});  % Assuming the first column contains unique identifiers

for i = 1:length(uniqueParams)
    paramValue = uniqueParams(i);  % Get the current parameter value
    
    % Find rows corresponding to this parameter value
    rows = data{:,1} == paramValue;
    
    % Extract remaining columns based on variable names
    for col = 2:width(data)
        varName = data.Properties.VariableNames{col};  % Get the variable name
        S(ii).(sprintf('P%d', paramValue)).(varName) = data{rows, col};
    end
end
end
end
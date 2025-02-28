function [lineStyles, plotColors, markerStyles] = plotSettingsArman(varargin)
% default settings
defaultInterpreter = 'latex';
defaultFontSize = 14;
defaultFontName = 'Helvetica';
defaultLineWidth = 1;
defaultFigureColor = 'w';
defaultPlotLineWidth = 1;
defaultlineStyles = {'-', '--', ':', '-.'};
defaultmarkerStyles = {'o', '+', '*', 'x', 'square', '^', 'v'};

    %[1 0 0], % red 
    %[0 1 0], % Green, 
    %[0 0 1], % Blue,
    %[0.9290 0.6940 0.1250], % yellow
    %[0 0 0]}; % Black

    % using the Colorblind color palette; check figure
    %     figure
    %     imageData = imread('C:\Users\BerghuisMW\surfdrive\Documents\MATLAB\myfunctions\plotting\colorblind_colomap_picture.jpg');
    %     imshow(imageData);

    S = load('colorblind_colormap.mat'); 
    CT = S.colorblind; % color table in matrix format
    plotColors = mat2cell(CT, ones(1, size(CT, 1)), size(CT, 2)); % cell format
    % Indices of the entries to swap
    idx1 = [1 2 3 4 5 6 7 8 9 10 11 12];  % original
    idx2 = [1 12 2 6 7 10 3 4 5 8 9 11];  % new (drag and drop)
    plotColors = plotColors(idx2); 

    % overwrite with (normalized) custom colors
    plotColors(1) = {[251 176 66]/255}; % oil 
    plotColors(2) = {[34 117 188]/255}; % water
    plotColors(3) = {[128 128 128]/255}; % air?
    plotColors(4) = {[89 187 223]/255}; % blue-ish
    plotColors(5) = {[237 91 93]/255}; % reddish 

defaultplotColors = plotColors;

% parse the input  
p = inputParser;
addParameter(p,'Interpreter',defaultInterpreter,@(x) ischar(x) || isstring(x));
addParameter(p, 'FontSize', defaultFontSize, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FontName', defaultFontName, @(x) ischar(x) || isstring(x));
addParameter(p, 'LineWidth', defaultLineWidth, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'FigureColor', defaultFigureColor, @(x) ischar(x) || isstring(x));
addParameter(p, 'PlotLineWidth', defaultPlotLineWidth, @(x) isnumeric(x) && isscalar(x));
addParameter(p,'lineStyles',defaultlineStyles,@(x) iscell(x));
addParameter(p,'markerStyles',defaultmarkerStyles,@(x) iscell(x));
addParameter(p,'plotColors',defaultplotColors,@(x) iscell(x));

% Handle reset condition
if nargin == 1 && strcmpi(varargin{1}, 'reset')
    resetDefaults();
    disp('Plot settings reset to default.');
    return;
end

% Parse input arguments
try
    parse(p, varargin{:});
catch ME
    error('Error parsing inputs: %s', ME.message);
end

% return some plotsettings
lineStyles = p.Results.lineStyles;
markerStyles = p.Results.markerStyles;
plotColors = p.Results.plotColors;

% Set default figure properties
set(0, 'defaultTextInterpreter', p.Results.Interpreter);
set(0, 'defaultLegendInterpreter', p.Results.Interpreter);
set(0, 'defaultColorbarTickLabelInterpreter', p.Results.Interpreter);
set(0, 'defaultColorbarTickLabelInterpreter', p.Results.Interpreter);
set(0, 'defaultLegendFontSize', p.Results.FontSize);
set(0, 'DefaultAxesFontSize', p.Results.FontSize);                  
set(0, 'DefaultAxesFontName', p.Results.FontName);        
set(0, 'DefaultAxesLineWidth', p.Results.LineWidth);      
set(0, 'DefaultLineLineWidth', p.Results.PlotLineWidth);   
set(0, 'DefaultFigureColor', p.Results.FigureColor);  

% not yet added to the parser
set(0, 'defaultAxesBox', 'on');
% set(0, 'DefaultFigurePBaspect', [5 3 1]); %pbaspect([5 3 1 ])
set(0, 'DefaultAxesDataAspectRatio', [5 3 1]);
set(0, 'defaultLegendBox', 'off');
end


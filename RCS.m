% computes and plots the radar cross section for output from the embedded
% boundary wave solver

clear all;
close all;

r = 0.4;
t = 800;
n_theta = 1000;
rcs = zeros(n_theta);

% upload a .json file and check its contents
fname = 'input.json';
val = jsondecode(fileread(fname));

k = sqrt(val.kx^2 + val.ky^2);
wavelength = 2 * pi / k;

M = val.M; Mx = M; My = M; L = val.L; Lx = L; Ly = L;

U = zeros(Mx,My);
x_space = linspace(0,Lx,Mx);
y_space = linspace(0,Ly,My);

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 1);

% Specify range and delimiter
opts.DataLines = [4, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = "VecObject8MPIprocesses";
opts.VariableTypes = "double";

% Specify file level properties
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
%%

% Import the data
dat = readtable("/home/peter/git/embedded-boundary-method-wave-solver/Circle/wave"+t+".out", opts);

% Convert to output type
tmp = table2array(dat);

U(:,:) = reshape(tmp(1:  end,:),Mx,My);

%% Calculate the radar cross section
for ii = 1:n_theta
    theta = 2 * pi * ii / n_theta;
    
    % calculate the root mean square at that angle and radius
    n_rms_calc = 50;
    rms = 0;
    for jj = 1:n_rms_calc
        x = (r - (jj - 1) * wavelength / n_rms_calc) * sin(theta) + Lx / 2;
        y = (r - (jj - 1) * wavelength / n_rms_calc) * cos(theta) + Ly / 2;
        rms = rms + interp2(x_space, y_space, U, x, y)^2;
    end
    rms = sqrt(rms);
    
    E_s = rms;
    E_i = 1;
    
    rcs(ii) = r * E_s / E_i;
end

figure
polarplot(2 * pi - linspace(0, 2 * pi, n_theta), rcs)
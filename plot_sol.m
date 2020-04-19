clear all; 
close all;
Mx = 200; My = 200; Lx = 1.0; Ly = 1.0; Nt = 500;

U = zeros(Mx,My,Nt);

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
x = linspace(0,Lx,Mx);
y = linspace(0,Ly,My);

for i = 1:Nt
  % Import the data
  dat = readtable("/home/spencerwalker/Documents/repos/embedded-boundary-method-wave-solver/wave"+(i-1)+".out", opts);

  %% Convert to output type
  tmp = table2array(dat);

  U(:,:,i) = reshape(tmp(1:4:end,:),Mx,My);
  imagesc(x,y,((U(:,:,i))));
  colormap('jet')
  colorbar
  axis([0.15 1-0.15 0.15 1-0.15])
  caxis([-1,1])
  pause(0.001)
end 
%% Clear temporary variables
clear opts tmp



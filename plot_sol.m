Mx = 50; My = 50; Lx = 1.0; Ly = 1.0; Nt = 500;

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

for i = 1:Nt
  % Import the data
  dat = readtable("/home/spencerwalker/Documents/repos/embedded-boundary-method-wave-solver/wave"+(i-1)+".out", opts);

  %% Convert to output type
  tmp = table2array(dat);

  U(:,:,i) = reshape(tmp(1:2:end,:),Mx,My);
  imagesc(abs(U(:,:,i)).^2);
  pause(0.1)
end 
%% Clear temporary variables
clear opts tmp



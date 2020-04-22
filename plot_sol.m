clear all; 
close all;
Mx = 400; My = 400; Lx = 1.0; Ly = 1.0; Nt = 500;

U = zeros(Mx,My);

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
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'nextplot','replacechildren'); 
v = VideoWriter('test.avi');
open(v);
for i = 1:1:Nt
  % Import the data
  dat = readtable("/home/spencerwalker/Documents/repos/embedded-boundary-method-wave-solver/wave"+(i-1)+".out", opts);

  %% Convert to output type
  tmp = table2array(dat);

  U(:,:) = reshape(tmp(1:  end,:),Mx,My);
      imagesc(x,y,(U(:,:)))  ;
  colormap('jet')
  caxis([-1,1])
  colorbar
  hold on
  %r=0.1;
  %x0=Lx/2;
  %y0=Ly/2;
  %fplot(@(t) r*cos(t)+x0,@(t) r*sin(t)+y0,[0,2*pi],'k','LineWidth',2)
  xline(0.1)
  yline(0.1)
  xline(0.9)
  yline(0.9)
  xline(0.4)
  yline(0.4)
  xline(0.6)
  yline(0.6)
  hold off
  axis equal
  frame = getframe(gcf);
  writeVideo(v,frame);
end 
%% Clear temporary variables
clear opts tmp
%%
close(v)
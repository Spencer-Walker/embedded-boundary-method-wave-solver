clear all; 
close all;
M = 400;

Mx = M;
My = M; Lx = 1.0; Ly = 1.0; Nt = 1000;

U = zeros(Mx,My);

%%


x = linspace(0,Lx,Mx);
y = linspace(0,Ly,My);
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gca,'nextplot','replacechildren'); 
v = VideoWriter('test.avi');
open(v);
for i = 1:1:Nt
  % Import the data
  %% Initialize variables.
filename = "/home/becker/spwa4419/Documents/repos/embedded-boundary-method-wave-solver/wave"+(i-1)+".out";
delimiter = ' ';
startRow = 3;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts text in the input cell array to numbers. Replaced non-numeric
% text with NaN.
rawData = dataArray{1};
for row=1:size(rawData, 1)
    % Create a regular expression to detect and remove non-numeric prefixes and
    % suffixes.
    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
    try
        result = regexp(rawData(row), regexstr, 'names');
        numbers = result.numbers;
        
        % Detected commas in non-thousand locations.
        invalidThousandsSeparator = false;
        if numbers.contains(',')
            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                numbers = NaN;
                invalidThousandsSeparator = true;
            end
        end
        % Convert numeric text to numbers.
        if ~invalidThousandsSeparator
            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
            numericData(row, 1) = numbers{1};
            raw{row, 1} = numbers{1};
        end
    catch
        raw{row, 1} = rawData{row};
    end
end


%% Exclude rows with non-numeric cells
I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
raw(I,:) = [];

%% Allocate imported array to column variable names
tmp = cell2mat(raw(:, 1));

%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp R;
  %% Convert to output type

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
  %xline(0.1)
  %yline(0.1)
  %xline(0.9)
  %yline(0.9)
  %xline(0.4)
  %yline(0.4)
  %xline(0.6)
  %yline(0.6)
  hold off
  axis equal
  frame = getframe(gcf);
  writeVideo(v,frame);
end 
%% Clear temporary variables
clear opts tmp
%%
close(v)

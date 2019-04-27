% This code uses the Transmission Line Method (TLM) to estimate the sheet
% resistance and contact resistance of a sample.
% 
% For reference see:
% http://tuttle.merc.iastate.edu/ee432/topics/metals/tlm_measurements.pdf
% 
% The data is read using the following directory tree:
% dataPath
% |-- L1um
% |     |--datafile1.mat
% |     |--datafile2.mat
% |     |-- .
% |     |-- .
% |     |-- .
% |-- L2um
% |     |--datafile1.mat
% |     |--datafile2.mat
% |     |-- .
% |     |-- .
% |     |-- .
% |-- .
% |-- .
% |-- .
%
% The electrode distance is read from the name of the subfolder. The
% following convetion is used:
% LXum
% X is an integer greater than zero.
%
% The output folder is named after the dataPath:
% output_dataPath
%
% Author: erickrmartinez@gmail.com

clear all; close all;


% The directory containing the data
dataPath = '1V';
% Specify the contact width in cm
contactWidth    = 10E-4;   % cm
contactWidthErr = 1E-4;     % cm

% Read the contents of the dataPath
contents = dir(dataPath);
% Construct the output folder name
outputPath = strcat('output_',dataPath);
% If the folder does not exist, create it
if 7 ~= exist(outputPath,'dir')
    mkdir(outputPath);
end

% Current working directory
cwd = pwd;
% The regular expression used to extract the distance between electrodes
cExpression = '^L(\d+)um';
% Count the number of data subfolders. Do not count the current (.) and
% upper (..) folder pointers. 
nFolders = 0;
for i=1:length(contents)
    if contents(i).isdir && ~strcmp(contents(i).name,'.') && ~strcmp(contents(i).name,'..')
        nFolders = nFolders + 1;
    end
end

% Prealocate arrays for the results
TL      = zeros(length(nFolders),1); % The length in cm
R       = zeros(length(nFolders),1); % The resistance in Ohms
Rerr    = zeros(length(nFolders),1); % The uncertainty of the resistance


% Iterate over all the subfolders
cnt = 1; % Define a counter 
for i=1:length(contents)
    if contents(i).isdir && ~strcmp(contents(i).name,'.') ...
            && ~strcmp(contents(i).name,'..')
        % Get the name of the folder
        dataFolderName = contents(i).name;
        % Get the local path to look for .mat files
        localPath = strcat(dataPath,'/',dataFolderName,'/');
        % Find all .mat files in the path
        files = dir(strcat(localPath,'*.mat'));
        % Count the numbe of files
        nMatFiles = length(files);
        % Match the regular expression
        [mat, tok] = regexp(dataFolderName,cExpression,'match','tokens');
        % Store the length in cm
        TL(cnt) = str2double(tok{1})*1E-4; 
        % Average the resistance from the .mat files
        Ri = 0;
        % Sum the errors by 'quadrature'
        Ri_err = zeros(length(files),1);
        for j=1:nMatFiles
            % Load the matfile
            mfile = matfile(strcat(localPath,files(j).name)); 
            % Extract the voltage and current
            V = mfile.voltage;
            I = mfile.current;
            % fit to a line
            fitPoly1 = fit(I,V,'poly1');
            % The resistance is the slope of the line
            Ri = Ri + fitPoly1.p1;
            % Extract the uncertainty from the 95% confidence interval
            ci = confint(fitPoly1,0.95);
            Ri_err(i) = mean(abs(ci(:,1)-fitPoly1.p1));
        end
        % Store the average R in the array
        R(cnt) = Ri / nMatFiles;
        % Store the errors in the array
        Rerr(cnt) = norm(Ri_err);
        % Increase the counter
        cnt = cnt + 1;
    end
end

% Transpose the vectors to feed them to the optimizer
TL = TL';
R = R';
% Get the weights for fitting (use the inverse of the error)
wt = Rerr.^(-2);
% Normalize the weights to 1
wt = wt./max(wt);
% Fit the equation
fitPoly1 = fit(TL,R,'poly1','Weight',wt)
slope = fitPoly1.p1;
intercept = fitPoly1.p2;
% Get the 95% confidence interval
ci = confint(fitPoly1,0.95);
% Find the prediction intervals
pi = predint(fitPoly1,TL,0.95,'functional','on');

% Estimate the errors in the slope and intercept using the confidence
% interval
slope_err = max(abs(ci(:,1)-slope));
intercept_err = max(abs(ci(:,2)-intercept));

% Estimate the error in Rs
Rs = slope*contactWidth;
Rs_err = Rs*norm([slope/slope_err, contactWidthErr/contactWidth]);
% Estimate the transfer length and its uncertainty
LT = intercept/(2*slope);
LT_err = LT*abs([intercept_err/intercept, slope_err/slope]);
% Estimate the contact resistance from the intercept
Rc = intercept/2;
Rc_err = Rc*abs(intercept_err/intercept);
% Estimate the contact resistivity from the contact resistance and transfer
% length
rho_c = Rc*LT*contactWidth;
rho_c_err = rho_c*norm([Rc_err/Rc, LT_err/LT, contactWidthErr/contactWidth]);

% Plot the results
fig1 = figure;
set(gca, 'FontSize',14);
hold on
errorbar(TL,R,Rerr,'s','LineWidth',2.5,'MarkerSize',14);
plot(TL,fitPoly1(TL),'LineWidth',2.5);
% plot(TL,pi,':','LineWidth',2.5,'color','r'); % << Prediction intervals
hold off

% Add info text to the plot
info_str = "R_s\t= %.1ge ± %.2g Ohm /sq\nR_c\t= %.3g ± %.2g Ohm\nrho_c\t= " ...
    + "%.1g ± %.2g Ohm· cm^2";
info_str = sprintf(info_str, ....
    Rs, Rs_err, Rc, Rc_err,rho_c,rho_c_err);
info_str_sci = "R_s = %s \\Omega /sq\nR_c = %s \\Omega\n\\rho_c = " ...
    + "%s \\Omega \\cdot cm^2";
info_str_sci = sprintf(info_str_sci, ....
    formatScientificWithError(Rs, Rs_err), ...
    formatScientificWithError(Rc, Rc_err), ...
    formatScientificWithError(rho_c,rho_c_err));
info_txt = text(0.05,0.95,info_str_sci,...
    'Units','normalized',...
    'HorizontalAlignment','left',...
    'VerticalAlignment','top');
info_txt.Color  = 'k';
info_txt.FontSize = 14;
box on
xlabel('Length (cm)','FontSize',20);
ylabel('Resistance (Ohm)','FontSize',20);
title('TLM Resistance','FontSize',20);

% Save the results
print(strcat(outputPath,'\TLM.png'),'-dpng','-r300');
savefig(strcat(outputPath,'\TLM.fig'));
save(strcat(outputPath,'\TLM_results.mat'),'TL','R','Rerr','Rs','Rc',...
    'contactWidth','contactWidthErr','Rs_err','Rc_err','rho_c','rho_c_err');

fprintf(strcat(info_str,'\n'));


function [scientific] = formatScientificWithError(number,err,digits)
% Returns the number and its uncertainty in a scientific notation for the
% plot interpreter.
%
% Parameters
% ----------
% number : double
%   The main quantity
% err : double
%   Its associated uncertainty
% digits : int
%   The significative digits to print
%
% Returns
% -------
%   The formatted string
    if ~exist('digits','var')
        digits = 2;
    end
    ss       = sprintf('(%%.%df \\\\pm %%.%df) \\\\times10^{%%d}',(digits-1),digits);
    % get the exponent
    exponent = round(log10(abs(number)))-1;
    argument = number./power(10,exponent);
    arg_err  = err./power(10,exponent);
    scientific = sprintf(ss,argument,arg_err,exponent);
end
    
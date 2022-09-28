function [timeOff,out] = analyzeMats()
[filenames, pathname] = uigetfile2('*.*','Select files to run through', 'MultiSelect','on');

rateExos = [1:NaN:length(filenames)];
processExos = [1:NaN:length(filenames)];
timeOn = [1:NaN:length(filenames)];
starting = [1:NaN:length(filenames)];
ending = [1:NaN:length(filenames)];


for j=1:length(filenames)
    molecules(j) = load(filenames{j},'-mat');
    rateExos(j) = -1000*molecules(1,j).molec.digest_rate;
    processExos(j) = molecules(1,j).molec.end_length - molecules(1,j).molec.start_length;
    timeOn(j) = molecules(1,j).molec.x(end);

end
rates = rateExos';
process = processExos';

%Outputs-Average of Rate and Process
avgrate = mean(rates)
avgprocess = mean(process)

%Outputs-MAT files 
save('Rate.mat','rates')
save('Process.mat','process')

end

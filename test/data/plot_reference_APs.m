close all
clear all

file_listing = dir(['reference_traces' filesep '*.dat']);

for i=1:length(file_listing)
    
    % Skip the files that are to do with tissue sims
    if (~isempty(strfind(file_listing(i).name, 'tissue')))
        continue
    end
    
    file_listing(i).name
    
    d = importdata(['reference_traces' filesep file_listing(i).name]);
    if (isstruct(d))
        data = d.data;
    else
        data = d;
    end

    figure(i)
    plot(data(:,1), data(:,2),'b-')
    hold on
    title(file_listing(i).name)    
    xlabel('Time (ms)')
    xlabel('Voltage (mV)')   
end


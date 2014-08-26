file_listing = dir(['reference_traces' filesep '*.dat']);

for i=1:length(file_listing)
    file_listing(i).name
    d = importdata(['reference_traces' filesep file_listing(i).name]);
    data = d.data;
    
    figure
    plot(data(:,1), data(:,2),'b-')
    title(file_listing(i).name)
    xlabel('Time (ms)')
    xlabel('Voltage (mV)')    
end
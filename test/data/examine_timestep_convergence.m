% You can run this after "CalculateTimesteps" to look at the resulting
% action potential traces relative to the definitive reference traces.

close all
clear all

file_listing = dir(['reference_traces' filesep '*.dat']);

required_steps_data = importdata('required_steps.txt');

% Index in this vector (minus 1) is the solver code in the results files.
solver_mapping = {'CVODE (analytic Jacobian)',...
                  'CVODE (numerical Jacobian)',...
                  'Forward Euler',...
                  'Backward Euler',...
                  'Runge-Kutta (2nd order)',...
                  'Runge-Kutta (4th order)',...
                  'Rush-Larsen',...
                  'Generalised Rush-Larsen 1',...
                  'Generalised Rush-Larsen 2'};

for i=1:length(file_listing)
    
    % Skip the files that are to do with tissue sims
    if (~isempty(strfind(file_listing(i).name, 'tissue')))
        continue
    end
    
    %file_listing(i).name
    k = strfind(file_listing(i).name, '.dat');
    model_name = file_listing(i).name(1:k-1);
    
    d = importdata(['reference_traces' filesep file_listing(i).name]);
    if (isstruct(d))
        data = d.data;
    else
        data = d;
    end
    
    figure(i)    
    
    % Now see how many different timestep approximations we have for this
    hardcoded_results_folder = '/export/testoutput/Frontiers/CalculateTimesteps/';
    
    for gbu = 1:2
        % Good, bad, ugly errors
        subplot(2,1,gbu)
        plot(data(:,1), data(:,2),'-')
        hold all
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        legend_entries{1} = 'reference';
        
        model_rows = find(strcmp(model_name, required_steps_data.textdata));
        
        for i=0:8
            
            solver_rows = find(i==required_steps_data.data(:,1));
            relevant_rows = intersect(model_rows,solver_rows);
            if (isempty(relevant_rows))
                continue
            end
            
            folder_name = [hardcoded_results_folder model_name filesep num2str(i) filesep];
            refine_list = dir([folder_name model_name '_deltaT_*.dat']);
            for j=1:length(refine_list)
                
                % Work out the timesteps for the legend
                start_idx = strfind(refine_list(j).name, '_deltaT_') + 8;
                last_idx = strfind(refine_list(j).name, '.dat') - 1;
                dt = refine_list(j).name(start_idx:last_idx);
            
                % Work out the error metric associated with this
                dt_number = str2num(dt);
                
                timestep_rows = find(dt_number==required_steps_data.data(:,3));
                this_row = intersect(timestep_rows,relevant_rows);
                
                % Only look at the non-lookup table cases
                lookup_rows = find(0==required_steps_data.data(:,2));
                this_row = intersect(this_row,lookup_rows);                
                
                if (isempty(this_row))
                    continue
                else
                    assert(length(this_row)==1)
                    mrms_error = required_steps_data.data(this_row,11);
                    if (gbu==1)
                       if mrms_error > 0.05
                           continue
                       end
                       title([model_name ' MRMS <= 0.05'])
                    else 
                        if mrms_error <= 0.05
                            continue
                        end
                        title('MRMS > 0.05')
                    end
                end
                               
                d = importdata([folder_name refine_list(j).name]);
                
                if length(legend_entries)<7
                    linestyle = '-';
                else
                    linestyle = '--';
                end
                plot(d.data(:,1), d.data(:,2), linestyle)
                
                if i<=1
                    legend_entries{end+1} = [solver_mapping{i+1} ' Tol ' num2str(10^(-(2+dt_number))) ' MRMS = ' num2str(mrms_error)];
                else
                    legend_entries{end+1} = [solver_mapping{i+1} ' dt = ' dt ' MRMS = ' num2str(mrms_error)];
                end
            end
        end
        legend(legend_entries)
        clear legend_entries
    end
end



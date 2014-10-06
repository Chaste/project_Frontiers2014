% You can run this after "MonodomainCalculateTimesteps" to look at the resulting
% action potential traces relative to the definitive reference traces.

close all
clear all

file_listing = dir(['reference_traces' filesep '*.dat']);

hardcoded_results_folder = '/export/testoutput/Frontiers/MonodomainCalculateTimesteps/';

solver_listing = dir([hardcoded_results_folder 'beeler_reuter_model_1977']);
solver_listing = solver_listing(4:end); % Get rid of '.', '..' and '.chaste_deletable_folder'.

required_steps_data = importdata('required_steps_tissue.txt');

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
    
    % Skip the files that are to do with non-tissue,
    % or the wrong kind of tissue sims    
    if (isempty(strfind(file_listing(i).name, 'tissue')))
        continue
    elseif (isempty(strfind(file_listing(i).name, 'Reset')))
        continue
    elseif (isempty(strfind(file_listing(i).name, 'h_0.01')))
        continue
    elseif (~isempty(strfind(file_listing(i).name, 'pde_1')))
        continue
    elseif (~isempty(strfind(file_listing(i).name, 'pde_0.5')))
        continue
    elseif (~isempty(strfind(file_listing(i).name, 'pde_0.001')))
        continue
    end

    file_listing(i).name
    k = strfind(file_listing(i).name, '_tissue');
    model_name = file_listing(i).name(1:k-1);
    
    k_1 = strfind(file_listing(i).name, '_pde_');
    k_2 = strfind(file_listing(i).name, '_h_');
    pde_step = str2num(file_listing(i).name((k_1+5):k_2-1));  
    
    d = importdata(['reference_traces' filesep file_listing(i).name]);
    if (isstruct(d))
        data = d.data;
    else
        data = d;
    end
    
    figure(i)    
        
    for gbu = 1:3
        % Good, bad, ugly errors
        subplot(3,1,gbu)
        plot(data(:,1), data(:,2),'-')
        hold all
        xlabel('Time (ms)')
        ylabel('Voltage (mV)')
        legend_entries{1} = 'reference';
        
        model_rows = find(strcmp(model_name, required_steps_data.textdata));
        
        
        for i=0:8
            
            solver_rows = find(i==required_steps_data.data(:,2));
            relevant_rows = intersect(model_rows,solver_rows);
            if (isempty(relevant_rows))
                continue                
            end
                        
            folder_name = [hardcoded_results_folder model_name filesep solver_mapping{i+1} filesep];
            refine_list = dir([folder_name 'results_pde_' num2str(pde_step) '_ode_*']);
            
            for j=1:length(refine_list)
                
                % Work out the timesteps for the legend
                start_idx = strfind(refine_list(j).name, '_ode_') + 5;
                dt = refine_list(j).name(start_idx:end);
            
                % Work out the error metric associated with this
                dt_number = str2num(dt);
                
                timestep_rows = find(dt_number==required_steps_data.data(:,3));
                this_row = intersect(timestep_rows,relevant_rows);
                pde_rows = find(pde_step==required_steps_data.data(:,1));
                this_row = intersect(this_row,pde_rows);
                
                if (isempty(this_row))
                    continue
                else
                    assert(length(this_row)==1)
                    mrms_error = required_steps_data.data(this_row,11);
                    if (gbu==1)
                        title('Good')
                        if mrms_error > 0.01
                            continue
                        end
                    elseif (gbu==2)
                        title(['Bad ' model_name ' PDE dt = ' num2str(pde_step)])
                        if mrms_error < 0.01 || mrms_error > 0.05
                            continue
                        end
                    else % gbu==3
                        title('Ugly')
                        if mrms_error < 0.05
                            continue
                        end
                    end
                end
                
                try
                    data = importdata([folder_name refine_list(j).name filesep 'last_node_trace.dat']);
                catch
                    continue
                end
                
                if length(legend_entries)<7
                    linestyle = '-';
                else
                    linestyle = '--';
                end
                plot(data(:,1), data(:,2), linestyle) 
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
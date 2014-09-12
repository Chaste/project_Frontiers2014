close all
clear all

file_listing = dir(['reference_traces' filesep '*_GccOpt_tissue_*.dat']);

% Pull out the model names
for i=1:length(file_listing)
    file = file_listing(i).name;
    k = strfind(file, '_GccOpt_tissue_');
    model{i} = file(1:k-1);
end

model_list = unique(model);

% Let's look at summary files first for each model
for m=1:length(model_list)
    figure
    fprintf('Processing model: %s\n', model_list{m});
    for compiler=1:2
    
        if compiler==1
            prefix = '_tissue_';
        elseif compiler==2
            prefix = '_GccOpt_tissue_';
        end
        file_listing = dir(['reference_traces' filesep model_list{m} prefix '*.summary']);

        % Look at the two mesh resolutions separately.    
        for h = [0.01 0.001]
            pde_timesteps = [];
            % Find all the traces to do with this h
            for i=1:length(file_listing)
                file = file_listing(i).name;
                j = strfind(file, '_pde_');
                k = strfind(file, ['_h_' num2str(h)]);
                if (isempty(k))
                    continue
                end
                pde_timesteps = [pde_timesteps; str2num(file(j+5:k-1))];
            end

            % Right, now we have all the available combinations, let's plot the
            % convergence of the summary metrics
            for pde_idx = 1:length(pde_timesteps)
                pde_step = pde_timesteps(pde_idx);
                file_of_interest = ['reference_traces' filesep model_list{m} ...
                    prefix 'pde_' num2str(pde_step)  '_h_' num2str(h) '.summary'];
                d = importdata(file_of_interest);
                metrics(:,pde_idx) = d.data;
            end

            num_metrics = size(metrics,1);

            for i=1:num_metrics
                subplot(1,num_metrics,i)
                if compiler==1
                    linestyle = '.-';
                else
                    linestyle = '.--';
                end
                semilogx(pde_timesteps,metrics(i,:),linestyle)
                hold all
                if i == 3
                    title(strrep(model_list{m}, '_', ' '))
                end
                xlim([min(pde_timesteps) max(pde_timesteps)])
                xlabel('PDE Timestep (ms)')
                set(gca,'XTick',[0.001 0.01 0.1 1])
                ylabel(strrep(d.textdata{i}, '_', ' '))
            end        
            clear metrics
        end
    end
    legend('h = 0.01 cm Intel travis','h = 0.001 cm Intel travis','h = 0.01 cm GccOpt Skip','h = 0.001 cm GccOpt skip','Location','SouthWest')
end
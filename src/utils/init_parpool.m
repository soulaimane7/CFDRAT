% manage_parallel_pool.m
function init_parpool(enable_parallel)
% Starts or stops the parallel pool based on user selection.
if strcmp(enable_parallel, 'on')
    if isempty(gcp('nocreate'))
        fprintf('Starting parallel pool, please wait...\n');
        parpool;
    end
    if ~isempty(gcp('nocreate'))
        fprintf('Parallel pool started successfully.\n');
    else
        fprintf('Failed to start parallel pool, will run in serial mode.\n');
    end
else
    if ~isempty(gcp('nocreate'))
        fprintf('Parallel computing disabled. Shutting down existing pool...\n');
        delete(gcp('nocreate'));
    else
        fprintf('Parallel computing disabled.\n');
    end
end
end

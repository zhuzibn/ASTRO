% Reproducible regression baseline for the current ASTRO simulation path.
%
% Quick mode (default):
%   run('benchmark/run_baseline_benchmark.m')
%
% Current-size mode:
%   benchmark_mode = 'current';
%   run('benchmark/run_baseline_benchmark.m')

if ~exist('benchmark_mode', 'var')
    benchmark_mode = 'quick';
end
benchmark_mode = lower(char(benchmark_mode));
if ~ismember(benchmark_mode, {'current', 'quick'})
    error('benchmark_mode must be ''current'' or ''quick''.');
end

benchmark_version = '1.0';
benchmark_script = mfilename('fullpath');
benchmark_dir = fileparts(benchmark_script);
repo_root = fileparts(benchmark_dir);
input_dir = fullfile(benchmark_dir, 'baseline', 'input');
if exist('benchmark_output_root', 'var') && ~isempty(benchmark_output_root)
    benchmark_output_root = char(benchmark_output_root);
    output_dir = fullfile(benchmark_output_root, 'output');
    figure_dir = fullfile(benchmark_output_root, 'figures');
    metadata_file = fullfile(benchmark_output_root, 'benchmark_metadata.md');
    output_file = fullfile(output_dir, 'benchmark_result.mat');
else
    output_dir = fullfile(benchmark_dir, 'baseline', 'output');
    figure_dir = fullfile(benchmark_dir, 'baseline', 'figures');
    metadata_file = fullfile(benchmark_dir, 'baseline', 'baseline_metadata.md');
    output_file = fullfile(output_dir, 'baseline_result.mat');
end
addpath(repo_root);
astro_defaults = astro_default_config();

if ~exist(input_dir, 'dir'), mkdir(input_dir); end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
if ~exist(figure_dir, 'dir'), mkdir(figure_dir); end

git_commit = benchmark_command_output(sprintf('git -C "%s" rev-parse HEAD', repo_root));
git_branch = benchmark_command_output(sprintf('git -C "%s" branch --show-current', repo_root));
matlab_version = version;
timestamp = char(datetime('now', 'TimeZone', 'local', ...
    'Format', 'yyyy-MM-dd HH:mm:ss Z'));
operating_system = system_dependent('getos');
if strcmp(benchmark_mode, 'quick')
    exact_command = 'run(''benchmark/run_baseline_benchmark.m'')';
else
    exact_command = sprintf( ...
        'benchmark_mode=''%s''; run(''benchmark/run_baseline_benchmark.m'')', ...
        benchmark_mode);
end
if exist('benchmark_output_root', 'var') && ~isempty(benchmark_output_root)
    escaped_output_root = strrep(benchmark_output_root, '''', '''''');
    exact_command = sprintf('benchmark_output_root=''%s''; %s', ...
        escaped_output_root, exact_command);
end

gpu_info = struct('available', false, 'name', '', 'compute_capability', '', ...
    'driver_version', '', 'toolkit_version', NaN, ...
    'total_memory_bytes', NaN);
try
    selected_gpu = gpuDevice;
    gpu_info.available = true;
    gpu_info.name = selected_gpu.Name;
    gpu_info.compute_capability = selected_gpu.ComputeCapability;
    gpu_info.driver_version = selected_gpu.GraphicsDriverVersion;
    gpu_info.toolkit_version = selected_gpu.ToolkitVersion;
    gpu_info.total_memory_bytes = selected_gpu.TotalMemory;
catch gpu_error
    gpu_info.error = gpu_error.message;
end

% Save both distributions so neither benchmark mode depends on future RNG
% implementation details. baseline_atomtype.mat is the primary 20x30 input.
current_input_file = fullfile(input_dir, 'baseline_atomtype.mat');
quick_input_file = fullfile(input_dir, 'baseline_atomtype_quick.mat');
benchmark_ensure_atomtype(current_input_file, astro_defaults.natomW, ...
    astro_defaults.natomL, astro_defaults.compositionn, 104729, ...
    'Fixed 20x30 ASTRO baseline distribution; 1=RE/Gd, 0=TM/FeCo.');
benchmark_ensure_atomtype(quick_input_file, 3, 5, ...
    astro_defaults.compositionn, 104730, ...
    'Fixed 3x5 ASTRO quick baseline distribution; 1=RE/Gd, 0=TM/FeCo.');

if strcmp(benchmark_mode, 'current')
    natomW = astro_defaults.natomW;
    natomL = astro_defaults.natomL;
    atomtype_file = current_input_file;
else
    natomW = 3;
    natomL = 5;
    atomtype_file = quick_input_file;
end

benchmark_ran = false;
benchmark_elapsed_seconds = NaN;
benchmark_warning = '';
benchmark_write_metadata(metadata_file, false, benchmark_elapsed_seconds, ...
    benchmark_warning, benchmark_mode, timestamp, git_commit, git_branch, ...
    matlab_version, gpu_info, operating_system, exact_command, natomW, ...
    natomL, atomtype_file, output_file, figure_dir);

try
    if ~gpu_info.available
        error(['No compatible MATLAB GPU is available. The current ASTRO ' ...
            'solver path requires gpuArray support.']);
    end

    %% Shared production defaults from astro_default_config.m
    constantfile;
    clear gam
    astro_default_fields = fieldnames(astro_defaults);
    astro_conditional_fields = {'natom_mc_W', 'natom_mc_L', ...
        'fixededgeL', 'mxleft', 'myleft', 'mzleft', 'mxright', ...
        'myright', 'mzright'};
    for astro_default_index = 1:numel(astro_default_fields)
        astro_default_name = astro_default_fields{astro_default_index};
        if ~ismember(astro_default_name, astro_conditional_fields)
            eval([astro_default_name ' = astro_defaults.(astro_default_name);']);
        end
    end
    clear astro_default_fields astro_default_index astro_default_name ...
        astro_conditional_fields
    if strcmp(benchmark_mode, 'current')
        natomW = astro_defaults.natomW;
        natomL = astro_defaults.natomL;
    else
        natomW = 3;
        natomL = 5;
    end
    natom = natomW*natomL;
    if dipolemode == 3
        natom_mc_W = astro_defaults.natom_mc_W;
        natom_mc_L = astro_defaults.natom_mc_L;
    end
    if enablefixedge
        fixededgeL = astro_defaults.fixededgeL;
        mxleft = astro_defaults.mxleft;
        myleft = astro_defaults.myleft;
        mzleft = astro_defaults.mzleft;
        mxright = astro_defaults.mxright;
        myright = astro_defaults.myright;
        mzright = astro_defaults.mzright;
    end

    % Run the existing generation script, then replace its random atom
    % distribution and rebuild the same initial-state formula deterministically.
    rng(104729, 'twister');
    systemgeneration;
    fixed_input = load(atomtype_file);
    if fixed_input.natomW ~= natomW || fixed_input.natomL ~= natomL || ...
            fixed_input.compositionn ~= compositionn
        error('Fixed atom distribution metadata does not match benchmark mode.');
    end
    atomtype_ = fixed_input.atomtype_;
    seed_used = fixed_input.seed_used;
    [mx_init, my_init, mz_init] = benchmark_initial_state( ...
        atomtype_, natomW, natomL, dwcalc);

    gpusteps = round(gpusave/tstep);
    runtime = gpurun_number*gpusave;
    totstep = round(runtime/tstep);
    t = (0:savetstep:totstep)*tstep;
    if mod(gpusteps, savetstep) ~= 0
        error('gpusteps must be an integer multiple of savetstep.');
    end
    final_m_savestep = numel(t);
    if (SOT_DLT || SOT_FLT) && rk4 ~= 1
        error('Only RK4 is implemented for spin-torque-driven runs.');
    end

    enabled_terms = struct( ...
        'exchange', true, ...
        'anisotropy', true, ...
        'DMI', logical(DMIenable), ...
        'external_field_nonzero', any(Hext ~= 0), ...
        'thermal', logical(thermalenable), ...
        'dipole', logical(dipolemode), ...
        'SOT_DLT', logical(SOT_DLT), ...
        'SOT_FLT', logical(SOT_FLT), ...
        'STT_DLT', logical(STT_DLT), ...
        'STT_FLT', logical(STT_FLT), ...
        'fixed_edge', logical(enablefixedge));
    parameter_snapshot = struct( ...
        'Ksim', Ksim, 'Jgdgd', Jgdgd, 'Jfefe', Jfefe, ...
        'Jfegd', Jfegd, 'gTM', gTM, 'gRE', gRE, ...
        'musRE', musRE, 'musTM', musTM, 'Dsim', Dsim, ...
        'alp', alp, 'T', T, 'Hext', Hext, 'jcSOT', jcSOT, ...
        'jcSTT', jcSTT, 'thetaSH', thetaSH, 'etaSTT', etaSTT, ...
        'psjSHE', psjSHE, 'psjSTT', psjSTT, 'dwcalc', dwcalc);

    solver_timer = tic;
    integrate_llg;
    benchmark_elapsed_seconds = toc(solver_timer);
    benchmark_ran = true;

    final_mx = mmx(:, :, end);
    final_my = mmy(:, :, end);
    final_mz = mmz(:, :, end);
    final_spin_norm = sqrt(final_mx.^2 + final_my.^2 + final_mz.^2);
    spin_norm_min = min(final_spin_norm(:));
    spin_norm_max = max(final_spin_norm(:));
    spin_norm_mean = mean(final_spin_norm(:));
    spin_norm_max_abs_error = max(abs(final_spin_norm(:) - 1));
    avg_mx_time = squeeze(mean(mean(mmx, 1), 2));
    avg_my_time = squeeze(mean(mean(mmy, 1), 2));
    avg_mz_time = squeeze(mean(mean(mmz, 1), 2));
    avg_mx_final = avg_mx_time(end);
    avg_my_final = avg_my_time(end);
    avg_mz_final = avg_mz_time(end);
    array_size = size(mmx);
    gd_site_count = nnz(atomtype_ == 1);
    tm_site_count = nnz(atomtype_ == 0);
    [avg_gd_spin_time, avg_tm_spin_time] = benchmark_sublattice_averages( ...
        mmx, mmy, mmz, atomtype_);

    save(output_file, 'benchmark_version', 'benchmark_mode', 'git_commit', ...
        'git_branch', 'matlab_version', 'gpu_info', 'timestamp', ...
        'operating_system', 'natomW', 'natomL', 'd', 'tstep', ...
        'gpusave', 'gpurun_number', 'runtime', 'savetstep', 'bc', 'rk4', ...
        'compositionn', 'seed_used', 'atomtype_', 'mx_init', 'my_init', ...
        'mz_init', 'mmx', 'mmy', 'mmz', 't', 'final_mx', 'final_my', ...
        'final_mz', 'spin_norm_min', 'spin_norm_max', 'spin_norm_mean', ...
        'spin_norm_max_abs_error', 'avg_mx_time', 'avg_my_time', ...
        'avg_mz_time', 'avg_mx_final', 'avg_my_final', 'avg_mz_final', ...
        'enabled_terms', 'parameter_snapshot', 'benchmark_elapsed_seconds', ...
        'array_size', 'gd_site_count', 'tm_site_count', ...
        'avg_gd_spin_time', 'avg_tm_spin_time');

    benchmark_create_figures(figure_dir, atomtype_, mx_init, my_init, ...
        mz_init, final_mx, final_my, final_mz, t, avg_mx_time, ...
        avg_my_time, avg_mz_time, final_spin_norm);
    benchmark_write_metadata(metadata_file, benchmark_ran, ...
        benchmark_elapsed_seconds, benchmark_warning, benchmark_mode, ...
        timestamp, git_commit, git_branch, matlab_version, gpu_info, ...
        operating_system, exact_command, natomW, natomL, atomtype_file, ...
        output_file, figure_dir);
    fprintf('ASTRO %s baseline completed in %.6f seconds.\n', ...
        benchmark_mode, benchmark_elapsed_seconds);
    fprintf('Result: %s\n', output_file);
catch benchmark_error
    benchmark_warning = benchmark_error.message;
    benchmark_write_metadata(metadata_file, false, benchmark_elapsed_seconds, ...
        benchmark_warning, benchmark_mode, timestamp, git_commit, git_branch, ...
        matlab_version, gpu_info, operating_system, exact_command, natomW, ...
        natomL, atomtype_file, output_file, figure_dir);
    rethrow(benchmark_error);
end

function text = benchmark_command_output(command)
[status, text] = system(command);
text = strtrim(text);
if status ~= 0
    text = 'unavailable';
end
end

function benchmark_ensure_atomtype(file_path, natomW, natomL, ...
        compositionn, seed_used, description)
if exist(file_path, 'file')
    fixed = load(file_path);
    required = {'atomtype_', 'natomW', 'natomL', 'compositionn', ...
        'seed_used', 'description'};
    if ~all(isfield(fixed, required)) || fixed.natomW ~= natomW || ...
            fixed.natomL ~= natomL || fixed.compositionn ~= compositionn
        error('Existing fixed atom distribution is incompatible: %s', file_path);
    end
    return
end
old_rng = rng;
rng(seed_used, 'twister');
natom = natomW*natomL;
gd_indices = randperm(natom, round(natom*compositionn));
atomtype_ = zeros(natomW, natomL);
atomtype_(gd_indices) = 1;
rng(old_rng);
save(file_path, 'atomtype_', 'natomW', 'natomL', 'compositionn', ...
    'seed_used', 'description');
end

function [mx_init, my_init, mz_init] = benchmark_initial_state( ...
        atomtype_, natomW, natomL, dwcalc)
mx_init = zeros(natomW, natomL);
my_init = zeros(natomW, natomL);
mz_init = zeros(natomW, natomL);
phi_ = 0;
for ctL = 1:natomL
    for ctW = 1:natomW
        is_tm = atomtype_(ctW, ctL) == 0;
        if dwcalc && ctL >= round(natomL/2)
            theta_degrees = 5 + 180*is_tm;
        else
            theta_degrees = 5 + 180*(~is_tm);
        end
        theta_ = theta_degrees/180*pi;
        mx_init(ctW, ctL) = sin(theta_)*cos(phi_);
        my_init(ctW, ctL) = sin(theta_)*sin(phi_);
        mz_init(ctW, ctL) = cos(theta_);
    end
end
end

function [gd_average, tm_average] = benchmark_sublattice_averages( ...
        mmx, mmy, mmz, atomtype_)
saved_count = size(mmx, 3);
gd_average = NaN(saved_count, 3);
tm_average = NaN(saved_count, 3);
gd_mask = atomtype_ == 1;
tm_mask = atomtype_ == 0;
for saved_index = 1:saved_count
    mx_slice = mmx(:, :, saved_index);
    my_slice = mmy(:, :, saved_index);
    mz_slice = mmz(:, :, saved_index);
    if any(gd_mask(:))
        gd_average(saved_index, :) = [mean(mx_slice(gd_mask)), ...
            mean(my_slice(gd_mask)), mean(mz_slice(gd_mask))];
    end
    if any(tm_mask(:))
        tm_average(saved_index, :) = [mean(mx_slice(tm_mask)), ...
            mean(my_slice(tm_mask)), mean(mz_slice(tm_mask))];
    end
end
end

function benchmark_create_figures(figure_dir, atomtype_, mx_init, ...
        my_init, mz_init, final_mx, final_my, final_mz, t, avg_mx_time, ...
        avg_my_time, avg_mz_time, final_spin_norm)
benchmark_matrix_figure(atomtype_, 'Atom type (1=RE/Gd, 0=TM/FeCo)', ...
    fullfile(figure_dir, 'atomtype_matrix.png'), [0, 1]);
benchmark_matrix_figure(mz_init, 'Initial m_z', ...
    fullfile(figure_dir, 'initial_mz.png'), [-1, 1]);
benchmark_matrix_figure(final_mz, 'Final m_z', ...
    fullfile(figure_dir, 'final_mz.png'), [-1, 1]);

fig = figure('Visible', 'off', 'Color', 'white');
plot(t*1e12, avg_mx_time, 'LineWidth', 1.5);
hold on;
plot(t*1e12, avg_my_time, 'LineWidth', 1.5);
plot(t*1e12, avg_mz_time, 'LineWidth', 1.5);
xlabel('Time (ps)');
ylabel('Average spin component');
legend('m_x', 'm_y', 'm_z', 'Location', 'best');
grid on;
title('Average spin versus saved time');
exportgraphics(fig, fullfile(figure_dir, 'average_spin_vs_time.png'), ...
    'Resolution', 150);
close(fig);

fig = figure('Visible', 'off', 'Color', 'white');
histogram(final_spin_norm(:), 40);
xlabel('Final spin norm');
ylabel('Site count');
grid on;
title('Final spin norm distribution');
exportgraphics(fig, fullfile(figure_dir, 'spin_norm_histogram.png'), ...
    'Resolution', 150);
close(fig);

% Keep these variables referenced to make the figure API explicit.
assert(isequal(size(mx_init), size(my_init)));
assert(isequal(size(final_mx), size(final_my)));
end

function benchmark_matrix_figure(matrix, plot_title, file_path, limits)
fig = figure('Visible', 'off', 'Color', 'white');
% flipud maps physical length to x and reverses width relative to MATLAB rows.
imagesc(1:size(matrix, 2), 1:size(matrix, 1), flipud(matrix));
axis image;
set(gca, 'YDir', 'normal');
xlabel('Length index L');
ylabel('Width index W (reversed from matrix row order)');
title(plot_title);
colorbar;
clim(limits);
exportgraphics(fig, file_path, 'Resolution', 150);
close(fig);
end

function benchmark_write_metadata(file_path, benchmark_ran, elapsed, ...
        warning_text, mode, timestamp, commit, branch, matlab_version, ...
        gpu_info, operating_system, exact_command, natomW, natomL, ...
        atomtype_file, output_file, figure_dir)
fid = fopen(file_path, 'w');
if fid < 0
    error('Could not write metadata file: %s', file_path);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# ASTRO Baseline Metadata\n\n');
fprintf(fid, '- Benchmark date/time: `%s`\n', timestamp);
fprintf(fid, '- Git commit: `%s`\n', commit);
fprintf(fid, '- Git branch: `%s`\n', branch);
fprintf(fid, '- MATLAB version: `%s`\n', matlab_version);
fprintf(fid, '- GPU available: `%d`\n', gpu_info.available);
fprintf(fid, '- GPU: `%s`\n', gpu_info.name);
fprintf(fid, '- GPU compute capability: `%s`\n', gpu_info.compute_capability);
fprintf(fid, '- GPU driver: `%s`\n', gpu_info.driver_version);
fprintf(fid, '- CUDA toolkit version reported by MATLAB: `%.4f`\n', ...
    gpu_info.toolkit_version);
fprintf(fid, '- Operating system: `%s`\n', operating_system);
fprintf(fid, '- Benchmark mode: `%s`\n', mode);
fprintf(fid, '- Exact MATLAB expression: `%s`\n', exact_command);
fprintf(fid, '- Lattice size: `%d x %d` (W x L)\n', natomW, natomL);
fprintf(fid, '- Lattice constant: `0.4e-9 m`\n');
fprintf(fid, ['- Time settings: `tstep=2e-16 s`, `gpusave=1e-12 s`, ' ...
    '`gpurun_number=2`, `runtime=2e-12 s`, `savetstep=100`\n']);
fprintf(fid, '- Composition: `compositionn=0.1` (10%% RE/Gd)\n');
fprintf(fid, '- Fixed atom distribution: `%s`\n', atomtype_file);
fprintf(fid, ['- Enabled: exchange, anisotropy, DMI, SOT damping-like ' ...
    'torque\n']);
fprintf(fid, ['- Disabled/zero: thermal, dipole, SOT field-like, both STT ' ...
    'terms, fixed edges, external field\n']);
fprintf(fid, '- Result file: `%s`\n', output_file);
fprintf(fid, '- Figure directory: `%s`\n', figure_dir);
fprintf(fid, '- Benchmark actually ran: `%d`\n', benchmark_ran);
if benchmark_ran
    fprintf(fid, '- Solver runtime: `%.6f seconds`\n', elapsed);
else
    fprintf(fid, '- Solver runtime: unavailable\n');
end
if isempty(warning_text)
    fprintf(fid, '- Warnings: none recorded\n');
else
    fprintf(fid, '- Warnings/errors: `%s`\n', strrep(warning_text, '`', ''''));
end
fprintf(fid, '- Missing quantities: energy is not computed by current ASTRO output.\n');
fprintf(fid, '\n## Deviations From `main.m`\n\n');
fprintf(fid, ['- The runner does not execute `main.m` because it starts with ' ...
    '`clear all` and `rng(''shuffle'')`.\n']);
fprintf(fid, ['- It loads `astro_default_config.m`, runs the existing ' ...
    '`systemgeneration.m`, then replaces the random `atomtype_` with the ' ...
    'fixed input and rebuilds the same initial-state formula.\n']);
fprintf(fid, ['- It calls the unchanged `integrate_llg.m`, `field_calc.m`, and ' ...
    '`dipole_calc.m` path.\n']);
fprintf(fid, ['- It saves an explicit benchmark MAT-file instead of the root-level ' ...
    '`final.mat`; no `final.mat` is required or overwritten.\n']);
fprintf(fid, ['- Quick mode changes only lattice size from `20 x 30` to `3 x 5`; ' ...
    'all listed physics and time settings remain unchanged.\n']);
fprintf(fid, '\n## Interpretation\n\n');
fprintf(fid, ['This is a regression baseline for current code behavior. It is not ' ...
    'a proof or full validation of the physics.\n']);
end

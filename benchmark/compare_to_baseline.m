function passed = compare_to_baseline(future_result_file)
%COMPARE_TO_BASELINE Compare a future ASTRO result with the stored baseline.
%
%   addpath('benchmark')
%   compare_to_baseline('path/to/future_result.mat')

benchmark_dir = fileparts(mfilename('fullpath'));
baseline_file = fullfile(benchmark_dir, 'baseline', 'output', ...
    'baseline_result.mat');
if nargin < 1 || isempty(future_result_file)
    error('Provide the path to a future benchmark result MAT-file.');
end
if ~exist(baseline_file, 'file')
    error('Baseline result does not exist: %s', baseline_file);
end
if ~exist(future_result_file, 'file')
    error('Future result does not exist: %s', future_result_file);
end

tol_avg_spin = 1e-8;
tol_spin_array = 1e-6;
tol_spin_norm = 1e-8;
baseline = load(baseline_file);
future = load(future_result_file);

required = {'atomtype_', 'mmx', 'mmy', 'mmz', 'avg_mx_final', ...
    'avg_my_final', 'avg_mz_final', 'spin_norm_min', 'spin_norm_max', ...
    'spin_norm_mean', 'spin_norm_max_abs_error'};
missing_baseline = required(~isfield(baseline, required));
missing_future = required(~isfield(future, required));
if ~isempty(missing_baseline) || ~isempty(missing_future)
    error('Missing required fields. Baseline: %s; future: %s', ...
        strjoin(missing_baseline, ', '), strjoin(missing_future, ', '));
end

sizes_match = isequal(size(baseline.mmx), size(future.mmx)) && ...
    isequal(size(baseline.mmy), size(future.mmy)) && ...
    isequal(size(baseline.mmz), size(future.mmz));
atomtype_matches = isequal(baseline.atomtype_, future.atomtype_);
baseline_final_average = [baseline.avg_mx_final, baseline.avg_my_final, ...
    baseline.avg_mz_final];
future_final_average = [future.avg_mx_final, future.avg_my_final, ...
    future.avg_mz_final];
average_difference = future_final_average - baseline_final_average;
max_average_difference = max(abs(average_difference));

final_spin_array_l2 = NaN;
final_spin_array_max_abs = Inf;
trajectory_max_abs = Inf;
if sizes_match
    dx = future.mmx - baseline.mmx;
    dy = future.mmy - baseline.mmy;
    dz = future.mmz - baseline.mmz;
    final_dx = dx(:, :, end);
    final_dy = dy(:, :, end);
    final_dz = dz(:, :, end);
    final_spin_array_l2 = sqrt(norm(final_dx(:))^2 + ...
        norm(final_dy(:))^2 + norm(final_dz(:))^2);
    final_spin_array_max_abs = max([abs(final_dx(:)); ...
        abs(final_dy(:)); abs(final_dz(:))]);
    trajectory_max_abs = max([abs(dx(:)); abs(dy(:)); abs(dz(:))]);
end

baseline_norm_stats = [baseline.spin_norm_min, baseline.spin_norm_max, ...
    baseline.spin_norm_mean, baseline.spin_norm_max_abs_error];
future_norm_stats = [future.spin_norm_min, future.spin_norm_max, ...
    future.spin_norm_mean, future.spin_norm_max_abs_error];
norm_stats_difference = future_norm_stats - baseline_norm_stats;
max_norm_stats_difference = max(abs(norm_stats_difference));

passed = sizes_match && atomtype_matches && ...
    max_average_difference <= tol_avg_spin && ...
    final_spin_array_max_abs <= tol_spin_array && ...
    max_norm_stats_difference <= tol_spin_norm;

fprintf('Array sizes match: %d\n', sizes_match);
fprintf('atomtype_ matches exactly: %d\n', atomtype_matches);
fprintf('Final average difference [mx my mz]: [% .6e % .6e % .6e]\n', ...
    average_difference);
fprintf('Final spin-array L2 difference: %.6e\n', final_spin_array_l2);
fprintf('Final spin-array maximum absolute difference: %.6e\n', ...
    final_spin_array_max_abs);
fprintf('Full saved trajectory maximum absolute difference: %.6e\n', ...
    trajectory_max_abs);
fprintf(['Spin norm statistic difference [min max mean max_abs_error]: ' ...
    '[% .6e % .6e % .6e % .6e]\n'], norm_stats_difference);
fprintf('Tolerances: avg=%.1e, array=%.1e, norm=%.1e\n', ...
    tol_avg_spin, tol_spin_array, tol_spin_norm);
fprintf('PASS: %d\n', passed);
end

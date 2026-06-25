% Deterministic ASTRO invariant smoke checks.
%
% Batch command:
%   matlab -batch "run('tests/run_smoke_tests.m')"

smoke_script = mfilename('fullpath');
tests_dir = fileparts(smoke_script);
repo_root = fileparts(tests_dir);
addpath(repo_root);

benchmark_mode = 'quick';
benchmark_output_root = tempname;
run(fullfile(repo_root, 'benchmark', 'run_baseline_benchmark.m'));

result_file = fullfile(benchmark_output_root, 'output', 'benchmark_result.mat');
fixed_atomtype_file = fullfile(repo_root, 'benchmark', 'baseline', ...
    'input', 'baseline_atomtype_quick.mat');

result = load(result_file);
fixed_atomtype = load(fixed_atomtype_file);

assert_smoke(isequal(result.atomtype_, fixed_atomtype.atomtype_), ...
    'Smoke atom distribution does not match fixed quick benchmark input.');
assert_smoke(result.natomW == fixed_atomtype.natomW && ...
    result.natomL == fixed_atomtype.natomL, ...
    'Smoke atom distribution dimensions do not match fixed input metadata.');
assert_smoke(result.compositionn == fixed_atomtype.compositionn, ...
    'Smoke composition does not match fixed input metadata.');
assert_smoke(all(result.atomtype_(:) == 0 | result.atomtype_(:) == 1), ...
    'Smoke atom distribution contains values other than 0 and 1.');

expected_gd_site_count = round(result.natomW*result.natomL* ...
    result.compositionn);
expected_tm_site_count = result.natomW*result.natomL - ...
    expected_gd_site_count;
assert_smoke(nnz(result.atomtype_ == 1) == expected_gd_site_count, ...
    'Smoke Gd site count does not match deterministic composition.');
assert_smoke(nnz(result.atomtype_ == 0) == expected_tm_site_count, ...
    'Smoke TM site count does not match deterministic composition.');
assert_smoke(result.gd_site_count == expected_gd_site_count && ...
    result.tm_site_count == expected_tm_site_count, ...
    'Smoke result site-count metadata does not match atom distribution.');

expected_saved_count = round(result.runtime/result.tstep/result.savetstep) + 1;
expected_trajectory_size = [result.natomW, result.natomL, ...
    expected_saved_count];
assert_smoke(isequal(size(result.mmx), expected_trajectory_size), ...
    'Smoke mmx trajectory dimensions are incorrect.');
assert_smoke(isequal(size(result.mmy), expected_trajectory_size), ...
    'Smoke mmy trajectory dimensions are incorrect.');
assert_smoke(isequal(size(result.mmz), expected_trajectory_size), ...
    'Smoke mmz trajectory dimensions are incorrect.');
assert_smoke(numel(result.t) == expected_saved_count, ...
    'Smoke saved time count does not match expected saved count.');
assert_smoke(numel(result.t) == size(result.mmx, 3) && ...
    numel(result.t) == size(result.mmy, 3) && ...
    numel(result.t) == size(result.mmz, 3), ...
    'Smoke saved time count does not match saved state count.');

assert_smoke(all(isfinite(result.mmx(:))) && ...
    all(isfinite(result.mmy(:))) && all(isfinite(result.mmz(:))), ...
    'Smoke trajectory contains non-finite spin values.');

spin_norm = sqrt(result.mmx.^2 + result.mmy.^2 + result.mmz.^2);
spin_norm_tolerance = 1e-10;
assert_smoke(max(abs(spin_norm(:) - 1)) <= spin_norm_tolerance, ...
    'Smoke spin normalization check failed.');

time_tolerance = eps(result.runtime);
assert_smoke(result.t(1) == 0, 'Smoke saved times do not start at t=0.');
assert_smoke(abs(result.t(end) - result.runtime) <= time_tolerance, ...
    'Smoke saved times do not end at runtime.');
assert_smoke(all(diff(result.t) > 0), ...
    'Smoke saved times are not strictly increasing.');

fprintf('ASTRO smoke checks passed.\n');
fprintf('Result checked: %s\n', result_file);
fprintf('Saved states: %d\n', expected_saved_count);
fprintf('Max spin-norm error: %.6e\n', max(abs(spin_norm(:) - 1)));

function assert_smoke(condition, message)
if ~condition
    error(message);
end
end

% Deterministic ASTRO boundary and field-term checks.
%
% Batch command:
%   matlab -batch "run('tests/run_field_tests.m')"

field_test_script = mfilename('fullpath');
tests_dir = fileparts(field_test_script);
repo_root = fileparts(tests_dir);
addpath(repo_root);

assert_exchange_neighbor_coefficients();
assert_nonperiodic_boundary_field_terms();
assert_periodic_boundary_field_terms();
assert_anisotropy_term();
assert_dmi_term();
assert_external_field_term();

fprintf('ASTRO deterministic field checks passed.\n');

function assert_exchange_neighbor_coefficients()
atomtype_ = [1 0 1; 0 1 0];
bc = 1;
Jgdgd = 11;
Jfefe = 22;
Jfegd = 33;

[AsimnextL, AsimpreviousL, AsimnextW, AsimpreviousW] = ...
    astro_exchange_neighbor_coefficients(atomtype_, bc, Jgdgd, Jfefe, Jfegd);

% W is the first dimension and L is the second.  With bc=1, L-edge links
% and W-edge links are zero; unlike neighbors use Jfegd.
assert_close(AsimnextL, [33 33 0; 33 33 0], ...
    'Non-periodic AsimnextL coefficients changed.');
assert_close(AsimpreviousL, [0 33 33; 0 33 33], ...
    'Non-periodic AsimpreviousL coefficients changed.');
assert_close(AsimnextW, [0 0 0; 33 33 33], ...
    'Non-periodic AsimnextW coefficients changed.');
assert_close(AsimpreviousW, [33 33 33; 0 0 0], ...
    'Non-periodic AsimpreviousW coefficients changed.');

bc = 0;
[AsimnextL, AsimpreviousL, AsimnextW, AsimpreviousW] = ...
    astro_exchange_neighbor_coefficients(atomtype_, bc, Jgdgd, Jfefe, Jfegd);

% Periodic boundaries wrap.  Same RE/Gd links use Jgdgd, same TM/FeCo links
% use Jfefe, and mixed links use Jfegd.
assert_close(AsimnextL, [33 33 11; 33 33 22], ...
    'Periodic AsimnextL coefficients changed.');
assert_close(AsimpreviousL, [11 33 33; 22 33 33], ...
    'Periodic AsimpreviousL coefficients changed.');
assert_close(AsimnextW, [33 33 33; 33 33 33], ...
    'Periodic AsimnextW coefficients changed.');
assert_close(AsimpreviousW, [33 33 33; 33 33 33], ...
    'Periodic AsimpreviousW coefficients changed.');
end

function assert_nonperiodic_boundary_field_terms()
mmx = [1 2; 3 4];
mmy = zeros(2, 2);
mmz = zeros(2, 2);
muigpu = [10 20; 40 80];
AsimnextL = 2*ones(2, 2);
AsimpreviousL = 3*ones(2, 2);
AsimnextW = 5*ones(2, 2);
AsimpreviousW = 7*ones(2, 2);

fields = astro_evaluate_field_terms(mmx, mmy, mmz, AsimnextL, ...
    AsimpreviousL, AsimnextW, AsimpreviousW, muigpu, 0, 0, [0 0 0], 1);

% For bc=1 at (W=1,L=1), previous-L and next-W neighbors are zeroed:
% hex_x = -(2*2 + 3*0 + 5*0 + 7*3)/10 = -2.5 T.
% The other entries follow the same edge-zeroing formula.
expected_exchange_x = [-25/10, -31/20; -13/40, -19/80];
assert_close(fields.exchange.x, expected_exchange_x, ...
    'Non-periodic exchange field changed.');
assert_close(fields.total.x, expected_exchange_x, ...
    'Non-periodic deterministic total field changed.');
end

function assert_periodic_boundary_field_terms()
mmx = [1 2; 3 4];
mmy = zeros(2, 2);
mmz = zeros(2, 2);
muigpu = [10 20; 40 80];
AsimnextL = 2*ones(2, 2);
AsimpreviousL = 3*ones(2, 2);
AsimnextW = 5*ones(2, 2);
AsimpreviousW = 7*ones(2, 2);

fields = astro_evaluate_field_terms(mmx, mmy, mmz, AsimnextL, ...
    AsimpreviousL, AsimnextW, AsimpreviousW, muigpu, 0, 0, [0 0 0], 0);

% For bc=0 at (W=1,L=1), all four neighbors wrap:
% hex_x = -(2*2 + 3*2 + 5*3 + 7*3)/10 = -4.6 T.
expected_exchange_x = [-46/10, -53/20; -32/40, -39/80];
assert_close(fields.exchange.x, expected_exchange_x, ...
    'Periodic exchange field changed.');
end

function assert_anisotropy_term()
mmx = zeros(2, 2);
mmy = zeros(2, 2);
mmz = [1 -2; 3 -4];
muigpu = [2 4; 5 10];

fields = astro_evaluate_field_terms(mmx, mmy, mmz, zeros(2, 2), ...
    zeros(2, 2), zeros(2, 2), zeros(2, 2), muigpu, 5, 0, [0 0 0], 1);

% hani_z = 2*Ksim/mui * m_z.
assert_close(fields.anisotropy.x, zeros(2, 2), ...
    'Anisotropy x field should remain zero.');
assert_close(fields.anisotropy.y, zeros(2, 2), ...
    'Anisotropy y field should remain zero.');
assert_close(fields.anisotropy.z, [5 -5; 6 -4], ...
    'Anisotropy z field changed.');
end

function assert_dmi_term()
mmx = [1 2 3; 4 5 6];
mmy = [7 8 9; 10 11 12];
mmz = [13 14 15; 16 17 18];
muigpu = [2 4 5; 8 10 20];
Dsim = 6;

fields = astro_evaluate_field_terms(mmx, mmy, mmz, zeros(2, 3), ...
    zeros(2, 3), zeros(2, 3), zeros(2, 3), muigpu, 0, Dsim, [0 0 0], 1);

% With bc=1: hdmi_x = D/mui*(-m_z,nextL + m_z,previousL),
% hdmi_y = D/mui*(-m_z,nextW + m_z,previousW), and
% hdmi_z = D/mui*(m_x,nextL - m_x,previousL + m_y,nextW - m_y,previousW).
assert_close(fields.dmi.x, [-42 -3 84/5; -102/8 -12/10 102/20], ...
    'DMI x field changed.');
assert_close(fields.dmi.y, [48 51/2 108/5; -78/8 -84/10 -90/20], ...
    'DMI y field changed.');
assert_close(fields.dmi.z, [-24 -54/4 -84/5; 72/8 60/10 24/20], ...
    'DMI z field changed.');
end

function assert_external_field_term()
mmx = zeros(2, 2);
mmy = zeros(2, 2);
mmz = zeros(2, 2);
muigpu = ones(2, 2);
Hext = [0.25 -0.5 0.75];

fields = astro_evaluate_field_terms(mmx, mmy, mmz, zeros(2, 2), ...
    zeros(2, 2), zeros(2, 2), zeros(2, 2), muigpu, 0, 0, Hext, 1);

assert_close(fields.external.x, 0.25*ones(2, 2), ...
    'External x field changed.');
assert_close(fields.external.y, -0.5*ones(2, 2), ...
    'External y field changed.');
assert_close(fields.external.z, 0.75*ones(2, 2), ...
    'External z field changed.');
assert_close(fields.total.x, fields.external.x, ...
    'External-only total x field changed.');
assert_close(fields.total.y, fields.external.y, ...
    'External-only total y field changed.');
assert_close(fields.total.z, fields.external.z, ...
    'External-only total z field changed.');
end

function assert_close(actual, expected, message)
tolerance = 1e-12;
if ~isequal(size(actual), size(expected)) || ...
        any(abs(gather(actual) - expected) > tolerance, 'all')
    error(message);
end
end

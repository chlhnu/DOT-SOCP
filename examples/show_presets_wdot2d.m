function show_presets_wdot2d(varargin)
%SHOW_PRESETS_WDOT2D Display and demonstrate available problem and algorithm presets
%
%   show_presets_wdot2d()              - Show all available presets
%   show_presets_wdot2d('problems')    - Show problem types
%   show_presets_wdot2d('algorithms')  - Show algorithms only
%   show_presets_wdot2d('usage')       - Show usage examples
%
% SEE ALSO: demo_wdot2d, get_example, solver_wdotsocp2d

if nargin == 0 || strcmp(varargin{1}, 'problems')
    show_problems();
end

if nargin == 0 || strcmp(varargin{1}, 'algorithms')
    show_algorithms();
end

if nargin == 0 || strcmp(varargin{1}, 'usage')
    show_usage();
end

end

function show_problems()
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Problem Presets (set in demo_wdot2d.m: problem = "...")\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('Type 1 - Barrier problems:\n');
fprintf('  "example8"      - Love heart obstacle\n');
fprintf('  "example9"      - Maze obstacle\n');
fprintf('  "circle-pillar" - Circle and pillars obstacle\n');
fprintf('  "maze14"        - Maze (Papadakis et al. 2014)\n\n');

fprintf('  Usage:\n');
fprintf('      [rho0, rho1] = get_example(problem, nx, ny);\n');
fprintf('      [weight, barrier] = get_weight_from_barrier(problem, nx, ny, nt);\n');
fprintf('      [rho0, rho1, barrierh] = ensure_barrier_validity(rho0, rho1, barrier);\n\n');

fprintf('────────────────────────────────────────────────────────────────\n\n');

fprintf('Type 2 - General weighted problems:\n');
fprintf('  Densities: "example1"-"example4", "example8", "example9", "circle"\n');
fprintf('  Weights:   "circle", "circleInv"\n\n');

fprintf('  Usage:\n');
fprintf('      [rho0, rho1] = get_example(problem, nx, ny);\n');
fprintf('      weight = get_weight(weight_name, nx, ny, nt);\n\n');
end

function show_algorithms()
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Algorithm Presets (set in demo_wdot2d.m: algorithm = "...")\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('  "inpALM"     - Inexact Proximal ALM (recommended)\n');
fprintf('  "ALG2"       - Augmented Lagrangian Algorithm (ALG2) splitting scheme\n');
fprintf('  "acc-pADMM"  - Accelerated pADMM\n\n');
end

function show_usage()
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Usage Guide\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('1. Edit demo_wdot2d.m, choose problem type:\n\n');

fprintf('   Type 1 (barrier problem):\n');
fprintf('       problem = "circle-pillar";\n');
fprintf('       [rho0, rho1] = get_example(problem, nx, ny);\n');
fprintf('       [weight, barrier] = get_weight_from_barrier(problem, nx, ny, nt);\n');
fprintf('       [rho0, rho1, barrierh] = ensure_barrier_validity(rho0, rho1, barrier);\n\n');

fprintf('   Type 2 (general weighted problem):\n');
fprintf('       problem = "example1"; weight_name = "circle";\n');
fprintf('       [rho0, rho1] = get_example(problem, nx, ny);\n');
fprintf('       weight = get_weight(weight_name, nx, ny, nt);\n\n');

fprintf('2. Run:   >> demo_wdot2d\n\n');

fprintf('3. Common parameters:\n');
fprintf('       nt = 33;           %% Time resolution (2^k+1)\n');
fprintf('       nx = 129;          %% Spatial resolution (2^k+1)\n');
fprintf('       levelN = 3;        %% Multilevel depth\n');
fprintf('   Note: Grid sizes 2^k+1 (k >= levelN) facilitate multilevel coarsening,\n');
fprintf('         which can be overcome with more flexible grid handling.\n\n');
end

function show_presets_dot1d(varargin)
%SHOW_PRESETS_DOT1D Display and demonstrate available problem and algorithm presets
%
%   show_presets_dot1d()              - Show all available presets
%   show_presets_dot1d('problems')    - Show problem types only
%   show_presets_dot1d('algorithms')  - Show algorithms only
%   show_presets_dot1d('usage')       - Show usage examples
%
% SEE ALSO: demo_dot1d, get_example, solver_dotsocp1d

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
fprintf('Problem Presets (set in demo_dot1d.m: problem = "...")\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('  "gaussian"   - Gaussian distribution transport\n');
fprintf('  "box"        - Box-shaped distribution transport\n\n');

fprintf('  Usage:\n');
fprintf('      [rho0, rho1] = get_example(problem, nx, lowerBound);\n');
fprintf('          lowerBound - Optional. Add lower bound to densities (default: 0)\n\n');
end

function show_algorithms()
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Algorithm Presets (set in demo_dot1d.m: algo = "...")\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('  "inpALM"     - inexact proximal ALM\n\n');
fprintf('  "ALG2"       - Augmented Lagrangian Algorithm (ALG2) splitting scheme\n');
end

function show_usage()
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Usage Guide\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

fprintf('1. Edit demo_dot1d.m, change these lines:\n');
fprintf('       problem = "gaussian";   %% Try: "gaussian" or "box"\n');
fprintf('       algo = "inpALM";\n\n');

fprintf('2. Run:   >> demo_dot1d\n\n');

fprintf('3. Common parameters to adjust:\n');
fprintf('       nt = 33;       %% Time resolution (2^k+1 for multilevel coarsening)\n');
fprintf('       nx = 1025;     %% Spatial resolution (2^k+1 for multilevel coarsening)\n');
fprintf('       levelN = 3;    %% Multilevel depth\n');
fprintf('   Note: Grid sizes of form 2^k+1 (k >= levelN) are only convenient for multilevel coarsening,\n');
fprintf('         which can be overcome with more flexible grid handling.\n\n');
end

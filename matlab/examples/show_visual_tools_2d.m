function show_visual_tools_2d()
%SHOW_VISUAL_TOOLS_2D Display available visualization tools for 2D problems
%
%   show_visual_tools_2d() displays all available visualization
%   functions and their usage for analyzing DOT 2D and WDOT 2D results.
%
% SEE ALSO: show_evolution_2d, show_movement_2d, show_residualCurve

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('  Visualization Tools for 2D Problems (DOT/WDOT)\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

show_evolution();
fprintf('────────────────────────────────────────────────────────────────\n\n');
show_movement();
fprintf('────────────────────────────────────────────────────────────────\n\n');
show_convergence();

end

function show_evolution()
fprintf('1. Density Evolution: show_evolution_2d(rho, showFunc, title, barrier)\n\n');

fprintf('   Visualization types (showFunc):\n');
fprintf('      "imshow"    - Image view (default)\n');
fprintf('      "mesh"      - 3D surface plot\n');
fprintf('      "contourf"  - Filled contour plot (recommended for visualizing barriers in WDOT)\n');
fprintf('      "contour"   - Contour lines\n');
fprintf('      "contour3"  - 3D contour plot\n\n');

fprintf('   Optional parameters:\n');
fprintf('      barrier     - Barrier matrix (for WDOT problems)\n\n');

fprintf('   Examples:\n');
fprintf('      %% DOT: Default image view\n');
fprintf('      show_evolution_2d(output.rho, "imshow", "Evolution");\n\n');

fprintf('      %% WDOT: With barrier visualization\n');
fprintf('      show_evolution_2d(output.rho, "contourf", "WDOT Evolution", barrier);\n\n');
end

function show_movement()
fprintf('2. Density Movement: show_movement_2d(rho, Ex, Ey, title, barrier)\n\n');

fprintf('   Description:\n');
fprintf('      Displays density evolution with velocity vector field overlay.\n');
fprintf('      Arrows indicate transport direction and magnitude.\n\n');

fprintf('   Optional parameters:\n');
fprintf('      barrier     - Barrier matrix (for WDOT problems)\n\n');

fprintf('   Examples:\n');
fprintf('      %% DOT: Basic movement visualization\n');
fprintf('      show_movement_2d(output.rho, output.Ex, output.Ey, ...\n');
fprintf('                       "Transport vectors");\n\n');

fprintf('      %% WDOT: Movement with barrier\n');
fprintf('      show_movement_2d(output.rho, output.Ex, output.Ey, ...\n');
fprintf('                       "WDOT Movement", barrier);\n\n');
end

function show_convergence()
fprintf('3. Convergence Analysis: show_residualCurve(data, title, names, ...)\n\n');

fprintf('   Display KKT errors or primal-dual gap over iterations or time.\n\n');

fprintf('   Key parameters:\n');
fprintf('      ''xIteration'', iter  - Plot vs. iteration numbers\n');
fprintf('      ''xTime'', time       - Plot vs. computation time\n\n');

fprintf('   Examples:\n');
fprintf('      %% KKT errors vs. iterations\n');
fprintf('      show_residualCurve(runHist.kkt, "KKT errors", ...\n');
fprintf('                         runHist.kktNames, ...\n');
fprintf('                         ''xIteration'', runHist.iter);\n\n');

fprintf('      %% Convergence over time (all levels)\n');
fprintf('      show_residualCurve(runHistML.kkt, "KKT errors (all levels)", ...\n');
fprintf('                         runHistML.kktNames, ...\n');
fprintf('                         ''xTime'', runHistML.time);\n\n');

fprintf('      %% Primal-dual gap\n');
fprintf('      show_residualCurve(runHist.pdGap, "Primal-dual gap", [], ...\n');
fprintf('                         ''xIteration'', runHist.iter, ...\n');
fprintf('                         ''yLabel'', ''Relative primal dual gap'');\n\n');
end

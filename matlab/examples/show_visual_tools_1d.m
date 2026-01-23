function show_visual_tools_1d()
%SHOW_VISUAL_TOOLS_1D Display available visualization tools for 1D problems
%
%   show_visual_tools_1d() displays all available visualization
%   functions and their usage for analyzing DOT 1D results.
%
% SEE ALSO: show_evolution_1d, show_residualCurve

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('  Visualization Tools for 1D Problems (DOT)\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

show_evolution();
fprintf('────────────────────────────────────────────────────────────────\n\n');
show_diagnostics();
fprintf('────────────────────────────────────────────────────────────────\n\n');
show_convergence();

end

function show_evolution()
fprintf('1. Density Evolution: show_evolution_1d(rho, showFunc, figName)\n\n');

fprintf('   Visualization types (showFunc):\n');
fprintf('      "join"   - Multiple curves in one plot (default)\n');
fprintf('      "tile"   - Separate subplots in tiled layout\n\n');

fprintf('   Examples:\n');
fprintf('      %% Default: joined curves\n');
fprintf('      show_evolution_1d(rho, "join");\n\n');

fprintf('      %% Tiled layout\n');
fprintf('      show_evolution_1d(rho, "tile");\n\n');

fprintf('      %% Custom figure name\n');
fprintf('      show_evolution_1d(rho, "join", "Evolution: Gaussian");\n\n');
end

function show_diagnostics()
fprintf('2. Diagnostic Histograms:\n\n');

fprintf('   a) Negative density check:\n');
fprintf('      hist_negative_density(rho, title)\n\n');

fprintf('      Displays histogram of density values less than 0.\n');
fprintf('      Example:\n');
fprintf('          hist_negative_density(rho, "Histogram: density < 0");\n\n');

fprintf('   b) Constraint violation check:\n');
fprintf('      hist_violation_q_1d(q0, bx, title)\n\n');

fprintf('      Displays histogram of f(q) values greater than 0.\n');
fprintf('      Example:\n');
fprintf('          hist_violation_q_1d(q0, bx, "Histogram: f(q) > 0");\n\n');

fprintf('   c) Mass conservation check:\n');
fprintf('      check_massConservation_1d(rho)\n\n');

fprintf('      Prints mass at each time step to verify conservation.\n');
fprintf('      Example:\n');
fprintf('          check_massConservation_1d(rho);\n\n');
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

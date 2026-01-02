function [dotcomplem, normRho, norm_rhoFq, mRhoB, normM, normRhoB] = compute_kkt_dot_complement(q, alpha, z2, sigma, h, nt, nx, qInd, cScale, dScale, D, E)
    rhoT = (sigma * cScale * D) * alpha(1 : qInd.bx-1);
    rhoFq = rhoT + (dScale / D) * q(1 : qInd.bx-1) + sum(((dScale / E) * z2(:, 2:5)).^2, 2) / 4;
    rhoFq(rhoFq < 0) = 0;

    dotcomplem = normL2(rhoT - rhoFq, h);
    normRho = normL2(rhoT, h);
    norm_rhoFq  = normL2(rhoFq, h);

    rho = movmean(cat(2, zeros(nx, 1), reshape(rhoT, nx, nt-1), zeros(nx, 1)), 2, 2, "Endpoints", "discard");
    rhoBx = (dScale / D) * ( reshape(movmean(rho, 2, 1, "Endpoints", "discard"), [], 1) .* q(qInd.bx : end) );
    mx = (sigma * cScale * D) * alpha(qInd.bx : end);
    normM = normL2(mx, h);
    normRhoB = normL2(rhoBx, h);
    mRhoB      = normL2(mx - rhoBx, h);
end

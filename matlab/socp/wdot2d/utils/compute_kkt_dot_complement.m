function [dotcomplem, normRho, norm_rhoFq, mRhoB, normM, normRhoB] = compute_kkt_dot_complement(q, alpha, z2, sigma, h, nt, nx, ny, qInd, weight, cScale, dScale, D, E)
    Dalpha = weight .* alpha;
    rhoT = (sigma * cScale * D) * Dalpha(1 : qInd.bx-1);
    rhoFq = rhoT + (dScale / D) * q(1 : (nt-1)*ny*nx) + sum(((dScale / E) * z2(:, 2:9)).^2, 2) / 4;
    rhoFq(rhoFq < 0) = 0;
    
    dotcomplem = normL2(rhoT - rhoFq, h);
    normRho = normL2(rhoT, h);
    norm_rhoFq = normL2(rhoFq, h);

    rho = movmean(cat(3, zeros(ny, nx), reshape(rhoT, ny, nx, nt-1), zeros(ny, nx)), 2, 3, "Endpoints", "discard");
    rhoBx = (dScale / D) * ( reshape(movmean(rho, 2, 2, "Endpoints", "discard"), [], 1) .* q(qInd.bx : qInd.by-1) );
    rhoBy = (dScale / D) * ( reshape(movmean(rho, 2, 1, "Endpoints", "discard"), [], 1) .* q(qInd.by : end) );
    mx = (sigma * cScale * D) * Dalpha(qInd.bx : qInd.by-1);
    my = (sigma * cScale * D) * Dalpha(qInd.by : end);

    mRhoB = sqrt(normL2(mx - rhoBx, h)^2 + normL2(my - rhoBy, h)^2);
    normM = sqrt(normL2(mx, h)^2 + normL2(my, h)^2);
    normRhoB = sqrt(normL2(rhoBx, h)^2 + normL2(rhoBy, h)^2);
end

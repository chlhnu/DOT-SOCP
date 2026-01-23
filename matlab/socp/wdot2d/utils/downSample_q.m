function qR = downSample_q(nt, nx, ny, q)
%% Calculate the downsampling of q-like variables located on a staggered grid

[nt2, nx2, ny2] = deal((nt+1)/2, (nx+1)/2, (ny+1)/2);

ProlongT = kron(kron(gene_prolongMat1dim_nearest(nt2-1), gene_prolongMat1dim_linear(nx2)), gene_prolongMat1dim_linear(ny2));
ProlongX = kron(kron(gene_prolongMat1dim_linear(nt2), gene_prolongMat1dim_nearest(nx2-1)), gene_prolongMat1dim_linear(ny2));
ProlongY = kron(kron(gene_prolongMat1dim_linear(nt2), gene_prolongMat1dim_linear(nx2)), gene_prolongMat1dim_nearest(ny2-1));

RestriT = transpose(ProlongT ./ sum(ProlongT, 1));
RestriX = transpose(ProlongX ./ sum(ProlongX, 1));
RestriY = transpose(ProlongY ./ sum(ProlongY, 1));

bxInd = (nt-1)*nx*ny + 1;
byInd = bxInd + nt*(nx-1)*ny;
qR = cat(1, ...
    RestriT * q(1 : bxInd-1), ...
    RestriX * q(bxInd : byInd-1), ...
    RestriY * q(byInd : end));

end

%% Prolongation matrix in 1dim
function P = gene_prolongMat1dim_linear(nC)
    nR   = 2*(nC - 1) + 1;
    iVec = cat(2, 1:2:nR, repmat(2:2:nR-1, 1, 2));
    jVec = cat(2, 1:nC, 1:nC-1, 2:nC);
    vVec = cat(2, ones(1, nC), repmat(.5, 1, 2*(nC-1)));
    P    = sparse(iVec, jVec, vVec);
end

function P = gene_prolongMat1dim_nearest(nC)
    nR   = 2*nC;
    iVec = cat(2, 1:2:nR, 2:2:nR);
    jVec = cat(2, 1:nC, 1:nC);
    vVec = cat(2, ones(1, nC), ones(1, nC));
    P    = sparse(iVec, jVec, vVec);
end

function phic = downSample_phi(phi)
%% Downsample phi (1 dimension)

len  = length(phi) - 1;
lenc = len / 2;
phic = zeros(lenc+1, 1);

ind = 3 : 2 : len - 1;

phic(2:lenc) = 0.5 * phi(ind) + 0.25 * (phi(ind-1) + phi(ind+1));
phic(1) = (2 / 3) * phi(1) + (1 / 3) * phi(2);
phic(end) = (1 / 3) * phi(end-1) + (2 / 3) * phi(end);
                
end
    
function [] = show_movement_2d(rho, Ex, Ey, figName, barrier)
%% Show movement of rho

if ~exist('figName', 'var') || isempty(figName)
    figName = "Density movement";
end

if ~exist('barrier', 'var') || isempty(barrier)
    barrier = [];
end

[ny, nx, nt] = size(rho);

time_display = 3;
num1 = 40;
num2 = 30;

if nx <= num1
    skipx = 1;
else
    skipx = floor(nx/num2);
end

if ny <= num1
    skipy = 1;
else
    skipy = floor(nx/num2);
end

%% Plot
zeroRho = (rho < max(rho, [], 'all') / 1e2);
Ex(zeroRho) = 0;
Ey(zeroRho) = 0;

if ~ isempty(barrier)
    if all(size(barrier) == [ny, nx])
        barrier = repmat(barrier, [1, 1, nt]);
        rho(barrier) = Inf;
    else
        error("Argument at posistion 5 is invalid");
    end
end

% skip
indx = 1:skipx:nx;
indy = 1:skipy:ny;

% plot
[xx, yy] = meshgrid(indx, indy);
max_value = max(rho, [], "all");

figure("Name", figName);
for t = 1 : nt
    imshow(rho(:,:,t), [0, max_value], 'InitialMagnification', 512);
    hold on;
    quiver(xx, yy, Ex(indy,indx,t), Ey(indy,indx,t), 'b', 'LineWidth', 1);
    hold off;
    pause(time_display/nt);
end

end

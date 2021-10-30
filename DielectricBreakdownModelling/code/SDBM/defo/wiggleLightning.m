clear; clc;

[V, E] = readOBJ("sampleLightning2D.obj");

fig = figure('Position', [100, 100, 700, 600]);
hold on;
aggPlotS = plot( [V(E(:, 1), 1)'; V(E(:, 2), 1)'],  [V(E(:, 1), 2)'; V(E(:, 2), 2)'], 'LineWidth', 2, 'Color', 'black');
xlim([0, 60]);
ylim([0, 60]);
dt = 0.5;
vobj = figmp4("wiggle", []);
while(true)
    vel = rand(length(V), 3) - 0.5;
    vel(:, 3) = 0;
    
    V = V + vel*dt; 
    
    delete(aggPlotS);
    aggPlotS = plot( [V(E(:, 1), 1)'; V(E(:, 2), 1)'],  [V(E(:, 1), 2)'; V(E(:, 2), 2)'], 'LineWidth', 2, 'Color', 'black');
    drawnow;
    vobj = figmp4("wiggle", vobj);
end

vobj.close;
function refineBoltUI(V0, E0)
%REFINEBOLTUI Summary of this function goes here
%   Detailed explanation goes here
refineFig = figure;
rf = 4;
hold on;
set(gca,'Color','k');
set(gcf,'KeyPressFcn',@keyboard_down_callback);
xlim([0, 1]); ylim([0, 1]);
title("Bolt Refinement");

ogPlot = plot_edges(V0, E0, 'Color', [0.9 0.1 0.1]);
rPlot = [];
rV = V0;
rE = E0;
rPlot = plot_edges(rV, rE, 'Color', [1 1 1]);


    function keyboard_down_callback(src, ev)
        if(strcmp(ev.Character,'r'))
            rf = rf * 2;
            [rV, rE] = refineBolt(V0, E0, rf);
            delete(rPlot)
            rPlot = plot_edges(rV, rE, 'Color', [1 1 1]);
        elseif(strcmp(ev.Character, 'c'))
            rf = rf / 2;
            [rV, rE] = refineBolt(V0, E0, rf);
            delete(rPlot)
            rPlot = plot_edges(rV, rE, 'Color', [1 1 1]);        
        elseif(strcmp(ev.Key, 'space'))
            [rV, rE] = refineBolt(V0, E0, rf);
            delete(rPlot)
            rPlot = plot_edges(rV, rE, 'Color', [1 1 1]);
        end
        
    end
end


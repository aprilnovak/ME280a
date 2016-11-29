function [none] = PlotInTime(param, func, physical_domain, y_label, discr, num_steps)

current_max = max(param .* func{1,1});
current_min = min(param .* func{1,1});

for m = 2:num_steps
    maximum = max(param .* func{1,m});
    minimum = min(param .*func{1,m});
    
    if maximum > current_max
        current_max = maximum;
    end
    
    if minimum < current_min
        current_min = minimum;
    end
end

figure
plot(physical_domain, param .* func{1, 1})
ylim([minimum, maximum])
ylabel(y_label)
xlabel('Problem Domain')
hold on

for n = [2:discr:num_steps, num_steps]
    if n == num_steps
        plot(physical_domain, param .* func{1, n}, 'k-', 'LineWidth', 4)
    else
        plot(physical_domain, param .* func{1, n})
    end
    drawnow 
end

none = 0;
end


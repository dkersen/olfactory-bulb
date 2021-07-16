% Plots sample probability curves for an MC and GC

num_trials = 5;

trial_results = cell(2,num_trials);


for i = 1:num_trials

    gr = granule();
    gr = gr.assignProperties(1000);
    mi = mitral();
    mi = mi.assignProperties(1,50,50,1000,100);

    mi.z = gr.z0 + (gr.zmax-gr.z0)/2;
    syn = 0;
    gr.x = 0;
    gr.y = 0;

    mi.x = 0;
    max_distance = mi.radius + gr.calculateRadius(mi.z);
    distance_range = max_distance:-0.1:0;
    results = zeros(1,length(distance_range));

    for y = 1:length(distance_range)
        mi.y = distance_range(y);
        results(y) = 1-exp(-synProb(gr,mi,0, gr.calculateRadius(mi.z), distance_range(y)));
    end
    trial_results{1,i} = distance_range;
    trial_results{2,i} = results;
end

for i =1:num_trials
    plot(trial_results{1,i}, trial_results{2,i})
    hold on
end
ylim([0,1])
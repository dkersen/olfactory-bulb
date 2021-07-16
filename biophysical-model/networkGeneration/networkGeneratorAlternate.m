% Generates an alternate configuration of GCs given a pre-existing set of
% mitral cells

% maximum radius of the OB space, modeling the space as a thin cylinder
rmax = 600;

load('mitralCells.mat')
load('glomeruli.mat')

mitralNum = length(mitralArray);

% set the target number of GCs
granPerMit = 15; 
granuleNum = granPerMit * mitralNum;

% create the array to store GC objects
granuleArray = [];

% create an empty matrix to record connections (a connnection between an
% MC and GC is recorded as 1 in the relevant entry in the network matrix,
% is 0 otherwise)
network = zeros(mitralNum, granuleNum);

% GC index 
granuleindex = 0;

% distance matrix for synapses during simulation (
distance = zeros(mitralNum, granuleNum);


while length(granuleArray) < granuleNum
    
    % generate a new GC and assign properties
    newGranule = granule();
    newGranule = newGranule.assignProperties(rmax);
    
    % generate an empty matrix representing potential connections between
    % the MCs of the network and the new GC
    tempNet = zeros(mitralNum,1);
    tempDistance = ones(mitralNum,1)*-1;

    % for each MC, determine the probability of synapse with the new GC
    for j = 1:mitralNum
         
        % set the MC
        mitralCell = mitralArray(j);
        
        % calculate the total pre-existing synapses 
        deadSpace = sum(network(j,:));
        
        % determine the location of the center of the GC cone at the height
        % of the MC
        [g_loc_x, g_loc_y] = newGranule.calculateLocation(mitralCell.z);
        
        % calculate the radius of the GC cone at the height of the MC
        g_radius = newGranule.calculateRadius(mitralCell.z);
        
        % calculate the distance between the GC cone and MC at the height
        % of the MC
        s = norm([mitralCell.x, mitralCell.y] - [g_loc_x, g_loc_y]);
        
        % if the MC and GC are in range of one another, calculate the
        % average number of synapses
        if s < mitralCell.radius + g_radius && mitralCell.z > newGranule.z0 && mitralCell.z < newGranule.zmax
            
            % the average number of synapses for the MC/GC pair
            lambda = synProb(newGranule, mitralCell, deadSpace, g_radius, s);
            
            % calculate the probability of synapse assuming a Poisson
            % distribution
            prob = 1 - exp(-lambda);
            
            % sample the given probability and assign synaptic distance
            if rand < prob
                tempNet(j,1) = 1;
                while tempDistance(j,1) == -1
                   r = rand*mitralCell.radius;
                   theta = rand * 2 * pi;
                   x_r = r*cos(theta);
                   y_r = r*sin(theta);
                   if norm([x_r,y_r] - [s, 0]) < g_radius
                       tempDistance(j,1) = r;
                   end
                end       
            end
        end
    end
    
    % calculate the total number of synapses the GC has made and
    % incorporate the new GC if it is greater than 0
    totalSynapses = sum(tempNet(:,1));
   
    if totalSynapses > 0
        % if the total number of synapses is greater than the number of
        % available spines, select a random subset of MCs equal to the
        % number of available spines to be incorporated into the network
        if totalSynapses > newGranule.availableSpines
           synmitrals = find(tempNet);
           synmitrals = synmitrals(randperm(length(synmitrals)));
           tempNet(:) = 0;
           tempNet(synmitrals(1:newGranule.availableSpines)) = 1;
           tempDistance(synmitrals(newGranule.availableSpines+1:totalSynapses)) = -1;
        end
        % update the GC index
        granuleindex = granuleindex + 1;
        
        % assign the new connections to the appropriate entry in the
        % network for the GC
        network(:, granuleindex) = tempNet;
        distance(:,granuleindex) = tempDistance;
        
        % update the GC array
        granuleArray = [granuleArray newGranule];
    end
    
    disp(length(granuleArray));
end

% save arrays and network matrix
save('mitralCellsAlt.mat', 'mitralArray','-v7.3');
save('fullNetworkAlt.mat', 'network', '-v7.3');
save('granuleCellsAlt.mat', 'granuleArray', '-v7.3');
save('glomeruliAlt.mat', 'glomXYarray','-v7.3');
save('distanceAlt.mat', 'distance','-v7.3');

        
    

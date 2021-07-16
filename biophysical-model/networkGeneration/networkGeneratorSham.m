% Given a set of mitral cells, generates a network with
% distance-independent (i.e. constant) probability of synapse

load('mitralCells.mat')
load('glomeruli.mat')

mitralNum = length(mitralArray);

% maximum radius of the OB space, modeling the space as a thin cylinder
rmax = 600;


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
    
    prob = 0.0230;
    
    for j = 1:mitralNum
         
        % set the MC
        mitralCell = mitralArray(j);
        
        
        % calculate the distance between the GC cone and MC at the height
        % of the MC
        s = norm([mitralCell.x, mitralCell.y] - [newGranule.x, newGranule.y]);

        % sample the given probability and assign synaptic distance
        if rand < prob
            tempNet(j,1) = 1;
            tempDistance(j,1) = rand*s;
        end
    end
    
    % calculate the total number of synapses the GC has made and
    % incorporate the new GC if it is greater than 0
    totalSynapses = sum(tempNet(:,1));
   
    if totalSynapses > 0
        % if the total number of synapses is greater than the number of
        % available spines, select a random subset of MCs equal to the
        % number of available spines to be incorporated into the network

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
save('mitralCellsSham.mat', 'mitralArray','-v7.3');
save('fullNetworkSham.mat', 'network', '-v7.3');
save('granuleCellsSham.mat', 'granuleArray', '-v7.3');
save('distanceSham.mat', 'distance','-v7.3');
save('glomeruliSham.mat', 'glomeruli', 'glomXYarray','-v7.3');

        
    

% Examines how the number of MCs and cortical cells in the firing rate
% model affects the slope of the relationship between initial similarity
% and change in similarity at high threshold
% as well as the r^2 value of a linear fit 

feedback_type = 'positive positive'; %'negative negative', 'positive positive'

%Odor distribution
k_odor = -0.281; sigma_odor = 3.331; theta_odor = -0.4;
odor_dist = makedist('GeneralizedPareto','k',k_odor,'sigma',sigma_odor,'theta',theta_odor);
odor_dist = truncate(odor_dist, 0, 13);

% positive feedback distribution
pos_range = -1.5:0.001:8.5;
alpha =4.232;
omega =2.687;
xi = 0.7785;
pos_pdf = zeros(1,length(pos_range));
for i = 1:length(pos_range)
    pos_pdf(i) = skewed_gaussian(pos_range(i),alpha,omega,xi);
end
pi_zero = 0.69;

% negative feedback distribution
mu_neg = -0.1086;
sigma_neg = 0.4477;
neg_dist = makedist('Normal','mu', mu_neg,'sigma', sigma_neg);
neg_dist = truncate(neg_dist,-2.5,1.5);

mu_neg = 0.733; sigma_neg = 0.253;
%k_G = -0.118552; sigma_G = 0.787966; theta_G = -0.3;
G_dist = makedist('lognormal','mu',mu_neg,'sigma',sigma_neg);
G_dist = truncate(G_dist, 1, 4);


totalMit = 500:500:15000;
totalCor = 1000:1000:100000;
slopeMatrix = zeros(length(totalMit),length(totalCor));
rMatrix = zeros(length(totalMit),length(totalCor));

results = cell(length(totalMit),length(totalCor));
%Number of mitral cells per glomerulus
mitPerGlom = 20;
targetAct = 0.01;
numTrials = 1;

for mit = 1:length(totalMit)
    for cor= 1:length(totalCor)
        for trials = 1:numTrials

            mitralNum = totalMit(mit);

            %Number of glomeruli  
            glomNum = round(mitralNum/mitPerGlom);

            %Number of mitral cells

            %Number of cortical cells
            moduleNum = totalCor(cor);


            %Number of mitral cells which link to each cortical cell
            p_cortical = 0.07;
            mitPerMod = round(p_cortical*mitralNum);

            %create the module matrix
            module_matrix = zeros(moduleNum, mitPerMod);

            %assign MCs to modules
            for i = 1:moduleNum
                G = randperm(mitralNum);
                module_matrix(i,:) = G(1:mitPerMod)';
            end

            %Fraction of glomeruli targeted by an odor
            p_odor = 0.12;
            glomPerOdor = round(p_odor*glomNum);

            %Number of overlaps
            overlap_fraction = 0:0.05:1;
            num_overlap = length(overlap_fraction);

            %Fraction of MCs receiving positive feedback
            p_feedback = 0.08;

            %positive feedback overlap
            fb_overlap = 1;

            %Resample fraction
            resample_frac = 1;

            %Number of MCs targeted by an odor
            odorMit = mitPerGlom*glomPerOdor;

            mit_scramble = randperm(mitralNum);
            odor_responses = zeros(1,mitralNum);

            %create firing rate vectors
            R = zeros(moduleNum,2);
            RFB = zeros(moduleNum,2);

            for i = 1:length(overlap_fraction)

                % Create odor response vector
                odor_responses = random(odor_dist, 1, mitralNum);
                odor_responses_2 = random(odor_dist, 1, mitralNum);


                %Assign odor vectors
                odor = zeros(2,mitralNum);
                odor_FB = zeros(2,mitralNum);
                odorReceivers1 = zeros(1,mitralNum);
                odorReceivers2 = zeros(1,mitralNum);

                front_glom = randi(glomNum-2*glomPerOdor+1)-1;
                front = round(front_glom/glomNum * mitralNum + 1);
                back = round(front + odorMit-1);
                glomRange = front:back;
                odor(1,glomRange) = odor(1,glomRange) + odor_responses(glomRange);
                odorReceivers1(glomRange) = 1;

                front_glom = front_glom + round((1-overlap_fraction(i))*glomPerOdor);
                front_2 = round(front_glom/glomNum * mitralNum + 1);
                back_2 = round(front_2 + odorMit-1);
                glomRange2 = front_2:back_2;
                odor(2,glomRange2) = odor(2,glomRange2) + odor_responses_2(glomRange2);
                odorReceivers2(glomRange2) = 1;


                %Add feedback
                odor_FB(1,:) = odor(1,:); 
                odor_FB(2,:) = odor(2,:);

                if strcmp(feedback_type, 'positive positive') == true %feedback targets same MCs with same positive strength

                    FB_num = round(p_feedback*mitralNum);
                    feedback = randpdf(pos_pdf,pos_range,[1,mitralNum]);
                    zero_sample = rand(1,mitralNum);
                    indices = find(zero_sample < pi_zero);
                    feedback(indices) = 0;

                    FB_set_1 = mit_scramble(1:FB_num);
                    feedback_1 = zeros(1,mitralNum);
                    feedback_1(FB_set_1) = feedback(FB_set_1);

                    overlap_num = round((1-fb_overlap)*FB_num);
                    FB_set_2 = mit_scramble(1+overlap_num:FB_num+overlap_num);
                    feedback_2 = zeros(1,mitralNum);
                    feedback_2(FB_set_2) = feedback(FB_set_2);

                    odor_FB(1,:) = odor_FB(1,:) + feedback_1;
                    odor_FB(2,:) = odor_FB(2,:) + feedback_2;

                elseif strcmp(feedback_type, 'negative negative') == true

                    O_vals = find(odor_responses ~= 0);    
                    feedback_1(O_vals) = random(neg_dist, 1, length(O_vals));

                    O_vals = find(odor_responses_2 ~= 0);
                    feedback_2(O_vals) = random(neg_dist, 1, length(O_vals));

                    odor_FB(1,:) = poslin(odor_FB(1,:) + feedback_1);
                    odor_FB(2,:) = poslin(odor_FB(2,:) + feedback_2);     
                end


                %assign sum of MC firing rate set to each module entry in the firing rate vector
                R = zeros(moduleNum,2);
                RFB = zeros(moduleNum,2);

                for x = 1:2
                    for j = 1:moduleNum
                        R(j,x) = sum(odor(x,module_matrix(j,:)));
                        RFB(j,x) = sum(odor_FB(x,module_matrix(j,:)));
                    end
                end

                %mean-center the inputs (argue that effective inputs are the same)
                k_a = 1;
                R = R - k_a*mean(R);
                RFB = RFB - k_a*mean(RFB);


                %set threshold range
                if i == 1
                    minimum =  min(round(min(min(RFB))),round(min(min(R))));
                    maximum = max(round(max(max(RFB))),round(max(max(R))))+50;
                    threshold = minimum:1:maximum;
                    activation = zeros(num_overlap,length(threshold));
                    for t = 1:length(threshold)
                        Z1 = sigmoid(R(:,1),threshold(t));
                        activation(1,t) = sum(Z1)/moduleNum;  
                        if activation(1,t) < targetAct
                            break
                        end
                    end
                    activation(activation == 0) = [];
                    [I,J] = min(abs(activation - targetAct));
                    threshold = threshold(J);
                    overlap = zeros(num_overlap,length(threshold));
                    overlap_FB = zeros(num_overlap,length(threshold));
                    slope = zeros(num_overlap, length(threshold));
                end

                for t = 1:length(threshold)
                    Z1 = sigmoid(R(:,1),threshold(t));
                    Z2 = sigmoid(R(:,2),threshold(t));                           
                    Z1f = sigmoid(RFB(:,1),threshold(t));
                    Z2f = sigmoid(RFB(:,2),threshold(t));

                    overlap(i,t) = cosine_distance(Z1,Z2); 
                    overlap_FB(i,t) = cosine_distance(Z1f,Z2f);
                end
            end

            overlap(isnan(overlap)) = 0;
            overlap_FB(isnan(overlap_FB)) = 0;
            coefs = polyfit(overlap,overlap_FB-overlap,1);
            slopeMatrix(mit,cor) = slopeMatrix(mit,cor)+ coefs(1);
            [C,~] = rsquared(overlap_FB - overlap, coefs(1)*overlap + coefs(2),2);
            rMatrix(mit,cor) = rMatrix(mit,cor) + C;
        end
        slopeMatrix(mit,cor) = slopeMatrix(mit,cor)/numTrials;
        rMatrix(mit,cor) = rMatrix(mit,cor)/numTrials;
        save('NEresults','slopeMatrix', 'rMatrix')
    end
    disp(mit)
end
        

%save('sigmoidsamepos.mat','activation','overlap','overlap_FB','threshold','feedback_corr','feedback_cosdis')
   
function output = sigmoid(rate_vector, theta)
    sigma = 0.5; % value I'm using is 0.5
    output = 1./(1+exp(-sigma*(rate_vector-theta)));
end

function dist = cosine_distance(vec_1, vec_2)
    dist = sum(vec_1 .* vec_2)/(norm(vec_1)*norm(vec_2)); 
end





    
    
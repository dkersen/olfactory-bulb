
% Plots the relationship between initial similarity and change in similarity for different feedback conditions 


load('samePos.mat')
t=141;
%t = 111-28; %194 for high, 141 for low (0.01, 0.3)
disp(threshold(t))
curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
hold on
disp(feedback_cosdis(1))
pause(1)
% 
%load('sigmoiddiffposmid.mat')

load('samePosLite.mat')
t=142;
%t = 111-28; %195 for high, 142 for low
disp(threshold(t))
curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
disp(feedback_cosdis(1))
pause(1)


load('diffPos.mat')
t=156;
%t = 111-28; %209 for high, 156 for low
disp(threshold(t))
curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
disp(feedback_cosdis(1))
pause(1)

load('diffPosZero.mat')
t=147;
%t = 111-28; %200 for high, 147 for low
disp(threshold(t))
curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
disp(feedback_cosdis(1))
pause(1)


% 
% load('negPos.mat')
% t=182;
% %t = 111-28; %182 for high, 129 for low
% disp(threshold(t))
% curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
% plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
% disp(feedback_cosdis(1))
% pause(1)
% 
% 
% load('diffNeg.mat')
% t=176;
% %t = 111-28; %176 for high, 125 for low
% disp(threshold(t))
% curvefit = fit(overlap(:,t),overlap_FB(:,t) - overlap(:,t),'poly1');
% plot(curvefit, overlap(:,t),overlap_FB(:,t) - overlap(:,t))
% disp(feedback_cosdis(1))
% pause(1)






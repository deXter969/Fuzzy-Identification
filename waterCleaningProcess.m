close all; clear; clc;

      % TUB1 PH TE AL PAC TUB2
data = [10 7.1 18.8 53 1300 1;
        17 7 18.6 50 1300 1;
        22 7.3 19.4 46 1400 2;
        50 7.1 19.5 40 1400 1;
        9 7.3 23.3 48 900 4;
        11 7.1 20.7 50 900 1;
        12 7.2 21.3 50 900 3;
        14 7.2 23.6 53 900 4;
        35 7 17.8 35 1200 1;
        20 7 16.6 40 1100 1;
        20 6.9 17.8 42 1100 1;
        18 7.1 17.3 40 1100 1;
        12 7.2 18.8 55 900 3;
        8 7.2 18 50 1000 1.5;
        11 7.1 19.2 49 1000 2;
        50 7 18 37 1200 1.5;
        35 7 17.7 42 1200 1.5;
        30 7 17.3 41 1100 1.5;
        16 7.1 19.3 42 1100 3];
    
input = [data(:,1:4) data(:,6)];
output = data(:,5);
PH = input(:,2);
AL = input(:,4);
TE = input(:,3);
    
k = size(data,2) - 1;   % Number of variables
m = size(data,1);       % Number of data points
n = 8;                  % Maximum number of rules
    
TUB1min = min(data(:,1)); TUB1max = max(data(:,1));
PHmin = min(data(:,2)); PHmax = max(data(:,2));
TEmin = min(data(:,3)); TEmax = max(data(:,3));
ALmin = min(data(:,4)); ALmax = max(data(:,4));
PACmin = min(data(:,5)); PACmax = max(data(:,5));
TUB2min = min(data(:,6)); TUB2max = max(data(:,6));

mPH = 1/(PHmin-PHmax); cPH = PHmin/(PHmin-PHmax); cPH1 = PHmax/(PHmax-PHmin);
mAL = 1/(ALmin-ALmax); cAL = ALmin/(ALmin-ALmax); cAL1 = ALmax/(ALmax-ALmin);
mTE = 1/(TEmin-TEmax); cTE = TEmin/(TEmin-TEmax); cTE1 = TEmax/(TEmax-TEmin);

%% Model 3-k

beta = zeros(m,n);
b = beta;

    for i = 1:n
        for j = 1:m
            smallPH(j) = -mPH*PH(j)+cPH1;
            smallAL(j) = -mAL*AL(j)+cAL1;
            smallTE(j) = -mPH*TE(j)+cTE1;
            bigPH(j) = mPH*PH(j)+cPH;
            bigAL(j) = mAL*AL(j)+cAL;
            bigTE(j) = mTE*TE(j)+cTE;
            switch i
                case 1
                    b(j,i) = min([smallPH(j) smallAL(j) smallTE(j)]);
                case 2
                    b(j,i) = min([smallPH(j) smallAL(j) bigTE(j)]);
                case 3
                    b(j,i) = min([smallPH(j) bigAL(j) smallTE(j)]);
                case 4
                    b(j,i) = min([smallPH(j) bigAL(j) bigTE(j)]);
                case 5
                    b(j,i) = min([bigPH(j) smallAL(j) smallTE(j)]);
                case 6
                    b(j,i) = min([bigPH(j) smallAL(j) bigTE(j)]);
                case 7
                    b(j,i) = min([bigPH(j) bigAL(j) smallTE(j)]);
                case 8
                    b(j,i) = min([bigPH(j) bigAL(j) bigTE(j)]);
            end
        end
        beta(:,i) = b(:,i)./sum(b(:,i));
    end
    
    X = beta;
    for i = 1:n
        X = [X repmat(beta(:,i),1,k).*input];
    end
    
    P = pinv(X'*X)*X'*output;
    
    p = P(1:n)';
    for i = 1:k
        p = [p; P(n*i+1:n*(i+1))'];
    end
    
    x = [ones(m,1) input];
    y = x*p;
    y = b.*y;
    Y = [];
    for i = 1:m
        Y = [Y; sum(y(i,:))./sum(b(i,:))];
    end
    Y = X*P;
    PI = sqrt(sum((Y - output).^2)/m)
    
    % Statistical model
    % PAC = 9.11*sqrt(TB1) - 79.8*PH + 12.7*CL + 1255.6
    Ystat = [994.7; 995.9; 1119.6; 1151.1; 1409.4; 1066.4; 1068.9; 1012.3; 1286.8; 1246.8; 1151.4; 1199.5; 1159.4; 985.7; 1009.3; 1038.2; 1398.3; 1290.6; 1038.5];
    PIstat = sqrt(sum((Ystat - output).^2)/m)
    
    %Plot
    figure(1);
    scatter(output,Y);
    grid on;
    hold on;
    plot([500 1500],[500 1500]);
    xlabel('Operator');
    ylabel('Model');
    
  
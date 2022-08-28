close all; clear; clc;

    %   HP   Mn   SG  MA
data = [93.9 1.11 0.3 14;
        94 14.52 -0.08 135;
        93.1 5.22 -0.05 54;
        93.6 5.36 -0.43 53;
        85.6 13.56 0.11 129;
        85.8 11.55 0.12 106;
        86.4 11.88 -0.06 113;
        93.1 13.85 0.16 119;
        88.5 12.27 -0.02 129;
        93 11.69 -0.03 127;
        93 9.86 -0.01 98;
        95 4.96 -0.13 46;
        94.7 2.19 0.01 25;
        85 1.89 -0.1 23;
        87.5 14.26 -0.09 121;
        90.5 11.52 -0.1 112;
        90.1 12.59 -0.06 116;
        90.2 5.67 -0.28 60];

input = [data(:,1) data(:,3:4)];
output = data(:,2);

HP = input(:,1);
SG = input(:,2);
MA = input(:,3);

k = size(data,2) - 1;   % Number of variables
m = size(data,1);       % Number of data points
n = 4;                  % Maximum number of rules

HPmin = min(data(:,1)); HPmax = max(data(:,1));
Mnmin = min(data(:,2)); Mnmax = max(data(:,2));
SGmin = min(data(:,3)); SGmax = max(data(:,3));
MAmin = min(data(:,4)); MAmax = max(data(:,4));

mHP = 1/(HPmin-HPmax); cHP = HPmin/(HPmin-HPmax); cHP1 = HPmax/(HPmax-HPmin);
mSG = 1/(SGmin-SGmax); cSG = SGmin/(SGmin-SGmax); cSG1 = SGmax/(SGmax-SGmin);
mMA = 1/(MAmin-MAmax); cMA = MAmin/(MAmin-MAmax); cMA1 = MAmax/(MAmax-MAmin);


%% MODEL 2-1
beta = zeros(m,n);
b = beta;

    for i = 1:n
        for j = 1:m
            smallHP(j) = -mHP*HP(j)+cHP1;
            smallSG(j) = -mSG*SG(j)+cSG1;
            smallMA(j) = -mMA*MA(j)+cMA1;
            bigHP(j) = mHP*HP(j)+cHP;
            bigSG(j) = mSG*SG(j)+cSG;
            bigMA(j) = mMA*MA(j)+cMA;
            switch i
                case 1
                    b(j,i) = min([smallHP(j) smallSG(j)]);
                case 2
                    b(j,i) = min([smallHP(j) smallSG(j)]);
                case 3
                    b(j,i) = min([smallHP(j) bigSG(j)]);
                case 4
                    b(j,i) = min([smallHP(j) bigSG(j)]);
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
    PI1 = sqrt(sum((Y - output).^2)/m);

    %% MODEL 2-2
    
    beta = zeros(m,n);
b = beta;

    for i = 1:n
        for j = 1:m
            smallHP(j) = -mHP*HP(j)+cHP1;
            smallSG(j) = -mSG*SG(j)+cSG1;
            smallMA(j) = -mMA*MA(j)+cMA1;
            bigHP(j) = mHP*HP(j)+cHP;
            bigSG(j) = mSG*SG(j)+cSG;
            bigMA(j) = mMA*MA(j)+cMA;
            switch i
                case 1
                    b(j,i) = min([smallSG(j) smallMA(j)]);
                case 2
                    b(j,i) = min([smallSG(j) bigMA(j)]);
                case 3
                    b(j,i) = min([bigSG(j) smallMA(j)]);
                case 4
                    b(j,i) = min([bigSG(j) bigMA(j)]);
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
    PI2 = sqrt(sum((Y - output).^2)/m);
    
    %% OUTPUT
    
    if PI1 < PI2
        fprintf('Stable: Model 2-2\nVariables: HP & SG\n');
    else
        fprintf('Stable: Model 2-2\nVariables: SG & MA\n');
    end
    
    %Plot
    figure(1);
    scatter(output,Y);
    grid on;
    hold on;
    plot([0 20],[0 20]);
    xlabel('Operator');
    ylabel('Model');

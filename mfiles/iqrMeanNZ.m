function [ iqrM, iqrStd, qtRange ] = iqrMeanNZ( X )
    % x can either be a vector of length M (column, preferrably)
    %       or a MxN matrix
    % output is a single value (if input is vector) or vector of length N
    [nS nT] = size(X);
    if nS == 1
        X = X';
        [nS nT] = size(X);
    end
   
    iqrM = zeros(nT,1);
    iqrStd = iqrM;
    for i=1:nT
        x = X(:,i);
        x(x<0) = [];
        qtRange = quantile(x,[.25 .75]);
        xMask = x>qtRange(1) & x<qtRange(2);
        iqrM(i) = mean(x(xMask));
        iqrStd(i) = std(x(xMask));
    end
end


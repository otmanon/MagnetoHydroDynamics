
function I=weightedSampling( W, numSamples)
% Given a vector of weights, sample indeces acording to corresponding weight
% W. Inspired by sampling mesh in gptoolbox by Alec Jacobson.
%Returns I, samples indices, where
    
    W0 = W./sum(W);
    zeroStart = zeros(1, size(W0, 2));
    W0 = [zeroStart; W0];
    
    C = cumsum(W0);
    
    R = rand(numSamples, size(W0, 2));
    
    I = zeros(size(W0, 2), numSamples);
    for i=1:size(W0, 2);
        [~, ind] = histc(R(:, i), C(:, i));    %Sees which "bins" the random samples fall into...
        I(i, :) = ind;                                %The larger the weight the larger the bin.
                                        %QED : weighted sampling.
    end;

end




  

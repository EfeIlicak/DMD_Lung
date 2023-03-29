function img = reconstructFreqImage(modeAmp,FreqComps,indxToUse)
% function img = reconstructFreqImage(modeAmp,modeFreqComp,indxToUse)
%
%
% Computes the reconstruction from given indices and their weights.
%
% [This is Equation 6 in the manuscript]
% Efe Ilicak, 30/10/2022.

[sx,sy,~] = size(FreqComps);

img_Components = zeros(sx,sy,length(indxToUse));
for indx = 1:length(indxToUse)
    img_Components(:,:,indx) = modeAmp(indxToUse(indx))*FreqComps(:,:,indxToUse(indx));
end
img = 2*abs(sum(img_Components,3));

% figure,imshow3(img)
end
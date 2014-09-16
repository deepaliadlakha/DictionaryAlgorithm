function [dictionaryPic]=display_dictionary(D,p,numDisplay,dictsize)

    dictionaryPic=zeros(p*floor(dictsize/numDisplay),numDisplay*p);
    for i=1:dictsize,
        dictionaryPic(floor((i-1)/numDisplay)*p+1:floor((i-1)/numDisplay)*p+p,int32(mod(i,numDisplay))*p+1:int32(mod(i,numDisplay))*p+p)=reshape(D(:,i),p,p);
    end
    imagesc (dictionaryPic); colorbar; axis equal tight;

end
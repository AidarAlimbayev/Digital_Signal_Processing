function img = edu_imgpreprocess6(I)

Igray = rgb2gray(I);

Ibw = im2bw(Igray,graythresh(Igray));

Iedge = edge(uint8(Ibw));

se = strel('square',3);

Iedge2 = imdilate(Iedge, se); 

Ifill= imfill(Iedge2,'holes');

[Ilabel num] = bwlabel(Ifill);
num
Iprops = regionprops(Ilabel);

Ibox = [Iprops.BoundingBox];
%Ibox = reshape(Ibox,[4 50]);
Ibox = reshape(Ibox,[4 66]);

Ic = [Iprops.Centroid];
%Ic = reshape(Ic,[2 50]);
Ic = reshape(Ic,[2 66]);
Ic = Ic';
Ic(:,3) = (mean(Ic.^2,2)).^(1/2);
%Ic(:,4) = 1:50;
Ic(:,4) = 1:66;

% Extra lines compare to example2 to extract all the components into an
% cell array
Ic2 = sortrows(Ic,2);

%for cnt = 1:5
for cnt = 1:6    
    Ic2((cnt-1)*11+1:cnt*11,:) = sortrows(Ic2((cnt-1)*11+1:cnt*11,:),4);
end

Ic3 = Ic2(:,1:2);
ind = Ic2(:,4);

%for cnt = 1:50
for cnt = 1:66
    img{cnt} = imcrop(Ibw,Ibox(:,ind(cnt)));
end




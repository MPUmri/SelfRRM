% Display (=affichage in French) of a NMF solution, for image datasets
%
% a = affichage(V,lig,Li,Co)
%
% Input.
%   V              : (m x r) matrix whose colums contains vectorized images
%  lig             : number of images per row in the display
%   (Co,Li)        : dimensions of images
%   bw             : if bw=1, reverse gray level
%
% Output.
%   Diplay columns of matrix V as images

function [a, Vaff] = affichage(V,lig,Li,Co,bw)

V = max(V,0); [m,r] = size(V); 
for i = 1 : r
    if max(V(:,i)) > 0
        V(:,i) = V(:,i)/max(V(:,i));
    end
end
Vaff = []; 
for i = 1 : r
        ligne = floor((i-1)/lig)+1; col = i - (ligne-1)*lig;
        Vaff((ligne-1)*Li+1:ligne*Li,(col-1)*Co+1:col*Co) = reshape(V(:,i),Li,Co);
end
[m,n] = size(Vaff);
for i = 1 : n/Co-1
        Vaff = [Vaff(:,1:Co*i+i-1) 0.5*ones(m,1) Vaff(:,Co*i+i:end)];
end
[m,n] = size(Vaff);
for i = 1 : m/Li-1
        Vaff = [Vaff(1:Li*i+i-1,:); 0.5*ones(1,n); Vaff(Li*i+i:end,:)];
end
warning('off'); figure; 
if nargin == 5 && bw == 1
    imshow(Vaff,[0 1]); 
else
    imshow(1-Vaff,[0 1]); 
end
colormap(gray); 
warning('on');
a = 1;
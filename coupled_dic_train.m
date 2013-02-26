function [Dh, Dl] = coupled_dic_train(Xh, Xl, codebook_size, lambda)
% Function for k-means clustering based on the codebook size

hDim = size(Xh, 1);
lDim = size(Xl, 1);

% joint learning of the dictionary
X = [1/sqrt(hDim)*Xh; 1/sqrt(lDim)*Xl];
Xnorm = sqrt(sum(X.^2, 1));

clear Xh Xl;

X = X(:, Xnorm > 1e-5);
X = X./repmat(sqrt(sum(X.^2, 1)), hDim+lDim, 1);

fprintf('size X\n');
size(X)

[Dtemp,sse] = vgg_kmeans(X, codebook_size);

D = Dtemp;
fprintf('size D\n');
size(D)

Dh = D(1:hDim, :);
Dl = D(hDim+1:end, :);

% normalize the dictionary
Dh = Dh./repmat(sqrt(sum(Dh.^2, 1)), hDim, 1);
Dl = Dl./repmat(sqrt(sum(Dl.^2, 1)), lDim, 1);

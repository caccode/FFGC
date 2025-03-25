function D = pdist2( X, Y )
% Calculates the distance between sets of vectors.
%
% Let X be an m-by-p matrix representing m points in p-dimensional space
% and Y be an n-by-p matrix representing another set of points in the same
% space. This function computes the m-by-n distance matrix D where D(i,j)
% is the distance between X(i,:) and Y(j,:).  This function has been
% optimized where possible, with most of the distance computations
% requiring few or no loops.

Yt = Y';
XX = sum(X.*X,2);
YY = sum(Yt.*Yt,1);
D = bsxfun(@plus,XX,YY)-2*X*Yt;
end

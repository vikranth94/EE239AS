function plotContour(mu,Sigma,colour)
    x1 = 0:.1:20; x2 = 0:.1:20;
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    contour(X1,X2,F,[0.007 0.007],colour)
end
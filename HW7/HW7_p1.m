%% EE239AS HW 7
clc
clear
close all
load('hw7_data.mat')

%% Problem 1
%% Part A: Eigenvalue Spectrum

Y = Xplan';

Y0 = bsxfun(@minus, Y, mean(Y,2));
S = (1/size(Y0,2))*(Y0*Y0');

[U,V] = eig(S);

V_vec = sqrt(diag(V));
[V_sort,I] = sort(V_vec,'descend');

dim_plot = 1:length(V_sort);

figure(1)
plot(dim_plot, V_sort)
title('Square Rooted Eigenvalue Spectrum')
xlabel('Dimension')
ylabel('Square Rooted Eigenvalues')

[U_D, V_D] = eigs(S, 3);
per_1 = 100*V_D(1,1)/sum(diag(V));
per_2 = 100*V_D(2,2)/sum(diag(V));
per_3 = 100*V_D(3,3)/sum(diag(V));

per = per_1+per_2+per_3;

fprintf('\nPart A: \n\nOverall Variance Captured By First Eigenvalue: %2.2f\n', per_1)
fprintf('\nOverall Variance Captured By Second Eigenvalue: %2.2f\n', per_2)
fprintf('\nOverall Variance Captured By Third Eigenvalue: %2.2f\n', per_3)
fprintf('\nOverall Variance Captured By First Three Eigenvalues: %2.2f\n', per)

%% Part B: 3D Representation 

Y_top = Y0'*U_D;

figure(2)
plot3(Y_top(1:91,1),Y_top(1:91,2),Y_top(1:91,3),'*',...
    Y_top(92:182,1),Y_top(92:182,2),Y_top(92:182,3),'*',...
    Y_top(183:273,1),Y_top(183:273,2),Y_top(183:273,3),'*',...
    Y_top(274:364,1),Y_top(274:364,2),Y_top(274:364,3),'*',...
    Y_top(365:455,1),Y_top(365:455,2),Y_top(365:455,3),'*',...
    Y_top(456:546,1),Y_top(456:546,2),Y_top(456:546,3),'*',...
    Y_top(547:637,1),Y_top(547:637,2),Y_top(547:637,3),'*',...
    Y_top(638:end,1),Y_top(638:end,2),Y_top(638:end,3),'k*')
title('Part B: Data Projected Onto 3D Space')
xlabel('Dimension 1')
ylabel('Dimension 2')
zlabel('Dimension 3')

%% Part C: Eigenvector Visualization

% U_D calculated in Part A

figure(3)
imagesc(U_D')
colorbar
title('Part A: Top Three Eigenvectors')
xlabel('Electrode')
ylabel('Principle Component')
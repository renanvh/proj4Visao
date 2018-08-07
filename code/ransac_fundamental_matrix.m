% Written by Henry Hu for CSCI 1430 @ Brown and CS 4495/6476 @ Georgia Tech

% Find the best fundamental matrix using RANSAC on potentially matching
% points

% 'matches_a' and 'matches_b' are the Nx2 coordnewInAtes of the possibly
% matching points from pic_a and pic_b. Each row is a correspondence (e.g.
% row 42 of matches_a is a point that corresponds to row 42 of matches_b.

% 'Best_Fmatrix' is the 3x3 fundamental matrix
% 'inliers_a' and 'inliers_b' are the Mx2 corresponding points (some subset
% of 'matches_a' and 'matches_b') that are inliers with respect to
% Best_Fmatrix.

% For this section, use RANSAC to find the best fundamental matrix by
% randomly sample interest points. You would reuse
% estimate_fundamental_matrix() from part 2 of this assignment.

% If you are trying to produce an uncluttered visualization of epipolar
% lines, you may want to return no more than 30 points for either left or

% RANSACFUNCTION

function [Best_Fmatrix, inliers_a, inliers_b] = ransac_fundamental_matrix(matches_a, matches_b)


%%%%%%%%%%%%%%%%
% Your code here
%%%%%%%%%%%%%%%%

%num iteracoes
numIter = 5000;

%auxiliar para carregar o melhor inlier obtido a cada iteracao
aux = 0;

%limite de distancia especificado pelo usuario
limit = 0.05;

sizeA = size(matches_a, 1);
sizeB = size(matches_b, 1);

%melhor matrix, a que obter a maior quantidade de inliers
Best_Fmatrix = [];
inliers_a = [];
inliers_b = [];

    for i = 1 : numIter
	      %coleta as amostras, sizeB ou sizeA podem ser utilizados, 10 amostras estao sendo coletadas
        n = randsample(sizeB, 10);
        pointA = matches_a(n, :);
        pointB = matches_b(n, :);
        
        %estima a matrix fundamental das amostras coletadas
        F = estimate_fundamental_matrix(pointA, pointB);
    
        newInA = [];
        newInb = [];
        
        for j = 1 : sizeB
            pointA = [matches_a(j, :) 1];
            pointB = [matches_b(j, :) 1];
            
            %calculando a metrica de distancia
            err = pointB * F * pointA';
            
            %se o valor obtido for menor que o limite especificado, os inliers sao atualizados
            if (abs(err) < limit)
                newInA = [newInA; pointA(:, 1:2)];
                newInb = [newInb; pointB(:, 1:2)];
            end
        end

        %seleciona os melhores inliers e a matrix fundamental
        if (size(newInA, 1) > aux)
            inliers_a = newInA;
            inliers_b = newInb;
            Best_Fmatrix = F;
            aux = size(newInA, 1);
        end
    
    end
    
end

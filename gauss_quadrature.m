function [gaussian_points, weights] =  gauss_quadrature(num_points) 
%**************************************************************************
% Returns the weight and location of the gaussian points using which
% numerical integral of an definite integral can be caluclated. 3D case.
%**************************************************************************
%
% Input parameter:
% num_points      - Number of gaussian points to be used. It is assumed that
%                   in all the directions same number of points will be 
%                   present.
%
% Output:
% gaussian_points - Location of the gaussian points for a given number of
%                   gaussian points. This gives a matrix of size m x 3
%                   where m is according to the number of gaussian points
%                   integral.
% weights         - Weights at all those gaussian points. By adding up the 
%                   function values multiplied with the wieghts at all the 
%                   gaussian points will give the definite integral value.
%

    switch num_points
        case 1
            weights = 8; % 2 x 2 x 2
            gaussian_points = [];  
        case 2
            weights = ones(1,8);
            temp = 1 / (sqrt(3));
            % Simply a 2 x 2 x 2 matrix with zeta, eta and nu varying as 
            % +temp or -temp. There combination will be our gaussian points 
            % or the 
            
            gaussian_points = [
                -temp, -temp, -temp;
                 temp, -temp, -temp;
                -temp,  temp, -temp;
                 temp,  temp, -temp;
                -temp, -temp,  temp;
                 temp, -temp,  temp;
                -temp,  temp,  temp;
                 temp,  temp,  temp;
                 ];
        case 3
            %  It will be a 27 x 1 matrix containing the weight of the
            %  gauusian points in case of 3-points 3D gaussian quadrature.
            weights = [(8/9)^3,];
            % 3 x 3 x 3 number of elements in the matrix
            
            % @TODO In case if use of 3-points gaussian quadrature comes
            % complete it and use this.
            gaussian_points = [];
    end
end
% Example 2.5
% Find the angle between two short signals represented as vectors. Give the angle in deg. The signals are:

% x  = [1.7, 3, 2.2] and y = [2.6, 1.6, 3.2]

% Since these signals are short, plot the two vector representations in three dimensions.

%  Solution
% Construct the two vectors and find the scalar product using Equation 2.29. When multiplying 
% the two vectors, be sure to use MATLAB’s .* operator to implement a point-by-point multiplica
% tion. 
clear all; 
close all; 
clc;

% Example 2.5 Find the angle between two vectors
x = [1.7, 3, 2.2];                  % vetor x
y = [2.6, 1.6, 3.2];                % vetor y   

sp = sum(x.*y);                     % calcula o produto escalar entre os vetores x e y
mag_x = sqrt(sum(x.^2));            % calcula a magnitude do vetor x       
mag_y = sqrt(sum(y.^2));            % calcula a magnitude do vetor y
cos_theta = sp/(mag_x*mag_y);       % calcula o cosseno do ângulo entre os vetores x e y (cosseno de theta) 
angle = acos(cos_theta);            % calcula o ângulo em radianos
angle = angle*360/(2*pi);           % converte o ângulo para graus

hold on;
plot3(x(1),x(2),x(3),'k*');                        % Plota o ponto representando o vetor x             
plot3([0 x(1)],[0 x(2)],[0 x(3)]);                 % Plota a linha representando o vetor x            
plot3(y(1),y(2),y(3),'k*');                        % Plota o ponto representando o vetor y             
plot3([0 y(1)],[0 y(2)],[0 y(3)]);                 % Plota a linha representando o vetor y             
title(['Angle = ', num2str(angle,2),' (deg)']);     % o título do gráfico mostra o ângulo entre os vetores x e y
xlabel('x-axis'); 
ylabel('y-axis'); 
zlabel('z-axis');               
grid on;

% The plot representing the two vectors. 
% The angle between the two vectors is calculated to be 26°. If these were signals, 
% the fairly small angle would indicate some correlation between them.
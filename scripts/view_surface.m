clc; clear; close all;

% T(r) = A*ln(r) + B

r_1 = 0.15;
r_2 = 0.45;
c_x = 0.5;
c_y = 0.5;

T_1 = 1;
T_2 = 200;

inner_neumann = true;

syms A B

if inner_neumann
    eqns = [A/r_1 == T_1, A*log(r_2) + B == T_2];
else
    eqns = [A*log(r_1) + B == T_1, A*log(r_2) + B == T_2];
end

S = vpasolve(eqns, [A B]);

A_sol = S.A;
B_sol = S.B;

phi_data = load("data_export.csv");
% Remove unwanted datapoints
phi_data(phi_data==0) = nan;

h = 1.0/(width(phi_data)-1);

error = zeros(size(phi_data));

for i = 1:size(error,1)
    for j = 1:size(error,2)
        p_x = (i-1)*h - c_x;
        p_y = (j-1)*h - c_y;
        
        r = sqrt(p_x^2 + p_y^2);
        if (isnan(phi_data(i,j)) == false)
            error(i,j) = ( A_sol*log(r) + B_sol );
        end
    end
end
error(error==0) = nan;

hold on;
%surf(phi_data);
surf(abs(error));
hold off


%% Time-evolution plot
data_array = [];

time_list = 0.001:0.001:9.899;
str_time_list = compose('%0.6f', time_list);

h5disp("export_data.h5");

for i = 1:length(str_time_list)
    dir = "/solution/time_data/" + str_time_list(i);
    mat = h5read("export_data.h5", dir);
    mat(mat==0) = nan;
    data_array(:,:,end+1) = mat;
end

myVideo = VideoWriter('myVideoFile', 'MPEG-4'); %open video file
myVideo.Quality = 100;
myVideo.FrameRate = 120;  %can adjust this, 5 - 10 works well for me
open(myVideo)

display_mat = data_array(:,:,1);
surface = surf(display_mat);
set(gcf,'position',[0, 0, 1024, 768])
zlim([0 250])

it = 1;
draw_it = 0;
while it < length(str_time_list)
    draw_it = draw_it + 1;
    if draw_it > 2
        surface.ZData = data_array(:,:,it);
        it = it + 1;
        draw_it = 0;
        frame = getframe(gcf);
        writeVideo(myVideo, frame);
    %else
    %    pause(0.000000001)
    end
    drawnow update
end
close(myVideo)
close all;
clear all;
clc;

%% Define the geometry of the spoke
r = 1.0;        % The initial apex height of the spoke, mm
a = 0.60;     % The thickness distribution coefficient
L = 75;       % The vertical height of the spoke, mm
t0 = 2.5;     % The initial honrizontal distance between the left boundary and right boundary, mm
N = 72;       % The number of spokes
%% Define the mesh space
nx = 3;                         % The numder of nodes in the x-direction
ny = 76;                        % The number of nodes in the y-direction
num_nodes = nx*ny;              % Total number of nodes
num_elems = (nx-1)*(ny-1);      % Total number of elements
plot_mesh = 'no';

%% Define the coordinate of the nodes
x0 = zeros(1,ny);
t = zeros(1,ny);
x_l = zeros(1,ny);
x_r = zeros(1,ny);
X = zeros(ny,nx);
y0 = linspace(0,L,ny);                                   % The partition of nodes in y-direction
for i=1:ny
    x0(i) = -r*(cos(2*pi*y0(i)/L)-1);                  % The backbone shape function
    t(i) = 2*a*t0*(cos(2*pi*y0(i)/L))^2+t0*(1-a);  % The thickness distribution fuction
    x_l(i) = x0(i)-t(i)/2;
    x_r(i) = x0(i)+t(i)/2;
    X(i,:) = linspace(x_l(i),x_r(i),nx);
end

%% Generate the node coordinate
nodes = zeros(num_nodes, 3);
z0 = zeros(ny,1);
for j=1:nx
    for k=1:ny
        nodes((j-1)*ny+k,:) = [X(k,j),y0(k),z0(k)];
    end
end

%% Translation
% Define translation vector
tx = 0;        % x-translation
ty = -361.9;   % y-translation
tz = 0;        % z-translation
translation_matrix = [tx; ty; tz];

for i=1:num_nodes
% Define initial coordinates of the node
xx0(i) = nodes(i,1); % initial x-coordinate
yy0(i) = nodes(i,2); % initial y-coordinate
zz0(i) = nodes(i,3); % initial z-coordinate

nodes_T(i,:) = [xx0(i); yy0(i); zz0(i)] + translation_matrix;
end

%% Rotation
% Deg_x = 0;  % The degree between each spoke (around x-axis)
% Deg_y = 0;  % The degree between each spoke (around y-axis)
Deg_z = 5;  % The degree between each spoke (around z-axis)
nodes_all = zeros(N*num_nodes,3);
for i=1:N
% Define rotation angles (in degrees)
alpha = 0;       % rotation around x-axis
beta = 0;        % rotation around y-axis
gamma = Deg_z*(i-1); % rotation around z-axis
% Convert rotation angles from degrees to radians
alpha = deg2rad(alpha);
beta = deg2rad(beta);
gamma = deg2rad(gamma);

% Define rotation matrices
Rx = [1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
Ry = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz = [cos(gamma) -sin(gamma) 0; sin(gamma) cos(gamma) 0; 0 0 1];
rotation_matrix = Rz * Ry * Rx;

nodes_R = zeros(num_nodes,3);
for j=1:num_nodes
    % Define initial coordinates of the node
Rx0(j) = nodes_T(j,1); % initial x-coordinate
Ry0(j) = nodes_T(j,2); % initial y-coordinate
Rz0(j) = nodes_T(j,3); % initial z-coordinate
% Apply rotation to the node coordinates
nodes_R(j,:) = rotation_matrix * [Rx0(j); Ry0(j); Rz0(j)];
end
nodes_all((i-1)*num_nodes+1:i*num_nodes,:)= nodes_R;
end

%% Generate the elements connected by nodes (4-nodes element)
elem_connect_all = zeros(N*num_elems, 4);
elem_connect = zeros(num_elems, 4);
for i = 1:N
  for j = 1:(nx-1)
    for k = 1:(ny-1)
        n1 = (j-1)*ny+k+(i-1)*num_nodes;          % Node 1
        n2 = n1 + ny;             % Node 2
        n3 = n2 + 1;              % Node 3
        n4 = n1 + 1;              % Node 4
        elem = (j-1)*(ny-1) + k;  % Element ID
        elem_connect(elem,:) = [n1 n2 n3 n4];
    end
  end
elem_connect_all((i-1)*num_elems+1:i*num_elems,:)=elem_connect;
end

%% export the nodes and elements informations
node_ids = (1:num_nodes*N)'; % Node IDs
elem_ids = (1:num_elems*N)'; % Element IDs
% node_name = sprintf('nodes-r%de-01-a%de-02.txt',r*10,a*100);
% elem_name = sprintf('elems-r%de-01-a%de-02.txt',r*10,a*100);
node_name = sprintf('nodes_insert.txt');
elem_name = sprintf('elems_insert.txt');
writematrix([node_ids, nodes_all], node_name); 
writematrix([elem_ids, elem_connect_all], elem_name); 

%% Plot the spoke mesh
% if strcmpi(plot_mesh,'yes')==1
%   for j=1:N*num_elems
%         n1 = elem_connect_all(j,1);
%         n2 = elem_connect_all(j,2);
%         n3 = elem_connect_all(j,3);
%         n4 = elem_connect_all(j,4);
%         mesh_x = [nodes_all(n1,1) nodes_all(n2,1) nodes_all(n3,1) nodes_all(n4,1) nodes_all(n1,1)];
%         mesh_y = [nodes_all(n1,2) nodes_all(n2,2) nodes_all(n3,2) nodes_all(n4,2) nodes_all(n1,2)];
%         plot(mesh_x,mesh_y,'k','LineWidth',0.5)
%         hold on;
%   end
%         axis equal;
% end

%% Write inp.files
% read the original .inp file
original_text = split(string(fileread('NPT-MAT-ABA-noSpokes.inp')), newline);

% read the node information and element information
insert1 = split(string(fileread(node_name)), newline);
insert2 = split(string(fileread(elem_name)), newline);
insert_text1 = insert1(1:N*num_nodes);
insert_text2 = insert2(1:N*num_elems);

% insert the node and element information into the designated row
insert_index = [4 8056];
insert_text = {insert_text1, insert_text2};
for ii = 1:length(insert_index)
    idx = insert_index(ii);
    original_text = [original_text(1:idx); insert_text{ii}; original_text(idx+1:end)];
    % update the row numbers after inserting new contents
    insert_index = insert_index + numel(insert_text{ii});
end

% output the modified .inp file
NPT_name = sprintf('NPT-vertical-%dunits-r%de-02-a%de-02-t%de-02.inp',N,r*100,round(a*100),t0*100);
fileID = fopen(NPT_name, 'w');
fprintf(fileID, '%s\n', original_text{:});
fclose(fileID);
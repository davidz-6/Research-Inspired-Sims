%Define the system
%Simple example f(x) = Ax
A = [1 0; 0 2];     
Gamma = eye(2);
K = A;                 %Assumption (9)
c = 2;

% Network topology:
edges = [1 2; 2 3; 3 4; 3 5];   %Each row represents one connection
N = max(edges(:));


%Build coupling matrix G
G = zeros(N);

for e = 1:size(edges,1)
    i = edges(e,1);
    j = edges(e,2);
    G(i,j) = 1;
    G(j,i) = 1;
end

k = sum(G,2);       
G = G - diag(k);      

%Calculate theta (Corollary 1)
theta = norm(K,2);

%Check C before pinning
D = zeros(N);
C = theta*eye(N) + c*(G-D);

fprintf('theta = %.2f\n',theta);
disp('G ='); disp(G);
disp('Degrees k ='); disp(k);
fprintf('Initial max eigenvalue(C) = %.4f\n',max(eig(C)));

%Select nodes to pin (Corollary 3)
%Unpinned nodes must satisfy k > theta/c
pinned = find(k <= theta/c);
unpinned = setdiff(1:N,pinned);

fprintf('theta/c = %.2f\n',theta/c);
disp('Pinned nodes ='); disp(pinned');

%Find gains (Theorem 3)
%Partition C into pinned and unpinned blocks
C0 = theta*eye(N) + c*G;

Ctilde = C0(unpinned,unpinned);

if ~isempty(unpinned) && max(eig(Ctilde)) >= 0
    error('Selected nodes are insufficient. Pin more nodes.');
end

if isempty(unpinned)
    gain_min = max(eig(C0));
else
    B = C0(pinned,unpinned);
    Ap = C0(pinned,pinned);
    gain_min = max(eig(Ap - B*(Ctilde\B')));
end

%Theorem 3 requires c*d > gain_min
d = max(0, gain_min/c) + 0.1;

D(pinned,pinned) = d*eye(length(pinned));

%Verify final C < 0 (Corollary 1)
C = theta*eye(N) + c*(G-D);

disp('D ='); disp(D);
disp('G-D ='); disp(G-D);
disp('Final C ='); disp(C);

fprintf('Selected gain d = %.4f\n',d);
fprintf('Max eigenvalue(C) = %.4f\n',max(eig(C)));

if max(eig(C)) < -1e-10
    disp('SUCCESS: C < 0, synchronization guaranteed.');
else
    disp('FAILED: C is not negative definite.');
end

%Plot network topology
network = graph(edges(:,1), edges(:,2));

figure('Color','w');
p = plot(network, 'Layout','force', ...
    'NodeLabel', 1:N, ...
    'MarkerSize', 14, ...
    'NodeColor', [0.3 0.6 0.9], ...
    'EdgeColor', [0.5 0.5 0.5], ...
    'LineWidth', 2, ...
    'NodeFontSize', 13);

% Highlight pinned nodes in red
highlight(p, pinned, ...
    'NodeColor', [0.9 0.2 0.2], ...
    'MarkerSize', 16);

title('Network Topology (Red = Pinned)', 'FontSize', 14);
axis off;
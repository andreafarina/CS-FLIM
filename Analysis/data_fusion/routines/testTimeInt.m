N = 16;
M = 256;

X = rand(N*M,1);
% create operator for integration
a = ones(1,N);
b = eye(M,M);
M = kron(a,b);
figure,imagesc(M)

%function I = time_integration_matrix(w, h, t)
w = 4; h = 4; t = 18;
i = kron(1 : w*h, ones(t, 1));
j = repmat(0 : w * h : w * h * t - 1, [1, w * h]) + kron(1 : w * h, ones(1, t));
I = sparse(i, j, ones(w*t*h, 1), w*h, w*h*t);
figure,imagesc(I)



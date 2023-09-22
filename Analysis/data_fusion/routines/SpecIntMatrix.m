function M = SpecIntMatrix(sizeVec)
%x(lambda,time,x,y)
Nl = sizeVec(1);
Nt = sizeVec(2);
Nx = sizeVec(3);
Ny = sizeVec(4);

%% direct operator
M = spalloc(Nt*Nx*Ny,Nl*Nt*Nx*Ny,Nt*Nx*Ny);
vint = sparse(ones(1,Nl));
% vrep = repmat(vint,[1,Nt]);
% row = spalloc(1,Nt*Nl,Nl);
% row(1:Nt*Nl) = vrep;
% % single-pixel integrator
% for i = 1:Nl
%     Msp(i,:) = row;
%     row = circshift(row,1);
% end
M = kron(speye(Nt*Nx*Ny,Nt*Nx*Ny),vint);
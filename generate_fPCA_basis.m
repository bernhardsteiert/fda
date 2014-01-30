% Run paper_figure_10_fPCA_doseeffect.m first

fdobj = getcoef(c_signal_pcastr.harmfd);
fdobj = repmat(flipharm,size(fdobj,1),1).*fdobj;
basisobj = getbasis(c_signal_pcastr.harmfd);

newcoef = fdobj;
for irow = 1:size(Rmat,1)
    newcoef(:,irow) = sum(repmat(Rmat(irow,:),size(fdobj,1),1) .* fdobj(:,1:size(Rmat,1)),2);
end

harm_fPCA = fd(newcoef,basisobj,{'time','Rotated harmonics','Rotated harmonics for variables'});
% plot(harm_fPCA)
% harm_basis_fPCA = create_fd_basis(harm_fPCA);
% plot(harm_basis_fPCA)
harm_basis_full = create_fd_basis(harm_fPCA);
plot(harm_basis_full)

% save('harm_basis_fPCA','harm_basis_fPCA')
save('harm_basis_full','harm_basis_full')
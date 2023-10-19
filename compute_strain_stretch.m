function [Strain_cells, RightCauchy_cells, LeftCauchy_cells, Stretch_cells, Eigenvector_cells] = compute_strain_stretch(F_cells, ...
    compute_stretch)
if isempty(compute_stretch)
    compute_stretch = false;
end

RightCauchy_cells = pagemtimes(F_cells, 'transpose', F_cells, 'none');
LeftCauchy_cells = pagemtimes(F_cells, 'none', F_cells, 'transpose');

I = zeros(size(F_cells));
I(1,1,:,:) = 1;
I(2,2,:,:) = 1;
I(3,3,:,:) = 1;
Strain_cells = 0.5 * (RightCauchy_cells - I);

Stretch_cells = zeros(3, size(F_cells, 3), size(F_cells, 4));
Eigenvector_cells = zeros(size(F_cells));

if compute_stretch
    for ii = 1:size(F_cells, 3)
        for jj = 1:size(F_cells, 4)
            [V,D] = eig(RightCauchy_cells(:,:,ii,jj));
            D = abs(D); %D should be positive (check this). Making sure for sqrt operation. Any negative values will be very close to 0.
            Stretch_cells(:,ii,jj) = sqrt(diag(D));
            Eigenvector_cells(:,:,ii,jj) = V;
        end
    end
end

end


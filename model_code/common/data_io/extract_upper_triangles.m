function triu_vect = extract_upper_triangles(cmat)

% Extract upper triangles and stack vectorized matrices to format connectomic data for modeling

% Joseph Griffis 2024

triu_inds = find(triu(ones(size(cmat,[2,3])), 1));
triu_vect = zeros(size(cmat, 1), length(triu_inds));
for i = 1:size(cmat,1)
    temp = squeeze(cmat(i, :,:));
    triu_vect(i,:) = temp(triu_inds);
end

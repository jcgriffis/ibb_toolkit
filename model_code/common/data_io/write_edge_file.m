function write_edge_file(t, tu_mask, X_ind, stat, direction, out_name)

% Format matrix and write .edge file
sdc_inds = find(tu_mask); % Get upper triangle indices of connectivity matrix
tu = zeros(size(tu_mask)); % empty container matrix
tu(sdc_inds(X_ind)) = stat; % fill container matrix with 
tl = rot90(flipud(tu),3);
wmat = tu + tl;
dlmwrite([out_name '.edge'], wmat, 'delimiter', '\t', 'precision', 4); %Write out .edge file

if ~isempty(t) % Write out .node file if parcellation table is present
   
    % Format node parameters    
    if ~istable(t)
        load(t);
    end
    node_label = cellstr(t.RegionName); %t.RegionName; % n_regions by 1 cell array of strings corresponding to node labels (i.e. parcel names)
    node_color = t.NetworkID; %t.NetworkID; % n_regions-by-1 array of integer values corresponding to e.g. network assignments or partitions (used to color nodes in external viewers)
    node_pos = [t.X, t.Y, t.Z]; % Node coordinates
    
    % Write out .node files
    if strcmp(direction, 'pos')
        wmat(wmat < 0) = 0;
        node_size = sum(wmat,2)./max(sum(wmat,2)); % size nodes according to sum of connection weights
    elseif strcmp(direction, 'neg')
         wmat = wmat .* -1;
         wmat(wmat < 0) = 0;
         node_size = sum(wmat,2)./max(sum(wmat,2)); % size nodes according to sum of connection weights 
    else
        node_size = ones(size(wmat,1),1);
    end
    
    % construct output filename and save
    gen_node_file(node_pos, node_color, node_size, node_label, [out_name '.node']);
end

end
function gen_node_file(node_pos, node_color, node_size, node_label, output_file)

% Function is copy of GRETNA's gretna_gen_node_file.m

C=[node_pos, node_color, node_size];
C=num2cell(C);
C=[C, node_label];

fid=fopen(output_file, 'w');
Tmp=C';
fprintf(fid, '%g\t%g\t%g\t%d\t%f\t%s\n', Tmp{:});
fclose(fid);
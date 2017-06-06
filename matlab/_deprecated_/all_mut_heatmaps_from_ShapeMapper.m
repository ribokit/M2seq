function all_mut_heatmap_from_ShapeMapper( SHAPEmapper_config_file )
% all_mut_heatmap_from_ShapeMapper( SHAPEmapper_config_file )
%
% INPUTS:
%  SHAPEmapper_config_file = SHAPEmapper .cfg file in same directory as
%                             output directory output_A/counted_mutations/
%
% (C) Rhiju Das, Stanford University, 2017

set(figure(1),'position',[158,119,948,836] );

[primer_tag,RNA_tag] = get_alignments( SHAPEmapper_config_file );
dirname = fileparts( SHAPEmapper_config_file );
outdir = [dirname, '/output_A/counted_mutations/' ];

% plot 8 at a time
numrows = 8;clf;
epsdir = [dirname,'/mut_heatmaps/'];
if ~exist( epsdir ); mkdir( epsdir ); end;
plot_count = 0;
for i = 1:length( RNA_tag )
    csv_name = [outdir,primer_tag{i},'_',RNA_tag{i},'.csv'];
    fasta_file = [dirname,'/',RNA_tag{i},'.fa'];
    sequence = fastaread(fasta_file);
    sequence = sequence.Sequence;
    mut_heatmap( csv_name,1:length(sequence),sequence,0,[],[numrows,mod(i-1,numrows)+1]);

    if mod( i, numrows ) == 0 | i == length( RNA_tag ); 
        plot_count = plot_count+1;
        eps_file = [epsdir,'mut_heatmap',num2str(plot_count),'.eps'];
        fprintf( 'Creating: %s\n',eps_file );
        print(eps_file,'-depsc2');
    end
    if mod( i, numrows ) == 0 & i ~= length( RNA_tag ); 
        %pause;
        clf; 
    end;
end


function [primer_tag,RNA_tag] = get_alignments( SHAPEmapper_config_file );
% get [alignments] information
fid = fopen( SHAPEmapper_config_file, 'r' );
READ_STUFF = 0;
primer_tag = {};
RNA_tag = {};
count = 0;
while 1
    line = fgetl( fid );
    if ~ischar( line ); break; end;
    if ~isempty( strfind( line, '[alignments]' ) );  
        READ_STUFF = 1;
    elseif ~isempty( strfind( line, '[profiles]' ) );
        READ_STUFF = 0;
    end
    if READ_STUFF 
        cols = strsplit( line, ':' );
        if length( cols ) == 2 && line(1) ~= '#'
            count = count+1;
            primer_tag{count} = cols{1};
            cols = strsplit( cols{2}, '=' );
            RNA_tag{count} = strrep(cols{2},' ','');
        end
    end
end
fclose( fid );


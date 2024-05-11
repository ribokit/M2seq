function run_m2seq_plots( r_dms, r_nomod, name, output_folder, iterations, offset, shapeknots, fiveprime_trim, threeprime_trim, new_struct,DMS);
% run_m2seq_plots( r_dms, r_nomod, name, output_folder, iterations, offset, shapeknots, fiveprime_trim, threeprime_trim, new_struct,DMS);
%
% Inputs
%  r_dms   = RDAT file with DMS data, output from M2-seq
%  r_nomod = RDAT file with no modification data, output from M2-seq
%  name    = tag to use in filename [default: 'test']
%  output_folder = name of output folder, will be created if it doesn't
%                  exist [default: 'm2seq_plots']
%  iterations = number of bootstraps for statistical confidence estimation
%                of stems [Default: 10, for publication increase to 100-1000]
%  offset   = offset to apply to sequence numbering for secondary structure
%              figures
%  shapeknots = run secondary structure  modeling with pseudoknots (slow!)
%                or not [Default: 0 ]
%  fiveprime_trim = starting position of sequence [Default: 1]
%  threeprime_trim = nucleotides to remove from 3' end of sequence [Default: 0]
%  new_struct = structure in dotbracket notation to make VARNA plots with
%                 [Default: '', use M2-seq inferred sructure]
%  DMS      = is a DMS-based M2-seq experiment [Default: 1]
%
% (C) Das laboratory, Stanford University, 2018-2024
% Adapted from original scripts by Ved Topkar
%

if ~exist('name','var') name = 'test'; end;
if ~exist('output_folder','var') output_folder = 'm2seq_plots'; end;
if ~exist('iterations','var') | isempty(iterations), iterations = 10; end;
if ~exist('offset','var') | isempty(offset), offset = 0; end;
if ~exist('shapeknots','var') shapeknots = 0; end;
if ~exist('fiveprime_trim','var') fiveprime_trim = 1; end;
if ~exist('threeprime_trim','var') threeprime_trim = 0; end;
if ~exist('new_struct','var') new_struct = ''; end;
if ~exist('DMS','var') DMS = 0; end;

if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

% Load in the rdat files
rna_dms = read_rdat_file(r_dms);
rna_nomod = read_rdat_file(r_nomod);

% Extract the reference sequence from the DMS rdat
sequence = rna_dms.sequence;
sequence_length = length(sequence);
trimmed_sequence = sequence(fiveprime_trim:sequence_length-threeprime_trim);

% Graph the raw data
graph_raw(rna_dms, name, 'DMS', output_folder);
graph_raw(rna_nomod, name, 'NM', output_folder);

% Compute the Zscore matrix
Z = output_Zscore_from_rdat( [], r_dms, r_nomod, [], 1, 1 );
% Subset Zscores to trim ends
Z = Z(fiveprime_trim:sequence_length-threeprime_trim, fiveprime_trim:sequence_length-threeprime_trim)
%show_2d_das()
show_2d_das(-Z, name, offset, output_folder);

if DMS
    % Get 1D dms profile (again, from dms and nomod rdats)
    one_dimensional_dms = get_DMS_profile(r_dms, r_nomod)
    one_dimentional_dms_normalized = DMS_normalize(one_dimensional_dms, trimmed_sequence)

    % Save 1D data to CSV
    csvwrite(fullfile(output_folder,'1D_DMS_normalized.csv'), one_dimentional_dms_normalized)
else
    % Get 1D dms profile (again, from dms and nomod rdats)
    one_dimensional_dms = rna_dms.reactivity(:,1) - rna_nomod.reactivity(:,1)
    one_dimentional_dms_normalized = SHAPE_normalize(one_dimensional_dms)

    % Save 1D data to CSV
    csvwrite(fullfile(output_folder,'1D_SHAPE_normalized.csv'), one_dimentional_dms_normalized)
end

% Predict helices with M2Net
clf
close
m2net_figure = figure
[Z_m2net, structure_m2net] = m2net( Z,  trimmed_sequence, '', 0);
saveas(gcf, fullfile(output_folder, [name, '_', 'M2net_plot', '.png']))
clf
close
show_2dmap(Z_m2net, '', 0);

% Run bootstrapped RNA_structure with X iterations
if DMS
    bootstrapped_figure = figure
    [final_structure, bpp, one_dimensional_filtered] = rna_structure(trimmed_sequence, '', '', '', Z, iterations, shapeknots, one_dimentional_dms_normalized);
    clf
    close
else
    bootstrapped_figure = figure
    [final_structure, bpp, one_dimensional_filtered] = rna_structure(trimmed_sequence, one_dimentional_dms_normalized, '', '', Z, iterations, shapeknots, '');
    clf
    close
end

% Generate figure for bootstraped pairing probabilities
bpp_figure = figure
show_2dmap(bpp, '', offset, 1, 0);
pbaspect([1 1 1]); % Square aspect ratio
set(gca, 'FontSize', 14);
xlabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
colorbar;
title('Bootstrapped Base Pairing Probabilities');
saveas(gcf, fullfile(output_folder, [name, '_', 'bpp', '.svg']));
savefig(fullfile(output_folder, [name, '_bpp','.fig']));
clf
close

% Generate VARNA fig
% We pass in bpp to display the bootstrapping-derived helix probabilities
figure
output_varna(fullfile(output_folder, [name, '_VARNA' ,'.svg']), trimmed_sequence, final_structure, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
output_varna(fullfile(output_folder, [name, '_VARNA']), trimmed_sequence, final_structure, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')

if ~strcmp(new_struct, '');
    figure
    output_varna(fullfile(output_folder, [name, '_VARNA' ,'.svg']), trimmed_sequence, new_struct, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
    output_varna(fullfile(output_folder, [name, '_VARNA']), trimmed_sequence, new_struct, '', '', 0, '', [], one_dimentional_dms_normalized, bpp, '')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_2d_das(bpp, name, offset, output_folder)
show_2dmap(bpp, '', offset, 1, 0);
pbaspect([1 1 1]); % Square aspect ratio
set(gca, 'FontSize', 14);
xlabel('Sequence Position', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Mutation Position', 'FontSize', 20, 'FontWeight', 'bold');
saveas(gcf, fullfile(output_folder, [name, '_Z.svg']))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function graph_raw(data, name, treatment, output_folder);

sequence = data.sequence;
sequence_length = length(sequence);

fig = figure

show_rdat(data);
pbaspect([1 1 1]); % Square aspect ratio
title([name, ' raw ', treatment])
set(gca, 'FontSize', 20)
xlabel('Sequence Position', 'FontSize', 30, 'FontWeight', 'bold')
ylabel('Mutation Position', 'FontSize', 30, 'FontWeight', 'bold')

set(gca, 'XTickLabel', (50:50:sequence_length))
set(gca, 'XTick', (50:50:sequence_length))
set(gca, 'YTickLabel', (50:50:sequence_length))
set(gca, 'YTick', (50:50:sequence_length))

saveas(gcf, fullfile(output_folder, [name, '_', treatment ,'_raw','.pdf']))
savefig(fullfile(output_folder, [name, '_', treatment ,'_raw','.fig']))
clf
close


function [mut_avg, mut_max, muts, muts_diff] = mut_heatmap(csv_path, seqpos, sequence, offset, titl, save_path)
% [mut_avg, mut_max, muts, muts_diff] = mut_heatmap(csv_path, seqpos, sequence, offset, titl, save_path)
%
% Generates both mutation spectra across full sequence and average/maximum mutation matrix
%
%  INPUTS 
%   csv_path = csv format output(s) from SHAPEmapper "counted_mutations" directory:
%                -- supply directory name to get mutation spectra for all
%                     .csv files in directory (default: './' )
%                -- supply a single file to get mutation spectrum
%                -- supply cell with more than one file to do the DIFF between the two.
%   seqpos   = sequence positions to plot (can use conventional numbering -- 
%               length of this vector determines how much to plot)
%   sequence = needs to be full sequence given as ShapeMapper input
%   offset   = number of residues into sequence that region of interest begins minus 1 (so if ROI begins at 14, enter 13 for offset)
%                Note: this is not the same as offset typically used in
%                RDATs, which shifts to conventional numbering
%   title    = title to put on figure [default: csv name]
%   save_path= where to save [default: no save]
%               if supplied as a *string*, save each figure to file with that name.
%               if supplied as a doublet of integers, like [9,1],
%                 show the plots on row 1 out of 9 rows in the figure.
%
% OUTPUTS
%   mut_avg  = average mutation matrix  [4x5]
%   mut_max  = max mutation  matrix [4x5]
%   muts     = mutation specturm per position [17 X N], as follows
%              'del A','del T','del G', 'del C', ...
%              'A->T','A->G','A->C', ...
%              'T->A','T->G','T->C', ...
%              'G->A','G->T','G->C', ...
%              'C->A','C->T','C->G', ...
%              'Total counts'
%
% (C) Clarence Cheng, Das laboratory, Stanford University 2016-2017
% (C) Rhiju Das, 2017
mut_avg = [];
mut_max = [];
muts = [];
muts_diff = [];

if ~exist( 'csv_path' ) | isempty( csv_path )
    csv_path = './';
end
if ischar(csv_path);
    DIFF = 0;
    if exist( csv_path, 'dir' )
        run_on_all_csvs( csv_path );
        return;
    end
elseif iscell(csv_path);
    DIFF = 1;
end

if ~exist('titl','var') || isempty(titl);
    if ischar( csv_path )
        tmp  = strsplit(csv_path,'/');
        titl = tmp{end};
    else
        for i = 1:length( csv_path )
            tmp  = strsplit(csv_path{i},'/');
            titl{i} = tmp{end};
        end
    end
end
if ~exist('offset','var') || isempty( offset );
    offset = 0;
end

SAVE_EPS   = exist( 'save_path', 'var' ) && ~isempty( save_path ) && ischar( save_path);
DO_SUBPLOT = exist( 'save_path', 'var' ) && isnumeric( save_path );
if DO_SUBPLOT
    assert( length( save_path ) == 2 );
    num_rows = save_path(1); row = save_path(2);
    assert( row <= num_rows );
end

mutation_list = {
'del A', ...
'del T', ...
'del G', ...
'del C', ...
'A->T', ...
'A->G', ...
'A->C', ...
'T->A', ...
'T->G', ...
'T->C', ...
'G->A', ...
'G->T', ...
'G->C', ...
'C->A', ...
'C->T', ...
'C->G', ...
};
set(gcf, 'PaperPositionMode','auto','color','white');

if DIFF == 0;
    muts_init = csvread(csv_path,2,3);
    for i = 1:17;
        for j = 1:size(muts_init,2)-1;
            muts(i,j) = muts_init(i,j)/muts_init(18,j);
        end
    end
    if ~exist('seqpos','var') || isempty(seqpos);
        seqpos = 1:size(muts,2)-1;
    end
    colmap = jet(99);
    startpos = offset+1;
    seqrange = [ startpos : (startpos+length(seqpos)-1) ];
    % figure; hist(reshape(muts(1:16,1:end-1),size(muts(1:16,1:end-1),1)*size(muts(1:16,1:end-1),2),1), 1000)
    clims = [0 0.004];                                 % specify constant limits for image color data mapping, so images of different datasets are comparable
    if DO_SUBPLOT; axes('Position', [0.05,1-(row/num_rows), 0.7, 1/num_rows] ); else; figure; end;
    imagesc(seqpos,1:16,muts(1:16,seqrange),clims); axis image; 
    set(gca,'tickdir','out','ytick',1:16,'yticklabel',mutation_list,'fontsize',10,'ticklength',[0.0025 0.025]); colormap(gca,colmap); hold on;
    xlabel('Sequence position'); ylabel('Mutation'); title([titl ': Mutation spectrum'],'fontsize',10,'interpreter','none');
    maxpos  = seqpos(max(find( mod(seqpos,10) == 0 )));
    minpos  = seqpos(min(find( mod(seqpos,10) == 0 )));
    make_lines(minpos-1:10:maxpos,'k',0.5); make_lines_horizontal(4:3:16,'k',0.5);
    % set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    colorbar('ytick',[0 0.004]);
    if SAVE_EPS; save_fig(save_path); end

    % would be good to just get sequence from file...
    if ~exist('sequence','var') || isempty( sequence );
        fid = fopen(csv_path,'r');
        seqcell = textscan(fid,'%s','delimiter','\n');
        seqarry = strsplit(seqcell{1}{2},',');
        sequence = strjoin(seqarry(4:end),'');
    end
    
    ind = find_idx(sequence, offset, seqpos);
    [mut_avg, mut_max] = mut_matrix(muts, ind);

    clims = [0 0.025];
    if DO_SUBPLOT; axes('Position', [0.77,1-(row/num_rows), 0.1, 1/num_rows] ); else; figure; end;
    imagesc(mut_avg,clims); axis image; set(gca,'tickdir','out','ytick',1:4,'xtick',1:5,'yticklabel',{'A','U','G','C'},'xticklabel',{'A','U','G','C','del'},'xaxisLocation','top','fontsize',30,'ticklength',[0.0025 0.025]); 
    if DO_SUBPLOT; set(gca,'fontsize',8); end;
    colormap(gca,1-gray(100)); hold on;
    if DO_SUBPLOT; title( 'Avg. mut. rates' ); else; set(gca,'LooseInset',get(gca,'TightInset')); title([titl ': Average mutation rates'],'fontsize',20,'interpreter','none');end;
    colorbar('ytick',[0 0.025]);
    if SAVE_EPS; save_fig([save_path '_muts_avg']); end
    
    if DO_SUBPLOT; axes('Position', [0.87,1-(row/num_rows), 0.1, 1/num_rows] ); else; figure; end;
    imagesc(mut_max,clims); axis image; set(gca,'tickdir','out','ytick',1:4,'xtick',1:5,'yticklabel',{'A','U','G','C'},'xticklabel',{'A','U','G','C','del'},'xaxisLocation','top','fontsize',30,'ticklength',[0.0025 0.025]); 
    if DO_SUBPLOT; set(gca,'fontsize',8); end;
    colormap(gca,1-gray(100)); hold on;
    set(gca,'LooseInset',get(gca,'TightInset')); 
    if DO_SUBPLOT; title('Max mut. rates');else; title([titl ': Maximum mutation rates'],'fontsize',20,'interpreter','none');end;
    colorbar('ytick',[0 0.025]);
    if SAVE_EPS; save_fig([save_path '_muts_max']); end

elseif DIFF == 1;
    for k = 1:length(csv_path);
        muts_init{k} = csvread(csv_path{k},2,3);
        muts{k} = muts_init{k};
        for i = 1:16;
            for j = 1:size(muts{k},2)-1;
                muts{k}(i,j) = muts_init{k}(i,j)/muts_init{k}(18,j);
            end
        end
        if ~exist('seqpos','var') || isempty(seqpos);
            seqpos = 1:size(muts{k},2)-1;
        end
    end
    cmap = redblue(99);
    for k = 1:length(csv_path)-1;
        muts_diff = (muts{1}(1:16,:)-muts{k+1}(1:16,:));
        clims = [-0.004 0.004];                                 % specify constant limits for image color data mapping, so images of different datasets are comparable
        startpos = offset+1;
        seqrange = startpos:startpos+length(seqpos)-1;
        % figure; hist(reshape(muts_diff(1:16,1:end-1),size(muts_diff(1:16,1:end-1),1)*size(muts_diff(1:16,1:end-1),2),1), 1000)
        figure(4); imagesc(seqpos,1:16,muts_diff(:,1:end-1),clims); axis image; 
        set(gca,'tickdir','out','ytick',1:16,'yticklabel',mutation_list,'fontsize',10,'ticklength',[0.0025 0.025]); colormap(cmap); hold on;
        xlabel('Sequence position'); ylabel('Mutation'); title(titl{k},'fontsize',10,'interpreter','none');
        maxpos  = seqpos(max(find( mod(seqpos,20) == 0 )));
        minpos  = seqpos(min(find( mod(seqpos,20) == 0 )));
        make_lines(minpos-1:10:maxpos,'k',0.5); make_lines_horizontal(4:3:16,'k',1);
        make_lines_horizontal(1:16,[0.6 0.6 0.6],0.25);
        if SAVE_EPS; save_fig(save_path{k}); end
    end
end

function cmap = redblue( N );
cmap = ones( N, 3 );
n = (N-1)/2;
cmap( 1:n, 1) = [1:n]/n;; % blue
cmap( 1:n, 2) = [1:n]/n;; % blue
cmap( (n+1) + [1:n], 2) = [n:-1:1]/n; % red
cmap( (n+1) + [1:n], 3) = [n:-1:1]/n; % red

function [ind] = find_idx(sequence, offset, seqpos)
seq = sequence;
ind = {};
ind{1} = findstr(seq, 'A');
ind{2} = findstr(seq, 'T');
ind{3} = findstr(seq, 'G');
ind{4} = findstr(seq, 'C');
for i = 1:length(ind);
    for j = length(ind{i}):-1:1;
        if ind{i}(j) <= offset || ind{i}(j) > offset+length(seqpos);
            ind{i}(j) = [];
        end
    end
end


function [mut_avg, mut_max] = mut_matrix(muts, ind)
% take average across mutation type (along row) over all native positions for each starting nucleotide (A, U, G, C)
mut_avg = zeros(4,5);
mut_max = zeros(4,5);
for i = 1:17;
    if i >= 1 && i <= 4;
        row = i;
        col = 5;
        mut_avg(row, col) = mean(muts(i, ind{i}));
        mut_max(row, col) = max( muts(i, ind{i}));
    elseif i >= 5 && i <= 7;
        row = 1;
        col = i - 3;
        mut_avg(row, col) = mean(muts(i, ind{1}));
        mut_max(row, col) = max( muts(i, ind{1}));
    elseif i >= 8 && i <= 10;
        row = 2;
        if i == 8; col = i - 7; else; col = i - 6; end
        mut_avg(row, col) = mean(muts(i, ind{2}));
        mut_max(row, col) = max( muts(i, ind{2}));
    elseif i >= 11 && i <= 13;
        row = 3;
        if i == 13; col = i - 9; else; col = i - 10; end
        mut_avg(row, col) = mean(muts(i, ind{3}));
        mut_max(row, col) = max( muts(i, ind{3}));
    elseif i >= 14 && i <= 16;
        row = 4;
        col = i - 13;
        mut_avg(row, col) = mean(muts(i, ind{4}));
        mut_max(row, col) = max( muts(i, ind{4}));
    elseif i == 17;
        for j = 1:4;
            mut_avg(j,j) = mean( 1 - muts(i, ind{j}) );
            mut_max(j,j) = max(  1 - muts(i, ind{j}) );
        end
    end
end


function save_fig(save_path)
print( gcf, '-depsc2', '-loose', '-r300', save_path);
fprintf( ['Created: ', save_path, '\n'] );
hgsave( save_path );


%%%%%%%%%%
function run_on_all_csvs( csv_path );

csv_names = dir( [csv_path, '/*.csv' ] );

% plot 8 at a time
numrows = 8;
clf;
epsdir = [csv_path,'/mut_heatmaps/'];
if ~exist( epsdir ); mkdir( epsdir ); end;
plot_count = 0;
for i = 1:length( csv_names )
    csv_name = csv_names(i).name;
    sequence = [];
    mut_heatmap( csv_name,[],[],[],[],[numrows,mod(i-1,numrows)+1]);
    if mod( i, numrows ) == 0 | i == length( csv_names ); 
        plot_count = plot_count+1;
        eps_file = [epsdir,'mut_heatmap',num2str(plot_count),'.eps'];
        fprintf( 'Creating: %s\n',eps_file );
        print(eps_file,'-depsc2');
    end
    if mod( i, numrows ) == 0 & i ~= length( csv_names ); 
        pause;
        clf; 
    end;
end

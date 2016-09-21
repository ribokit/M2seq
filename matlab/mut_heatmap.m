function [mut_avg, mut_max, muts, muts_diff] = mut_conv_heatmap(csv_path, seqpos, sequence, offset, titl, save_path)
% Generates both mutation spectra across full sequence and average/maximum mutation matrix
%   sequence needs to be full sequence given as ShapeMapper input
%   offset is number of residues into sequence that region of interest begins minus 1 (so if ROI begins at 14, enter 13 for offset)

seq = sequence
ind = {};
ind{1} = findstr(seq, 'A')
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
ind{1}
ind{2}
ind{3}
ind{4}

if ischar(csv_path);
    DIFF = 0;
elseif iscell(csv_path);
    DIFF = 1;
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

if DIFF == 0;
    muts_init = csvread(csv_path,2,3);
    muts = muts_init;
    for i = 1:16;
        for j = 1:size(muts,2)-1;
            muts(i,j) = muts_init(i,j)/muts_init(18,j);
        end;
    end;
    if ~exist('seqpos','var') || isempty(seqpos);
        seqpos = 1:size(muts,2)-1;
    end;
    colmap = jet(99);
    % figure; hist(reshape(muts(1:16,1:end-1),size(muts(1:16,1:end-1),1)*size(muts(1:16,1:end-1),2),1), 1000)
    clims = [0 0.004];                                 % specify constant limits for image color data mapping, so images of different datasets are comparable
    figure; imagesc(seqpos,1:16,muts(1:16,1:end-1),clims); axis image; set(gca,'tickdir','out','ytick',1:16,'yticklabel',mutation_list,'fontsize',10,'ticklength',[0.0025 0.025]); colormap(colmap); hold on; % image color data mapping is direct, so images of different datasets are comparable
    xlabel('Sequence position'); ylabel('Mutation'); title(titl,'fontsize',10);
    maxpos  = seqpos(max(find( mod(seqpos,20) == 0 )));
    minpos  = seqpos(min(find( mod(seqpos,20) == 0 )));
    make_lines(minpos:10:maxpos,'k',0.5); make_lines_horizontal(4:3:16,'k',0.5);
    % set(gca, 'Position', get(gca, 'OuterPosition') - get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    if exist( 'save_path', 'var' ) && ~isempty( save_path ); save_fig(save_path); end;

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

    clims = [0 0.025];
    figure; imagesc(mut_avg,clims); axis image; set(gca,'tickdir','out','ytick',1:4,'xtick',1:5,'yticklabel',{'A','U','G','C'},'xticklabel',{'A','U','G','C','del'},'xaxisLocation','top','fontsize',30,'ticklength',[0.0025 0.025]); colormap(1-gray(100)); hold on;
    set(gca,'LooseInset',get(gca,'TightInset')); title(titl,'fontsize',20);
    if exist( 'save_path', 'var' ) && ~isempty( save_path ); save_fig([save_path '_muts_avg']); end;
    
    figure; imagesc(mut_max,clims); axis image; set(gca,'tickdir','out','ytick',1:4,'xtick',1:5,'yticklabel',{'A','U','G','C'},'xticklabel',{'A','U','G','C','del'},'xaxisLocation','top','fontsize',30,'ticklength',[0.0025 0.025]); colormap(1-gray(100)); hold on;
    set(gca,'LooseInset',get(gca,'TightInset')); title(titl,'fontsize',20);
    if exist( 'save_path', 'var' ) && ~isempty( save_path ); save_fig([save_path '_muts_max']); end;

elseif DIFF == 1;
    for k = 1:length(csv_path);
        muts_init{k} = csvread(csv_path{k},2,3);
        muts{k} = muts_init{k};
        for i = 1:16;
            for j = 1:size(muts{k},2)-1;
                muts{k}(i,j) = muts_init{k}(i,j)/muts_init{k}(18,j);
            end;
        end;
        if ~exist('seqpos','var') || isempty(seqpos);
            seqpos = 1:size(muts{k},2)-1
        end;
    end;
    cmap = redblue(99);
    for k = 1:length(csv_path)-1;
        muts_diff = (muts{1}(1:16,:)-muts{k+1}(1:16,:));
        clims = [-0.004 0.004];                                 % specify constant limits for image color data mapping, so images of different datasets are comparable
        figure; hist(reshape(muts_diff(1:16,1:end-1),size(muts_diff(1:16,1:end-1),1)*size(muts_diff(1:16,1:end-1),2),1), 1000)
        figure; imagesc(seqpos,1:16,muts_diff(:,1:end-1),clims); axis image; set(gca,'tickdir','out','ytick',1:16,'yticklabel',mutation_list,'fontsize',10,'ticklength',[0.0025 0.025]); colormap(cmap); hold on;
        xlabel('Sequence position'); ylabel('Mutation'); title(titl{k},'fontsize',10);
        maxpos  = seqpos(max(find( mod(seqpos,20) == 0 )));
        minpos  = seqpos(min(find( mod(seqpos,20) == 0 )));
        make_lines(minpos:10:maxpos,'k',0.5); make_lines_horizontal(4:3:16,'k',0.5);
        if exist( 'save_path', 'var' ) && ~isempty( save_path ); save_fig(save_path{k}); end;
    end;
end

function save_fig(save_path)
print( gcf, '-depsc2', '-loose', '-r300', save_path);
fprintf( ['Created: ', save_path, '\n'] );
hgsave( save_path );

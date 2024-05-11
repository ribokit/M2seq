function scalefactor = m2seqplot( R, save_path, scalefactor, titl, secstr, fntsize, pdb, c2 )
% m2seqplot( R, save_path, scalefactor, titl, secstr, fntsize, pdb, c2 )
%
% For plotting 2D mutational profiling data (M2-seq):
% Example command: map2dplot(rdat_var, '/path/to/filename', '', 'WK 10/2015 AMPure +DMS, <__ muts, ___ reads');
%
% REQUIRED input
% R           = RDAT file or RDAT object  with M2-seq data.
%
% OPTIONAL inputs
% save_path   = .eps file in which to save data (default: [], no output)
% scalefactor = scaling of data (default: autoscaling based on 10/median)
% title       = plot title (default: R.name )
% secstr      = secondary structure to show as squares. Give 1 to use RDAT.structure, or specify in dot-bracket notation (default 0). 
% fntsize     = font size in points (default: 15 )
% pdb         = PDB file from which to estimate 15 and 30- Angstrom distance contours.
% c2          = coloring for 3D contours (default: [.5 0 .5; 0.3 0.5 1],
%                  i.e., purple, marine ] ); 
%
% For pdb contours, you'll need to install MAPseeker (available on github).
%
% (C) C. Y. Cheng and Das lab, 2015-2017

if ~exist( 'fntsize', 'var' ) | isempty(fntsize); fntsize = 15; end
if ~exist( 'c2','var' ) | isempty( c2 ); c2 = [.5 0 .5; 0.3 0.5 1]; end;

if ischar( R );
    r = read_rdat_file( R );
else
    r = R;
end

reactivity = r.reactivity;
seqpos = r.seqpos;
if ~exist( 'scalefactor' , 'var' ) | isempty(scalefactor); scalefactor = 10/median(median(nonzeros(r.reactivity))); end
xrng = [seqpos(1)-1 seqpos];

% scalefactor

set(gcf, 'PaperPositionMode', 'Manual','PaperOrientation', 'Portrait','PaperPosition', [-0.65 0.15 12 8]);%,'color','white','Position', [0 0 800 600]);
image(seqpos,xrng,reactivity'*scalefactor); colormap(1-gray(100)); axis image;
% image(1:size(reactivity,2),1:size(reactivity,1),reactivity*scalefactor); colormap(1-gray(100)); axis image;
% image(xrng,seqpos,reactivity*scalefactor); colormap(1-gray(100)); axis image;

% Label x and y axes
axlabel = {'WT'};
maxpos  = seqpos(max(find( mod(seqpos,20) == 0 )));
minpos  = seqpos(min(find( mod(seqpos,20) == 0 )));
tickpos = fliplr(maxpos:-20:minpos);
count = 1;
for i = tickpos;
    count = count + 1;
    % the sequence characters are not useful and crowd the axis labels -- rhiju:
    %axlabel{count} = [r.sequence(i-r.offset), num2str(i)]; 
    axlabel{count} = [num2str(i)];
end;
set(gca,'xtick',tickpos, 'fontsize',fntsize, 'XTickLabel', axlabel(2:end),'xticklabelrotation',90 )
set(gca,'ytick',[r.offset tickpos], 'fontsize',fntsize, 'YTickLabel', axlabel )
set(gca,'TickDir','out');
xlabel( 'Sequence Position','fontsize',fntsize+5,'fontweight','bold' );
ylabel( 'Mutation Position','fontsize',fntsize+5,'fontweight','bold' );
hold on;

% Add title
if ~exist( 'titl', 'var' ) | isempty( titl ); titl = r.name; end;
title( titl, 'fontsize', fntsize+5, 'fontw', 'bold','interp','none' );

% Overlay secondary structure
if exist( 'secstr', 'var' ) & ~isempty( secstr );
    if ( secstr == 1 ) secstr = r.structure; end;
    if ischar( secstr )
        % show base pairs as red squares
        hold on;
        bps = convert_structure_to_bps( secstr );
        for i = 1:size( bps, 1 );
            rectangle( 'position', [r.seqpos(bps(i,1))-0.5, r.seqpos(bps(i,2))-0.5, 1, 1], 'edgecolor','r' )
            rectangle( 'position', [r.seqpos(bps(i,2))-0.5, r.seqpos(bps(i,1))-0.5, 1, 1], 'edgecolor','r' )
        end
        hold off;
    else
        % is this in use? -- rhiju, 2017
        seq = secstr{1};
        str = secstr{2};
        offset = secstr{3};
        for i = 1:length(secstr{1})
            data_types{i} = num2str(i);
        end
        area_pred = generate_area_pred(seq, str, 0, data_types, length(secstr{1}));
        % in future, use sequence and secstruct from rdat and crop to correct size
        [x_pred, y_pred] = find(area_pred);
        x_pred = x_pred + offset;
        y_pred = y_pred + offset;
        plot(x_pred, y_pred, 'o', 'color', [1 0 1]); hold on;
    end
end

% Overlay tertiary structure contours
% Is this in use? Rhiju -- 2017
if exist( 'pdb', 'var' ) && ~isempty( pdb )
    pdbvar = pdb{1}; contours = pdb{2};
    ligpos = xrng;
    [D_sim, res_rad, res_hit, dist_matrix, pdbstruct] = get_simulated_data( pdbvar );
    
    if contours ~= 0
        % Define contour levels and colors
        contour_levels = [15,30];
        
        % Add legends (NOTE: Å = char(197) for Ångstrom?)
        legends = cell(1, length(contour_levels));
        if max(res_rad)-min(res_rad) == max(seqpos)-min(seqpos) && ...
                max(res_hit)-min(res_hit) == max(ligpos)-min(ligpos);
            res_rad = seqpos;
            res_hit = ligpos;
        end;

        if ~isempty( dist_matrix );
            dist_smooth = smooth2d( dist_matrix );
            for i = 1:length( contour_levels )
                contour(res_rad, res_hit, tril(dist_smooth), ...
                    contour_levels(i) * [1 1],...
                    'color',c2(i,:),...
                    'linewidth', 1.5 );
                legends{i} = [num2str(contour_levels(i)),' Ångstrom'];
                hold on;
            end
        end
        if ~isempty( legends ); legend( legends, 'location', 'BestOutside' ); end;
    end
    % Set axis limits (crop to ROI)
    axis( [min(res_rad)-0.5 max(res_rad)+0.5 min(res_hit)-0.5 max(res_hit)+0.5 ]);
else
    % axis( [min(seqpos)-0.5 max(seqpos)+0.5 min(ligpos)-0.5 max(ligpos)+0.5 ]);
end;

% % Rotate xticklabels and reposition
% xticklabel = get(gca,'XTickLabel');
% set(gca,'XTickLabel','');
% hxLabel = get(gca,'XLabel');
% set(hxLabel,'Units','data');
% xLabelPosition = get(hxLabel,'Position');
% y = xLabelPosition(2)-2;
% XTick = str2num(xticklabel)+1;
% y = repmat(y,length(XTick),1);
% hText = text(XTick,y,xticklabel,'fonts',fntsize);
% set(hText,'Rotation',90,'HorizontalAlignment','right');
% xlab = get(gca,'XLabel');
% set(xlab,'Position',get(xlab,'Position') + [0 8 0]);

if exist( 'save_path', 'var' ) && ~isempty( save_path );
    print( gcf, '-depsc2', '-loose', '-r300', save_path);
    fprintf( ['Created: ', save_path, '\n'] );
    hgsave( save_path );
end;

hold off;
set(gcf, 'PaperPositionMode','auto','color','white');

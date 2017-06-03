function [r,Z] = m2seq_summary( r, summary_title );
% m2seq_summary( r, summary_title );
%
%  useful plots for m2seq data:
%
% 1. raw data      
% 2. ref. secondary structure with 1D
% 3. Z-score      
% 4. Z-score with reference bps
%
% r             = [REQUIRED] RDAT filename or RDAT object. (can also be a cell of {r, r_nomod})
% summary_title = [OPTION] title for first plot (default: RDAT filename or name)
%
% (C) Rhiju Das, Stanford University, 2017

if ~exist( 'r', 'var' ) help( mfilename); return; end;
if ~exist( 'varna_fig', 'file' ); 
     fprintf( 'Please set up VARNA output via Biers: https://ribokit.github.io/Biers/\n' );
     return;
end
if iscell( r )
    r_nomod = r{2};
    r = r{1};
end
if ~isobject( r ) & ischar( r )
    r_name = r;
    r = read_rdat_file( r_name );
end
if ~exist( 'r_nomod', 'var' ) | length( r_nomod ) == 0
    r_nomod = [];
elseif ~isobject( r_nomod ) & ischar( r_nomod )
    r_nomod = read_rdat_file( r_nomod );
end
if ~exist( 'summary_title', 'var' )
    if exist( 'r_name', 'var' ) summary_title = r_name; else; r_name = r.name; end
end
clf
subplot( 2, 2, 1 );
m2seqplot( r );
subplot( 2, 2, 2 );
output_varna( '/tmp/temp.png', r.sequence, r.structure,[],[],r.offset,[],[],SHAPE_normalize(r.reactivity(:,1)) );
imshow( imread( '/tmp/temp.png' ) );

subplot( 2, 2, 3);
Z = output_Zscore_from_rdat( [], r, r_nomod, [], 1, 1 );
rZ = get_Zscore_rdat( r, Z );
scalefactor = [];
m2seqplot( rZ, [], scalefactor );

subplot( 2, 2, 4);
%print_mode = 0;
%cluster_z_scores( Z, r.structure, r.offset, print_mode );
m2net( Z, r );

subplot(2,2,1);
title( summary_title );
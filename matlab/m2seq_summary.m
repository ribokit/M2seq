function m2seq_summary( r, summary_title );
% m2seq_summary( r, summary_title );
%
%  useful plots for m2seq data:
%
% raw data      ref. secondary structure with 1D
%
% Z-score       Z-score with reference bps
%
%
%

if ~isobject( r ) & ischar( r )
    r_name = r;
    r = read_rdat_file( r_name );
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
Z = output_Zscore_from_rdat( [], r );
rZ = get_Zscore_rdat( r, Z );
scalefactor = [];
m2seqplot( rZ, [], scalefactor );
subplot( 2, 2, 4);
m2seqplot( rZ, [], scalefactor, [],  1 );

subplot(2,2,1);
title( summary_title );
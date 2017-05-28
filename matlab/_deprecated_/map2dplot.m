function scalefactor = map2dplot( R, save_path, scalefactor, titl, secstr, fntsize, pdb, c2 )
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

fprint( 'map2dplot is deprecated; use m2seqplot instead.' );
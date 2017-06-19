function r_new = resample_rdat( r, frac_downsample );
% r_new = resample_rdat( r, frac_downsample );
% Assumes M2-seq formatted RDAT.
%
% INPUTS 
%  r = RDAT object holding M2seq data. must have 'reads' line (only
%      available in M2seq output after June 2017).
%  frac_downsample
%     = number from 0.0 to 1.0 to downsample. [Default: 1.0, which
%         simulates bootstrapping, i.e., resample with replacement]
%
% OUTPUT
%  r_new = new RDAT object with resampled data.
%
% (C) Rhiju Das, Stanford University

r_new = [];
if ( nargin < 1 ) help( mfilename); return; end; 
if ~exist( 'frac_downsample', 'var' ); frac_downsample = 0.01; end;

reads = str2num( strjoin( get_tag( r, 'reads' ) ) ) * frac_downsample;
counts = r.reactivity * diag( reads  ) * frac_downsample;
r_new = r;

% resample from Poisson distribution:
d = poissrnd( counts ) * diag(1./reads)/frac_downsample;
d( find( isnan( d ) ) ) = 0;

% force 'diagonal' to be exactly one.
% alternatively should divide the 'row' by value at that position?
d( find( abs(r.reactivity - 1.000 ) < 1.0e-4 ) ) = 1;

r_new.reactivity = d;


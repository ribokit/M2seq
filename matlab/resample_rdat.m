function r_new = resample_rdat( r, frac_downsample );
% assumes M2-seq formatted RDAT
% (C) Rhiju Das, Stanford University
frac = 0.01;
reads = str2num( strjoin( get_tag( r, 'reads' ) ) ) * frac_downsample;
r_new = r;
counts = r.reactivity * diag( reads  ) * frac_downsample;

% resample from Poisson distribution:
d = poissrnd( counts ) * diag(1./reads)/frac_downsample;
d( find( isnan( d ) ) ) = 0;
d( find( abs(r.reactivity - 1.000 ) < 1.0e-4 ) ) = 1;
r_new.reactivity = d;


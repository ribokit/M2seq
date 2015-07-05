% trying to reproduce my 'impromptu' Map2D analysis when
% we read the Homan paper for journal club.
%  -- rhiju

pos_cutoff = 10; % nucleotide position to which alignment must extend
hit_cutoff = 10; % maximum number of hits to allow before recording.

% this is so ad hoc:
if ~exist( 'D', 'var' ) 
  %filename = 'PCR_pAadaptBp/out.simple';
  %filename = 'PCR_pAadaptBp/out_no11.simple';
  %filename = 'Ligate_TruSeq/out_no11.simple';
  %filename = 'homan.simple'; offset = 102;
  %filename = '../cDNA_Ligate_TruSeq/2_Map2D/RTB001/out.simple'; % clarence's run.
  %filename = '../RTB001/out_no11.simple'; % clarence's run.
  %filename = '../../..//cDNA_PCR_pAadaptBp/2_Map2D/RTB001/out_no11.simple';
  %filename = 'out.simple';   offset = 89;
  filename = 'TEST/out.simple';   offset = 89;
  fid = fopen( filename );
  D = textscan( fid, '%d%d%s' );
  fclose( fid );
  fprintf( 'Read %s\n', filename );
end

F0 = zeros( 202,202 );
coverage = zeros( 202, 202 );
num_hits = [];

% Record 2D hits.
for i = 1:length( D{3} );
  if ( mod( i, 10000 ) == 0 ); fprintf( 'Doing %d of %d\n',i,length(D{3}) ); end;
  idx = strfind( char( D{3}( i ) ), '1' );
  pos = int32( idx ) + int32( D{1}(i) ) - 1;
  if D{1}(i) <= pos_cutoff
    if length( pos ) <= hit_cutoff 
      num_hits = [num_hits, length( pos ) ];
      covered_pos = [D{1}(i) : D{2}(i)];
      % keep track of which residue pairs had a chance of being covered.
      coverage( covered_pos, covered_pos ) = coverage( covered_pos, covered_pos ) + 1; 
      for m = pos
	F0( m, pos ) = F0( m, pos ) + 1;
      end
    end
  end
end
% Normalize to coverage
F = F0./coverage;
F( isnan( F ) ) = 0;

fprintf( 'Mean number of hits: %5.1f\n', mean( num_hits ) );

resnum = [1:202]+offset;
image( resnum, resnum,  F/100 )
colormap( 1 - gray(100 ) );

%F = a.reactivity;
for k = 1:length( F )
  A(:,k) = F(:,k ) / sum( F([1:(k-3) (k+3):end],k) );
end


image( 1:length( A ), 1:length( A ),  A' *10000 );

rhiju_page_setup
gp = find( mod(resnum,10) == 0 );
set( gca,'xtick',gp, 'xticklabel',resnum( gp ) );
set( gca,'ytick',gp, 'yticklabel',resnum( gp ) );
title( filename, 'interp','none')
export_fig( [filename,'.pdf' ] );


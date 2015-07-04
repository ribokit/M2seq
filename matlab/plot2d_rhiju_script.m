% trying to reproduce my 'impromptu' Map2D analysis when
% we read the Homan paper for journal club.
%  -- rhiju
if ~exist( 'D', 'var' ) 
  %filename = 'PCR_pAadaptBp/out.simple';
  %filename = 'PCR_pAadaptBp/out_no11.simple';
  %filename = 'Ligate_TruSeq/out_no11.simple';
  %filename = 'homan.simple';
  %filename = '../cDNA_Ligate_TruSeq/2_Map2D/RTB001/out.simple'; % clarence's run.
  %filename = '../RTB001/out_no11.simple'; % clarence's run.
  %filename = '../../..//cDNA_PCR_pAadaptBp/2_Map2D/RTB001/out_no11.simple';
  filename = 'out.simple';
  fid = fopen( filename );
  D = textscan( fid, '%d%d%s' );
  fclose( fid );
  fprintf( 'Read %s\n', filename );
end

F = zeros( 202,202 );
num_hits = [];
pos_cutoff = 1; % nucleotide position to which alignment must extend
hit_cutoff = 10; % maximum number of hits to allow before recording.
for i = 1:length( D{3} );
  if ( mod( i, 10000 ) == 0 ); fprintf( 'Doing %d of %d\n',i,length(D{3}) ); end;
  idx = strfind( char( D{3}( i ) ), '1' );
  pos = int32( idx ) + int32( D{1}(i) ) - 1;
  if D{1}(i) <= pos_cutoff
    if length( pos ) <= hit_cutoff 
      num_hits = [num_hits, length( pos ) ];
      for m = pos
	F( m, pos ) = F( m, pos ) + 1;
      end
    end
  end
end
fprintf( 'Mean number of hits: %5.1f\n', mean( num_hits ) );

resnum = [1:202]+90;
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


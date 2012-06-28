#test case for span parameter

use strict;
use warnings;
use Test::More tests => 12;
use Test::Exception;

BEGIN { 
	use_ok('Bio::Draw::FeatureStack'); 
	use_ok('Bio::DB::SeqFeature::Store'); 
	use_ok('Bio::Graphics::Glyph::decorated_gene'); 
	use_ok('Bio::Graphics::Glyph::decorated_transcript'); 
	use_ok('Bio::Graphics'); 
};

my $gff = 't/data/gene_models.gff3';

lives_ok { figure1() }  'Generation of span test';

sub figure1
{
	my $output_basename = "t/images/span";
	
	my @gene_names = (qw (PF11_0023 PF10_0392));

	my @features = load_features(\@gene_names);

	my $feature_stack = new Bio::Draw::FeatureStack
	(
		-features => \@features,
		-glyph => 'decorated_gene',
		-panel_params => {
							-width => 1024,
							-pad_left => 80,
							-pad_right => 20
		},
		-span => 800,
		-verbose => 0,
		-seperator => 1,
		-glyph_params => {
							-description => sub {"span set to 800 to only show 5' end" },
							-utr_color   => 'white',
							-pad_bottom  => 10,
							-height      =>  12,
							-label_transcripts => 1,
							-label => 1,
							-decoration_label_color => "black",
							-decoration_visible => 1,
							-decoration_border => 1,
						 }
	);

	my $png = $feature_stack->png;
	ok ($png, "PNG $output_basename" );
		
	my $png_file = $output_basename.".png";
	system("rm $png_file") if (-e $png_file);
	open(IMG,">$png_file") or die "could not write to file $png_file";
	print IMG $png;
	close(IMG);		
	ok (-e $png_file, "$png_file" );
	
		(my $svg, my $map) = $feature_stack->svg(-image_map => 1);

		ok ($svg, "SVG $output_basename" );
		ok ($map, "image map $output_basename");
		
		my $svg_file = $output_basename.".svg";
		system("rm $svg_file") if (-e $svg_file);
		open(IMG,">$svg_file") or die "could not write to file $svg_file";
		print IMG $svg;
		close(IMG);		

		ok (-e $svg_file, "$svg_file" );
}

sub load_features
{
	my $gene_names = shift;
	
	my $store = Bio::DB::SeqFeature::Store->new
	(
	   	-adaptor => 'memory',
		-dsn => $gff 
	);    				
	
	my @features;
	foreach my $name (@$gene_names)
	{
		my ($f) = $store->features(-name => $name, -aliases => 1, -type => "gene:test");
		if (!defined $f)
		{
			die "could not load feature $name from gff file $gff";
			return ();
		}
		
		push(@features, $f);	
	}
	
	return @features;
}


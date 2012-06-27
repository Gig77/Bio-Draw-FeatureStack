##test case for genbank data source

use strict;
use warnings;
use Test::More tests => 9;
use Test::Exception;

BEGIN { 
	use_ok('Bio::Draw::FeatureStack'); 
	use_ok('Bio::DB::SeqFeature::Store'); 
	use_ok('Bio::Graphics::Glyph::decorated_gene'); 
	use_ok('Bio::Graphics::Glyph::decorated_transcript'); 
	use_ok('Bio::Graphics'); 
	use_ok('Bio::DB::GenBank');
	use_ok('File::Basename');
	use_ok('Bio::SeqFeature::Tools::Unflattener'); 
};

my $gff = 't/data/gene_models.gff3';

lives_ok { figure1() }  'Generation of genbank figure';

sub figure1
{
	my $output_basename = "t/images/genbank";
	
	my @gene_names = (qw (AH013929.2 Z48783.1));

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
		-glyph_params => {
							-description => 1,
							-label => 1,
							-sub_part => 'exon',  # Unflattener creates exon subparts, not CDS
							-decoration_visible => 1,
							-decoration_color => 'yellow'
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
	
}

sub load_features
{
	my $gene_names = shift;

	my @features;
	
	my $db_obj = Bio::DB::GenBank->new;
	my $unflattener = Bio::SeqFeature::Tools::Unflattener->new;

	foreach my $name (@$gene_names)
	{
		my $seq = $db_obj->get_Seq_by_acc($name);

		if (!defined $seq)
		{
			die "could not load sequence $name from genbank; maybe no internet connection?\n";
			return ();
		}

		# unflatten flat genbank features into hierarchical gene structure		
	  	$unflattener->unflatten_seq(-seq => $seq, -use_magic => 1);
	  	
	  	# grep gene and mRNA top-level features
		my @f = grep { $_->primary_tag =~ /(gene|mRNA)/ } $seq->top_SeqFeatures;
		
		# let's add some dummy decorations to first transcript...
		my $mrna = $f[0]->primary_tag eq 'mRNA' ? $f[0] : (grep { $_->primary_tag eq 'mRNA' } $f[0]->get_SeqFeatures)[0];
		$mrna->add_tag_value('protein_decorations', 'dummy:dummy:10:50:0') if ($mrna); 
		
		push(@features, @f);
	}
	
	return @features;
}		


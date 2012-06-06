# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl Bio-Draw-FeatureStack.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More tests => 1;
BEGIN { use_ok('Bio::Draw::FeatureStack') };

#########################

use strict;
use warnings;

use lib "$ENV{HOME}/perl";
use ChenLab::Constants;
use ChenLab::DB::DataProvider;
use ChenLab::DB::DataProviderBundle;
use ChenLab::DB::DataProviderFasta;
use ChenLab::DB::Util;
use ChenLab::GFC::WebPageWriter;
use Bio::Draw::FeatureStack;
use Bio::DB::GFF;
use Bio::Graphics;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Data::Dumper;
use ChenLab::GFC::WebPageWriter;
use File::Basename;
use Bio::Graphics::Glyph::decorated_transcript 0.10;

use constant DEBUG => 1;

my $dp_cele = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => 'data/gene_models.gff3/home/cfa24/scripts_data/feature_stack/gff/c_elegans.WS231.gff3.annotated',
	-default_track => 'Coding_transcript',
	-species_long_name => Constants::SPECIES_LONG_NAME_CELE,
	-species_short_name =>  Constants::SPECIES_NAME_CELE
);
my $dp_cbri = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/c_briggsae.WS231.gff3.annotated',
	-default_track => 'curated',
	-species_long_name => Constants::SPECIES_LONG_NAME_CBRI,
	-species_short_name =>  Constants::SPECIES_NAME_CBRI
);
my $dp_human = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/Homo_sapiens.GRCh37.66.gff3.annotated',
	-default_track => 'protein_coding',
	-species_long_name => Constants::SPECIES_LONG_NAME_HSAP,
	-species_short_name =>  Constants::SPECIES_NAME_HSAP
);
my $dp_dmel = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/Drosophila_melanogaster.BDGP5.66.gff3.annotated',
	-default_track => 'protein_coding',
	-species_long_name => Constants::SPECIES_LONG_NAME_DMEL,
	-species_short_name =>  Constants::SPECIES_NAME_DMEL
);
my $dp_scer = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/saccharomyces_cerevisiae.gff3.annotated',
	-default_track => 'SGD',
	-species_long_name => Constants::SPECIES_LONG_NAME_SCER,
	-species_short_name =>  Constants::SPECIES_NAME_SCER
);
my $dp_acas = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/Acanthamoeba_castellani.gff3.annotated',
	-default_track => 'TIGR',
	-species_long_name => Constants::SPECIES_LONG_NAME_ACAS,
	-species_short_name =>  Constants::SPECIES_NAME_ACAS
);
my $dp_amac = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/allomyces_macrogynus_atcc_38327_3_transcripts.gff3.annotated',
	-default_track => 'A_macrogynus_V3_CALLGENES_FINAL_2',
	-species_long_name => Constants::SPECIES_LONG_NAME_AMAC,
	-species_short_name =>  Constants::SPECIES_NAME_AMAC
);
my $dp_mbre = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/Monbr1_best_models.gff3.annotated',
	-default_track => 'JGI',
	-species_long_name => Constants::SPECIES_LONG_NAME_MBRE,
	-species_short_name =>  Constants::SPECIES_NAME_MBRE
);

my $dp_nematode = new ChenLab::DB::DataProvider
(
	-adaptor => 'memory', 
	-dsn => '/home/cfa24/scripts_data/feature_stack/gff/nematode.RNApol2genes.gff3.annotated',
	-default_track => 'Coding_transcript',
	-species_long_name => 'Nematodes'
);


my $dpb = new ChenLab::DB::DataProviderBundle(-data_providers => [$dp_cele, $dp_cbri, $dp_human, $dp_dmel, $dp_scer, $dp_acas, $dp_amac, $dp_mbre, $dp_nematode]);

t1();
t2();

sub render_figure_1
{
	my %param = @_;
	my $output_basename = "figure1";

	my @gene_names = (qw (RFX3 RFX2 RFX1 dRFX ceDAF-19 Cbr-daf-19 cRFX1 RFX4 RFX6 cRFX2 RFX8 ACA1_270030 YLR176C RFX5 RFX7 dRFX1 AMAG_11601));
	my @features = load_features(\@gene_names);
	#my @gene_names = (qw (dRFX dRFX1 ceDAF-19 YLR176C ACA1_270030 AMAG_10999 AMAG_11601 cRFX1 cRFX2));
	#my @gene_names = (qw (ceDAF-19));
	#my @gene_names = (qw (YLR176C));

	
	Log::log("render(): Initializing FeatureStack...") if (DEBUG);
	
	my $feature_stack = new Bio::Draw::FeatureStack
	(
		-features => \@features,
		-glyph => 'decorated_gene',
		-flip_minus => 1,
		-ignore_utr => 1,
		-panel_params => {
							-width => 1024,
							-pad_left => 80,
							-pad_right => 20
		},
		-intron_size => 1,
		-verbose => DEBUG,
		-feature_offsets => "DBD",
#		-feature_offsets => {
#			'daf-19' => 800,
#			'Cbr-daf-19' => 570,
#			'RFX1' => 190,
#			'RFX2' => 1060,
#			'RFX3' => 1160,
#			'RFX4' => 1500,
#			'RFX5' => 1450,
#			'RFX6' => 1400,
#			'RFX7' => 1400,
#			'RFX8' => 1750,
#			'dRFX1' => 950,
#			'YLR176C' => 1080,
#			'ACA1_270030' => 1480,
#			'AMAG_10999' => 1200,
#			'AMAG_11601' => 1200,
#			'cRFX1' => 950,
#			'cRFX2' => 640
#		},
		-glyph_params => {
							-bgcolor     => 'lightgrey',
							-fgcolor     => 'black',
							-fontcolor   => 'black',
							-font2color  => 'blue',
							-utr_color   => 'white',
							-bump        =>  +1,
							-height      =>  12,
							-label_position => 'left',
							-label_transcripts => 1,
							-label => \&get_gene_label, #1, #$gene_models_visible ? \&get_feature_name : 0,
							-description => 0, #$gene_models_visible ? sub { Util::getGeneDescription($_[0]) || "n/a" } : " ",
							-decoration_position => \&get_decoration_position,
							-decoration_color => \&get_color_for_decoration,
							-decoration_label => \&get_decoration_label,
							-decoration_label_color => "black", #\&get_label_color_for_decoration,
							-decoration_visible => \&is_decoration_visible,
							-decoration_border => \&get_decoration_border,
							-decoration_level => \&get_decoration_level,
							-decoration_label_position => \&get_decoration_label_position,
							-decoration_height => \&get_decoration_height,
							-additional_decorations => \&get_additional_decorations,
#							-sorted_decorations => \&get_sorted_decorations,
							-box_subparts => 3,
							-title => '$name', # \&ChenLab::GFC::WebPageWriter::get_decoration_title,
							-link => \&ChenLab::GFC::WebPageWriter::get_decoration_link
						 }
	);
	
	# write PNG
	my ($png, $map) = $feature_stack->svg(-image_map => 1);
	my $output_file = $output_basename.".svg";
	open(IMG,">$output_file") or Bio::Root::Exception("could not write to file $output_file");
	print IMG $png;
	close(IMG);
	Log::log("Stacked feature image successfully written to $output_file.");
	
	# write HTML (including image map)	
	my $img_file_base = basename($output_file);
	my $html_file =  $output_basename.".html";
	open(HTML, ">$html_file") or Bio::Root::Exception("could not write to file $html_file");
	print HTML "<html>\n<body>\n";
	print HTML "<img src=\"$img_file_base\" usemap=\"#map\" />\n";
	print HTML "$map";
	print HTML "</body>\n</html>\n";
	close(HTML);

	Log::log("HTML file written to $html_file");	
}

sub render_figure_2
{
	my %param = @_;
	my $output_basename = "figure2";
	
#	my @gene_names = (qw (che-13 xbx-1 dylt-2 xbx-3 xbx-4 xbx-5 xbx-6 mks-1 ZK328.7 bbs-9 che-11 odr-4 osm-5 nhr-44 nphp-1 nphp-4 nud-1 dyf-2 osm-6 dyf-3 che-2 osm-1 bbs-1 bbs-2 bbs-5 osm-12 bbs-8 tub-1 che-12 dyf-5));
	my @gene_names = (qw (ceDAF-19 xbx-1 dylt-2 xbx-3 xbx-4 xbx-5 xbx-6 mks-1 ZK328.7 bbs-9 che-11 odr-4 osm-5 nhr-44 nphp-1 nud-1 dyf-2 osm-6 dyf-3 che-2 osm-1 bbs-1 bbs-2 bbs-5 osm-12 bbs-8 tub-1 dyf-5));
#	my @gene_names = (qw (xbx-1 dylt-2 xbx-3 xbx-4 mks-1 bbs-9 osm-5 nhr-44 dyf-5));
	my @features = load_features(\@gene_names);

	Log::log("render(): Initializing FeatureStack...") if (DEBUG);
	
	my $feature_stack = new Bio::Draw::FeatureStack
	(
		-features => \@features,
		-transcripts_to_skip => [qw(Transcript:F33H1.1b Transcript:F33H1.1c Transcript:F33H1.1d 
									Transcript:M04C9.5a Transcript:R148.1b 
									Transcript:F40F9.1b.1 Transcript:F40F9.1a.2 Transcript:F40F9.1b.2 Transcript:F40F9.1a.3 Transcript:F40F9.1b.3 Transcript:F40F9.1a.4 Transcript:F40F9.1b.4
									Transcript:ZK328.7b
									Transcript:Y102E9.1a.2 Transcript:Y102E9.1a.3 Transcript:Y102E9.1a.4 Transcript:Y102E9.1b Transcript:Y102E9.1c
									Transcript:F53A2.4.2
									Transcript:C04C3.5b Transcript:C04C3.5c
									Transcript:F38G1.1.2
									Transcript:ZK520.3b.1 Transcript:ZK520.3b.2)],
		-alt_feature_type => 'xbox:Coding_transcript',
		-flip_minus => 1,
		-ignore_utr => 1,
		-panel_params => {
							-width => 1024,
							-pad_left => 170,
							-pad_right => 20,
							-grid => 1							
		},
		-span => 800,
		-separator => 0,
#		-intron_size => 50,
		-verbose => DEBUG,
		-feature_offsets => 'start_codon',
		-glyph => 'decorated_gene',
		-glyph_params => {
							-bgcolor     => 'lightgrey',
							-fgcolor     => 'black',
							-fontcolor   => 'black',
							-font2color  => 'blue',
							-utr_color   => 'white',
							-bump        =>  +1,
							-height      =>  10,
							-label_position => 'left',
							-label_transcripts => 1,
							-pad_left => 200,
							-label => sub { my $f = shift; $f->primary_tag eq 'mRNA' ? $f->name : 0 }, #\&get_gene_label, #1, #$gene_models_visible ? \&get_feature_name : 0,
							-description => 0, #$gene_models_visible ? sub { Util::getGeneDescription($_[0]) || "n/a" } : " ",
							-decoration_position => \&get_decoration_position,
							-decoration_color => \&get_color_for_decoration,
							-decoration_label => \&get_decoration_label,
							-decoration_label_color => "black", #\&get_label_color_for_decoration,
							-decoration_visible => 0,
							-decoration_border => \&get_decoration_border,
							-decoration_level => 0, #\&get_decoration_level,
							-decoration_label_position => \&get_decoration_label_position,
							-decoration_height => \&get_decoration_height,
							-additional_decorations => \&get_additional_decorations,
#							-sorted_decorations => \&get_sorted_decorations,
							-box_subparts => 3,
							-title => '$name', # \&ChenLab::GFC::WebPageWriter::get_decoration_title,
							-link => \&ChenLab::GFC::WebPageWriter::get_decoration_link
						 },
		-alt_glyph => 'generic', 
		-alt_glyph_params => {
							-bgcolor     => 'red',
							-fgcolor     => 'black',
							-fontcolor   => 'black',
							-font2color  => 'blue',
							-bump        =>  +1,
							-height      =>  10,
							-label_position => 'left',
							-label => \&get_xbox_label,
							-description => 0
						 }
	);
	
	# write PNG
	my ($png, $map) = $feature_stack->svg(-image_map => 1);
	my $output_file = $output_basename.".svg";
	open(IMG,">$output_file") or Bio::Root::Exception("could not write to file $output_file");
	print IMG $png;
	close(IMG);
	Log::log("Stacked feature image successfully written to $output_file.");
	
	# write HTML (including image map)	
	my $img_file_base = basename($output_file);
	my $html_file =  $output_basename.".html";
	open(HTML, ">$html_file") or Bio::Root::Exception("could not write to file $html_file");
	print HTML "<html>\n<body>\n";
	print HTML "<img src=\"$img_file_base\" usemap=\"#map\" />\n";
	print HTML "$map";
	print HTML "</body>\n</html>\n";
	close(HTML);

	Log::log("HTML file written to $html_file");	
}

sub get_gene_label
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	return $feature->name if ($feature->primary_tag eq 'mRNA');
	return $feature->name if ($feature->name and $feature->name eq 'ceDAF-19');
	return 0;
}

sub get_xbox_label
{
	my $feature = shift;
	my ($seq,$len,$score) = $feature->desc =~ /(.*)-Length:(.*)-Score:(.*)/; 
	my ($dist) = $feature->get_tag_values('start_dist');  # dynamically computed by FeatureStack; no need to set in GFF file
	return $seq."(".$dist."bp;".$score.")";
}

sub get_decoration_level
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	
	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	return "mRNA" if ("$type:$name" eq "hmmer3:Pox_D5");
#	return "mRNA" if ("$type:$name" eq "hmmer3:BCD-domain");
	return "CDS";
}

sub get_additional_decorations
{
	my $feature = shift;	
	my %additional_decorations =
	(
#		"MAL7P1.59" => "test:test:24:24:0",
	);
	
	my $feature_id = get_feature_name($feature);
	$feature_id =~ s/rna_//;
	$feature_id =~ s/-1$//;
	
	my $add_decorations = $additional_decorations{$feature_id};
	
#	if ($hmmer{$feature_id})
#	{
#		$add_decorations .= "," if ($add_decorations);
#		$add_decorations .= $hmmer{$feature_id};
#	}
	
	return $add_decorations;
}

sub get_color_for_decoration
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	
	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	my %decoration_colors =
	(
		'SEG:LCR' => 'white',
		'hmmer3:DBD' => 'black',
		'hmmer3:A-domain' => 'red',
		'hmmer3:B-domain' => 'yellow',
		'hmmer3:C-domain' => 'blue',
		'hmmer3:D-domain' => 'lawngreen',
		'hmmer3:BCD-domain' => 'darkgray',
		'hmmer3:RFX1_trans_act' => 'darkslateblue',
		'hmmer3:Pox_D5' => 'darkgreen'
	);

	return $decoration_colors{"$type:$name"} 
		if (exists $decoration_colors{"$type:$name"});

	return $decoration_colors{$name} 
		if (exists $decoration_colors{$name});

	return $decoration_colors{$type} 
		if (exists $decoration_colors{$type});
		
	if ($type eq 'blastp')
	{
		if ($h->score < 1e-100)
		{
			return "red";
		}
		elsif ($h->score < 1)
		{
			return "blue";
		}
		else
		{
			return "black";
		}
	}
}

sub get_decoration_border
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

#	if ("$type:$name" eq "hmmer3:BCD-domain")
#	{
#		return "dashed";
#	}
	
	return "none";
}

sub get_decoration_height
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	return 5 if ("$type:$name" eq "hmmer3:Pox_D5");
	return 6 if ("$type:$name" =~ /hmmer3:[ABCD]-domain/);
#	return 2 if ("$type:$name" eq "hmmer3:BCD-domain");
#	return 10;
}

sub get_decoration_position
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	if ("$type:$name" eq "hmmer3:Pox_D5")
	{
		return 14;
	}
#	if ("$type:$name" eq "hmmer3:BCD-domain")
#	{
#		return 5;
#	}

	return "inside";
}

sub get_decoration_label_position
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

#	return "above" if ("$type:$name" eq "hmmer3:BCD-domain");
	return "inside" if ($name eq "LCR");
	return "below";
}

sub get_decoration_label
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	return 0 if ("$type:$name" eq "hmmer3:DBD");
	return "A" if ("$type:$name" eq "hmmer3:A-domain");
	return "B" if ("$type:$name" eq "hmmer3:B-domain");
	return "C" if ("$type:$name" eq "hmmer3:C-domain");
	return "D" if ("$type:$name" eq "hmmer3:D-domain");
	return 0 if ("$type:$name" eq "hmmer3:BCD-domain");
	return 0 if ("$type:$name" eq "hmmer3:RFX1_trans_act");
	return 0 if ("$type:$name" eq "hmmer3:Pox_D5");
	
}

sub get_label_color_for_decoration
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name) = ($h->type, $h->name);

	my %decoration_label_colors =
	(
		'SP' => 'black',
		'TM' => 'white',
		'VTS' => 'white',
		'Phobius:TM' => 'white'
	);	
		
	return $decoration_label_colors{"$type:$name"} 
		if (exists $decoration_label_colors{"$type:$name"});

	return $decoration_label_colors{$name} 
		if (exists $decoration_label_colors{$name});

	return $decoration_label_colors{$type} 
		if (exists $decoration_label_colors{$type});
}

sub is_decoration_visible
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my $h = $glyph->active_decoration;
	my ($type, $name, $start, $end, $score) = ($h->type, $h->name, $h->start, $h->end, $h->score);

	return 0 if ("$name" eq "RFX_DNA_binding");
#	return 0 if ("$name" eq "RFX1_trans_act");
	
	return 1;
}

sub get_sorted_decorations
{
	my ($feature, $option_name, $part_no, $total_parts, $glyph) = @_;	

	my @sorted_decorations = sort { $a->length <=> $b->length } (@{$glyph->mapped_decorations});
	return \@sorted_decorations;
}

sub get_feature_name {
	my $feature = shift;

	my $name = Util::getFeatureName($feature);
	$name = Util::getFeatureName(($feature->get_SeqFeatures)[0]) if (!$name and $feature->get_SeqFeatures);
	if (defined $name)
	{
		$name = $3 if ($name =~ /^(m?(RNA|CDS)_)?(.*?)(-\d+)?$/i);		
	}
	else
	{
		$name = "";	
	}

	return $name;
}

sub add_species_prefix
{
	my $id = shift;
	return "$id";
}

sub load_features
{
	my $gene_names = shift;

	my @features;
	foreach my $name (@$gene_names)
	{
		if (!defined $name)
		{
			push(@features, undef);
			next;
		}
		
		my $f = $dpb->containsFeature($name);
		die "feature $name not found in data provider. Skipped\n" if (!$f);
		push(@features, $f);	
	}
	
	return @features;
}



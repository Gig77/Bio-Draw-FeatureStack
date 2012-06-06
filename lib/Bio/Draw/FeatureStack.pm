package Bio::Draw::FeatureStack;

use 5.008008;
use strict;
use warnings;

#require Exporter;

#our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use Bio::Draw::FeatureStack ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
#our %EXPORT_TAGS = ( 'all' => [ qw(
#	
#) ] );
#
#our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
#
#our @EXPORT = qw(
#	
#);

our $VERSION = '0.01';

use Bio::Graphics::Feature;
use Bio::Graphics::Panel;
use List::Util qw[min max];

use base qw(Bio::Root::Root);

use constant DEBUG => 0;

sub new
{
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
	    
	my %param = @args;
	my $features = $param{'-features'} or throw Bio::Root::BadParameter("gene models not specified");  # Bio::DB::SeqFeature array ref
	my $transcripts_to_skip = $param{'-transcripts_to_skip'}; # optional
	my $alt_feature_type = $param{'-alt_feature_type'}; # optional
	my $glyph = defined $param{'-glyph'} ? $param{'-glyph'} : "generic";
	my $alt_glyph = defined $param{'-alt_glyph'} ? $param{'-alt_glyph'} : "generic";
	my $glyph_params = defined $param{'-glyph_params'} ? $param{'-glyph_params'} : {};
	my $panel_params = defined $param{'-panel_params'} ? $param{'-panel_params'} : {};
	my $alt_glyph_params = defined $param{'-alt_glyph_params'} ? $param{'-alt_glyph_params'} : {};
	my $ignore_utr = $param{'-ignore_utr'};
	my $separator = $param{'-separator'};
	my $intron_size = $param{'-intron_size'};
	my $flip_minus = $param{'-flip_minus'};
	my $feature_offsets = $param{'-feature_offsets'};
	my $ruler = defined $param{'-ruler'} ? $param{'-ruler'} : 1;
	my $span = $param{'-span'}; # optional
				
	$self->{'features'} = $features;
	$self->{'transcripts_to_skip'} = $transcripts_to_skip;
	$self->{'glyph'} = $glyph;
	$self->{'alt_glyph'} = $alt_glyph;
	$self->{'glyph_params'} = $glyph_params;
	$self->{'panel_params'} = $panel_params;
	$self->{'alt_glyph_params'} = $alt_glyph_params;
	$self->{'alt_feature_type'} = $alt_feature_type;
	$self->{'ignore_utr'} = $ignore_utr;
	$self->{'separator'} = $separator;
	$self->{'intron_size'} = $intron_size;
	$self->{'flip_minus'} = $flip_minus;
	$self->{'feature_offsets'} = undef;
	$self->{'ruler'} = $ruler;
	$self->{'span'} = $span;

	# set glyph param defaults
	$glyph_params->{'-height'} = 12 if (!defined $glyph_params->{'-height'});

	my @glyph_params;
	map { push(@glyph_params, ($_, $self->{'glyph_params'}->{$_})) } keys(%{$self->{'glyph_params'}});
	$self->{'glyph_param_array'} = \@glyph_params;
		
	my @alt_glyph_params;
	map { push(@alt_glyph_params, ($_, $self->{'alt_glyph_params'}->{$_})) } keys(%{$self->{'alt_glyph_params'}});
	$self->{'alt_glyph_param_array'} = \@alt_glyph_params;

	my @panel_params;
	map { push(@panel_params, ($_, $self->{'panel_params'}->{$_})) } keys(%{$self->{'panel_params'}});
	$self->{'panel_param_array'} = \@panel_params;

	# feature transformation: adjust intron sizes, remove UTRs, flip if on negative strand, remove unwanted transcripts (isoforms)
	my @f = @$features;
	if ($ignore_utr or (defined $intron_size and $intron_size > 1))
	{
		my @transformed_features;
		foreach my $feature (@{$features})
		{
			$self->throw("feature undefined") if (!defined $feature);
			my $transformed = $self->_transform_feature($feature);
			push(@transformed_features, $transformed);
		}
		@f = @transformed_features;
	}

	# calculate feature offsets if not left-aligned
	$self->_calc_feature_offsets(\@f, $feature_offsets)
		if (defined $feature_offsets);
			
	# determine coordinate span
	$span = $self->_calc_feature_span(\@f)
		if (!defined $self->{'span'}); 
	
	# align features
	my @aligned_features;
	for (my $i = 0; $i < @f; $i ++)
	{
		my $feature = $f[$i];
		
		my $offset = (defined $self->{'feature_offsets'} and exists $self->{'feature_offsets'}->{$feature->id}) ? $self->{'feature_offsets'}->{$feature->id} : 0;
		print $feature->id.": feature offset=$offset\n" if ($offset and DEBUG);
		
		my $aligned_feature = $self->_align_feature
		(
			$feature, 
			-$feature->start+$offset+1
		);
		$self->throw("aligned_feature undefined: ".$feature->name."\n") if (!$aligned_feature);
		push(@aligned_features, $aligned_feature);
	}

	$self->{'aligned_features'} = \@aligned_features;
	$self->{'aligned_alt_features'} = $alt_feature_type ? $self->_get_alt_features(\@aligned_features) : []; 
		
	return bless($self, $class);
}

sub _get_alt_features
{
	my $self = shift;
	my $features = shift;
	
	my ($alt_type, $alt_source) = split(':', $self->{'alt_feature_type'});
	
	$self->throw("could not parse type of alternative feature: ".$self->{'alt_feature_type'})
		if (!$alt_type or !$alt_source);

	my @alt_features;	
	foreach my $feature (@$features)
	{
		my @afs = grep {$_->primary_tag eq $alt_type and $_->source eq $alt_source} $feature->get_SeqFeatures();
		foreach my $af (@afs)
		{
			my $dist = $self->_calc_start_dist($feature, $af);
			$af->add_tag_value('start_dist', $dist);
			
		}
		push(@alt_features, \@afs);
	}	
	
	return \@alt_features;
}

# computes distance of alternative feature (e.g. DNA-binding site) to start of nearest transcript, in bp
sub _calc_start_dist
{
	my $self = shift;
	my $f = shift;
	my $af = shift;
	
	return $af->start-$f->start
		if ($f->primary_tag eq 'mRNA');

	my @transcripts = grep {$_->primary_tag eq 'mRNA'} $f->get_SeqFeatures();
	if (@transcripts == 0)
	{
		print "WARNING: No transcripts found for feature $f (".$f->name.")\n" if (DEBUG);
		return undef;
	}	
		
	my $nearest_transcript;
	foreach my $t (@transcripts)
	{
		$nearest_transcript = $t
			if (!$nearest_transcript or abs($af->start-$nearest_transcript->start)>abs($af->start-$t->start));
	}
	
	return $af->start-$nearest_transcript->start;
}

# removes UTRs, shrinks introns, flip negative strand, remove unwanted transcripts (isoforms)
sub _transform_feature
{
	my $self = shift;
	my $feature = shift;
	
	my @transcripts;
	if ($feature->primary_tag eq 'gene')
	{
		# fixed intron size mode currently only functional with mRNA features (not promoters)
		if ($self->{'intron_size'})
		{
			push(@transcripts, grep {$_->primary_tag =~ /mRNA/i} $feature->get_SeqFeatures());		
		}
		else
		{
			push(@transcripts, $feature->get_SeqFeatures());					
		}
	}
	elsif ($feature->primary_tag eq 'mRNA')
	{
		push(@transcripts, $feature);
	}
	else
	{
		$self->throw("invalid feature type: ".$feature->primary_tag);
	}
	
	# remove unwanted isoforms
	my %transcripts_to_skip;
	map {$transcripts_to_skip{$_} = 1} @{$self->{'transcripts_to_skip'}}
		if ($self->{'transcripts_to_skip'});
		
	@transcripts = grep {!$transcripts_to_skip{_get_id($_)}} @transcripts;
	
	my @shifted_transcripts;
	foreach my $transcript (@transcripts)
	{
		my $shift_by = 0;
		my $last_end;
		my @shifted_parts;
		my @parts;
		
		if ($self->{'ignore_utr'})
		{
		 	@parts = sort {$a->start <=> $b->start} grep {$_->primary_tag !~ /utr/i} $transcript->get_SeqFeatures();			
		}
		else
		{
		 	@parts = sort {$a->start <=> $b->start} $transcript->get_SeqFeatures();
		}
		
		$self->throw("not subfeature found for transcript ".$transcript->id."\n") 
			if ($transcript->primary_tag eq 'mRNA' and @parts == 0);
			
		foreach my $part (@parts)
		{
			if ($self->{'intron_size'} and $part->primary_tag =~ /utr|cds/i)
			{
				print "shift_by=$shift_by\n" if (DEBUG == 2);
				if (defined $last_end)
				{
					my $intron_size = $part->start-$last_end;
					$shift_by += $intron_size - $self->{'intron_size'} if ($intron_size > 1);  # shift by difference to maximum allowed intron size
				}				
				$last_end = $part->end;
			}
							
			my ($f, $id) = $self->_clone_feature
			(
				$part, 
				($self->{'flip_minus'} and $feature->strand < 1) ? $feature->start + $feature->end - $part->end + $shift_by : $part->start - $shift_by, 
				($self->{'flip_minus'} and $feature->strand < 1) ? $feature->start + $feature->end - $part->start + $shift_by : $part->end - $shift_by,
				$self->{'flip_minus'} ? 1 : undef
			);		
			push(@shifted_parts, $f);
		}
				
		@shifted_parts = sort {$a->start <=> $b->start} (@shifted_parts);
		my ($shifted_transcript, $id) = $self->_clone_feature
		(
			$transcript, 
			@shifted_parts > 0 ? $shifted_parts[0]->start 
							   : ($self->{'flip_minus'} and $feature->strand < 1) ? $feature->start + $feature->end - $transcript->end 
							   													: $transcript->start, 
			@shifted_parts > 0 ? $shifted_parts[@shifted_parts-1]->end 
							   : ($self->{'flip_minus'} and $feature->strand < 1) ? $feature->start + $feature->end - $transcript->start 
							                                                      : $transcript->end,
			$self->{'flip_minus'} ? 1 : undef
		);
		foreach my $c (@shifted_parts)
		{
			$shifted_transcript->add_SeqFeature($c);
		}
		push(@shifted_transcripts, $shifted_transcript);
	}

	if ($feature->primary_tag eq 'gene')
	{
		@shifted_transcripts = sort {$a->start <=> $b->start} (@shifted_transcripts);
		my ($gene, $id) = $self->_clone_feature
		(
			$feature, 
			$shifted_transcripts[0]->start,
			$shifted_transcripts[@shifted_transcripts-1]->end, 
			$self->{'flip_minus'} ? 1 : undef
		);
		foreach my $i (@shifted_transcripts)
		{
			$gene->add_SeqFeature($i);
		}
		return $gene;
	}
	
	return $shifted_transcripts[0];
}

sub _get_id
{
	my $feature = shift;
	
	my ($id) = $feature->get_tag_values('ID');
	($id) = $feature->get_tag_values('load_id') if (!$id);
	$id = $feature->id if ($feature->can('id') and !$id);
	
	return $id;
}

sub _calc_feature_span
{
	my $self = shift;
	my $features_ref = shift;

	my $max_span = 0;
	foreach my $feature (@$features_ref)
	{
		next if (!defined $feature);
  		my $span = abs($feature->end - $feature->start + 1);	
		$span +=  $self->{'feature_offsets'}->{$feature->id} 
			if (defined $self->{'feature_offsets'} and exists $self->{'feature_offsets'}->{$feature->id});	
		$max_span = $span if ($span > $max_span);		
	}

	$self->{'span'} = $max_span;

	return $max_span;	
}

sub _calc_feature_offsets
{
	my $self = shift;
	my $features_ref = shift;
	my $feature_offsets = shift;

	# user-defined feature offsets? 
	if (ref($feature_offsets) eq "HASH")
	{
		$self->{'feature_offsets'} = $feature_offsets;
	}
	else
	{
		my @features = @$features_ref;
		$self->{'max_offset'} = 0;
		my %f_offset;
				
		# determine maximum offset by decoration position
		foreach my $feature (@$features_ref)
		{
			my @transcripts;
			if ($feature->primary_tag eq 'gene')
			{
				push(@transcripts, grep {$_->primary_tag eq 'mRNA'} $feature->get_SeqFeatures());
			}
			else
			{
				push(@transcripts, $feature);
			}
	
			foreach my $t (@transcripts)
			{
				if ($feature_offsets eq "start_codon") 
				{
					# align by start codon
					my @cds = sort {$a->start <=> $b->start} grep {$_->primary_tag eq 'CDS'} $t->get_SeqFeatures();
					$f_offset{$feature->id} = $cds[0]->start-$feature->start
						if (!defined $f_offset{$feature->id} or $f_offset{$feature->id} > $cds[0]->start-$feature->start);
				}
				else
				{
					# align features by decoration
					# requires Bio::Graphics::Glyph::decorated_transcript for coordinate mapping
					use Bio::Graphics::Glyph::decorated_transcript;
					
					my @decorations = Bio::Graphics::Glyph::decorated_transcript::get_decorations_as_features($t);
					foreach my $decoration (@decorations)
					{
						next if ($decoration->name ne $feature_offsets);
						$f_offset{$feature->id} = $decoration->start - $feature->start; 
					}																			
				}
				$self->{'max_offset'} = $f_offset{$feature->id} if ($f_offset{$feature->id} > $self->{'max_offset'});
			}
		}
		
		# set offset for transcripts with this decoration
		foreach my $f_id (keys(%f_offset))
		{
			$self->{'feature_offsets'}{$f_id} = $self->{'max_offset'} - $f_offset{$f_id}; 
		}
	}	
}

sub _align_feature
{
	my $self = shift;
	my $feature = shift or throw Bio::Root::BadParameter("feature not specified");
	my $offset = shift;
	my $parent_id = shift;  # internal parameter for recursive calls

	$self->throw("offset not specified") 
		if (!defined $offset);

	my ($start, $end) = ($feature->start, $feature->end);	
	my ($aligned_feature, $id) = $self->_clone_feature($feature, $start + $offset, $end + $offset);	
	$parent_id = $id if (!$parent_id);
#	$aligned_feature->add_tag_value('parent_id', $parent_id) if ($parent_id);	

	# copy subfeatures recursively
	foreach my $subfeature ($feature->get_SeqFeatures())
	{
		$aligned_feature->add_SeqFeature
		(
			$self->_align_feature($subfeature, $offset, $parent_id) 
		);
	}

	return $aligned_feature;
}

sub _clone_feature
{
	my $self = shift;
	my $feature = shift;
	my $start = shift;
	my $end = shift;
	my $strand = shift;
	
	my $clone = Bio::Graphics::Feature->new
	(
		-name => $feature->can("display_name") ? $feature->display_name : $feature->name,
		-seq_id => $feature->seq_id,
		-source => $feature->source,
		-primary_tag => $feature->primary_tag,
		-start => $start ? $start : $feature->start,
		-end => $end ? $end : $feature->end,
		-score => $feature->score,
		-strand => defined $strand ? $strand : $feature->strand,
		-frame => $feature->phase,
		-phase => $feature->phase
	);			

	# copy tags
	foreach my $tagname ($feature->get_all_tags())
	{
		next if (uc($tagname) eq 'LOAD_ID');
#		next if (uc($tagname) eq 'NAME');
#		next if (uc($tagname) eq 'PARENT_ID');
		next if (uc($tagname) eq 'ID');
		foreach my $value ($feature->get_tag_values($tagname))
		{
			$clone->add_tag_value($tagname, $value);
		}
	}

	my ($id) = $feature->get_tag_values('ID');
	($id) = $feature->get_tag_values('load_id') if (!$id);
	$clone->add_tag_value('ID', $id) if ($id);

	return ($clone, $id);	
}

sub _render_panel
{
	my $self = shift;
	my $panel = shift;

	# add ruler (or spacer)
	if ($self->{'ruler'})
	{
		$panel->add_track
		(
			Bio::Graphics::Feature->new(-start => 1, -end => $self->{'span'}),
			-glyph  => 'arrow',
			-fgcolor => 'black',
			-bump   => 0,
			-double => 1,
			-tick   => 2,
			-relative_coords => $self->{'max_offset'} ? 1 : 0,
			-relative_coords_offset => -$self->{'max_offset'}
		);		
	}
	else
	{
		$panel->add_track
		(
			Bio::Graphics::Feature->new(-start => 0, -end => 0),
			-glyph  => 'segments',
			-fgcolor => 'white',
			-height => 19
		);				
	}

	# add tracks for all features
	for (my $i = 0; $i < @{$self->{'aligned_features'}}; $i ++)
	{
		# render track with alternative feature above main track (if specified)
		my $alt_f = $self->{'aligned_alt_features'}->[$i];
		if ($alt_f)
		{
			$panel->add_track
			(
				$alt_f,
				-glyph => $self->{'alt_glyph'},
				@{$self->{'alt_glyph_param_array'}}
			);						
		}		

		# add main feature track
		my $f = $self->{'aligned_features'}->[$i];
		print "adding track ".$f->name."\n" if (DEBUG);		
		$panel->add_track
		(
			$f,
			-glyph => $self->{'glyph'},
			@{$self->{'glyph_param_array'}}
		);
		
		# separator
		if ($self->{'separator'} and $i < @{$self->{'aligned_features'}}-1)
		{	
			$panel->add_track
			(
				Bio::SeqFeature::Generic->new(-start => 1, -end => $self->{'span'}),
				-glyph  => 'line',
				-height => 1,
				-fgcolor => 'black',
				-bgcolor => 'black'
			);
		}		
	}	
}

sub png
{
	my $self = shift;
	my %param = @_;
	my $image_map = defined $param{'-image_map'} ? $param{'-image_map'} : 0;
	
	my $panel = Bio::Graphics::Panel->new
	(
		-length    => $self->{'span'}+1,
		-key_style => 'between',
		@{$self->{'panel_param_array'}}
	);
	$self->_render_panel($panel);
	
	my $png = $panel->png;
	
	if ($image_map)
	{
		my $map = $panel->create_web_map();
		$panel->finished();
		return ($png, $map);
	}
	
	$panel->finished();
	return $png;
}

sub svg
{
	my $self = shift;
	my %param = @_;
	my $image_map = defined $param{'-image_map'} ? $param{'-image_map'} : 0;
	
	use GD::SVG 0.32;  # minimum version 0.32 for correct color management (sub colorAllocateAlpha)
	
	my $panel = Bio::Graphics::Panel->new
	(
		-image_class => "SVG",
		-length    => $self->{'span'}+1,
		@{$self->{'panel_param_array'}}
	);

	$self->_render_panel($panel);
	
	my $svg = $panel->gd;
	$svg = _consolidate_svg($svg);

	if ($image_map)
	{
		my $map = $panel->create_web_map();
		$panel->finished();
		return ($svg, $map);
	}
	
	$panel->finished();
	return $svg;
}

# agreed, this is a bit of a hack to fix font alignment problems, but GBrowse does the same
# see GBrowse/cgi-bin/gbrowse_img
sub _consolidate_svg
{
	my $g = shift;

	my $height    = ($g->getBounds)[1];
	my $width    += ($g->getBounds)[0];

    my $image_height = $height;
    
    my $svg = qq(<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n\n);
    $svg   .= qq(<svg height="$image_height" width="$width" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n);

	my $offset = 0;
    my $s = $g->svg;
    my $current_width = 0;
    foreach (split "\n",$s) {
		if (m!</svg>!) {
		    last;
		}
		elsif (/<svg.+width="([\d.]+)"/) {
		    $current_width = int($1+0.5);
		    my $height     = $height - ($g->getBounds)[1];
			    $svg .= qq(<g transform="translate($offset,$height)">\n);
		}
		elsif ($current_width) {
		    $svg .= "$_\n";
		}
	}
	$svg .= "</g>\n" if $current_width;
	$offset += $current_width;
	$svg   .= qq(</svg>\n);

    # munge fonts slightly for systems that don't have Helvetica installed
    $svg    =~ s/font="Helvetica"/font="san-serif"/gi;
    $svg    =~ s/font-size="11"/font-size="9"/gi;  
    $svg    =~ s/font-size="13"/font-size="12"/gi;  

    return $svg;	
}


1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Bio::Draw::FeatureStack - Perl extension for blah blah blah

=head1 SYNOPSIS

  use Bio::Draw::FeatureStack;
  blah blah blah

=head1 DESCRIPTION

Stub documentation for Bio::Draw::FeatureStack, created by h2xs. It looks like the
author of the extension was negligent enough to leave the stub
unedited.

Blah blah blah.

=head2 EXPORT

None by default.



=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

A. U. Thor, E<lt>cfa24@localdomainE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by A. U. Thor

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.8 or,
at your option, any later version of Perl 5 you may have available.


=cut

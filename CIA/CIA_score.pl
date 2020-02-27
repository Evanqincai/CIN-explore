#! /usr/bin/perl

use strict;
use warnings;
use File::Basename qw(basename dirname);
use Statistics::Descriptive;
use Math::Round;
##########################################################
if(@ARGV != 3){print "perl CIA.pl sample_id *cna.seg.txt output_dir\n";die;};
if (!-d  $ARGV[2]){mkdir $ARGV[2]};

########################################################## Z_value
my @z_score;
open (CNA,$ARGV[1]) or die $!;
open (OUT,">$ARGV[2]/$ARGV[0].Z_value.txt") or die $!;
while (<CNA>){
	chomp;
	if (/^chr/){print OUT "$_\tZ_value\n";next;};
	my @tmp = split(/\s+/,$_);
	next if ($tmp[9] eq 'NA');
	my $copy_num_log2 = log($tmp[9]/2)/log(2);
	my $copy_num_log2_abs = abs($copy_num_log2);
	my $copy_num_log2_abs_sqrt = sqrt($copy_num_log2_abs);
	push @z_score ,$copy_num_log2_abs_sqrt;
	print OUT "$_\t$copy_num_log2_abs_sqrt\n";

}
close CNA;
close OUT;

########################################################## CIA score
my $qualitle = &mid(@z_score);
my ($qualitle95,$qualitle99) = split('_',$qualitle);
#print($qualitle95,$qualitle99);
`less $ARGV[2]/$ARGV[0].Z_value.txt |awk  '\$11 <= $qualitle99 && \$11 >=$qualitle95' |awk \'{a+=\$11}END{print \"$ARGV[0]\t\"a}\' >$ARGV[2]/$ARGV[0].CIA_score.txt`;
##########################################################

sub mid {
        my @list = sort @_;
        my $count = @list;
        if ($count == 0){
                return 0;die;
        }
        my $up = $list[round(($count-1)*0.95)];
        my $down = $list[round(($count-1)*0.99)];
        my $th = "$up\_$down";
        return $th;
}




#!/usr/bin/perl

my @pos;
my @name;
my @ORDER;

#warn "format_primers parameters: ",join("\t",@ARGV),"\n";

my $tail2="";
my $tail3="";
my $basename;
if (@ARGV > 3) {
   $basename = "$ARGV[3]$ARGV[4]";
   $tail2=$ARGV[6];
   $tail3=$ARGV[7];
}

$header = <STDIN>;
$header =~s/>//;

local $/="//";
$seq = <STDIN>;

exit unless ( $seq=~/\[/ );

my @dir = split/\//, $ARGV[0];
pop(@dir);
my $output = join("/", @dir);
$output = './';
open OUT, ">>$output/primers.txt" or die "cannot open outfile\n$!\n";

print "#############################################\n";
print "\n        $header\n";
print "#############################################\n\n";

my $column1 = substr $header, 0, -1 if $header =~ /\n$/;

print OUT $column1 . "\t";

$header=~s/\.matches\.Contig/\_/;
$header=~s/ .*//;
$prefix=$header;
$prefix=~s/[> \n\r]//g;
$prefix=$basename if $basename;

$seq=~s/[ 0-9\n\[\]\/\r]//g;

#for (0..length($seq)){ $pos[$_]=" "; $name[$_]=" ";}
#$seq=~s/(.{70})/$1\n/g;
#print $header,$seq;

local $/="\n";
<STDIN>;
<STDIN>;

 print "\n--------------\n\n"; 

while ($amp=<STDIN>){
 if ($basename) {
   $amp=~/AMPLICON ([A-Z]) /;
   $prefix="$basename-$1";
 }
 <STDIN>; 
 $pcr1size = <STDIN>;
 $f1=<STDIN>;
 @pcr1F = split(/ +/,$f1);
 $r1=<STDIN>;
 @pcr1R = split(/ +/,$r1);
 $pcr1par = <STDIN>;
 $trash=<STDIN>; 
 $trash=<STDIN>;
 $pcr2size = <STDIN>;
 $f2=<STDIN>;
 @pcr2F = split(/ +/,$f2);
 $r2=<STDIN>;
 @pcr2R = split(/ +/,$r2);
 $pcr2par = <STDIN>;
 <STDIN>;
 <STDIN>;


 if (($f1=~/^NO/) || ($f2=~/^NO/)){
   print $amp,"\n",$pcr1size,$f1,$r1,$pcr1par,"\n",$pcr2size,$f2,$r2,$pcr2par,"\n";
   print "\n--------------\n\n"; 
   next;
 }

 $from=$pcr1F[3]-20;
 $from=0 if $from < 0;
 $ampseq=substr($seq,$from,$pcr1R[3]-$pcr1F[3]+$pcr1R[4]+40);
 undef @pos;
 undef @name;
 for (0..length($ampseq)){ $pos[$_]=" "; $name[$_]=" "; }

 $from=$pcr1F[3]-19;
 $from=1 if $from < 0;
 fill_pos($pcr2F[4]-1,$pcr2F[3]-$from,">") if $f2 ne "\n";
 fill_pos($pcr2R[4]-1,$pcr2R[3]-$from,"<") if $r2 ne "\n"; 
 fill_pos($pcr1F[4]-1,$pcr1F[3]-$from,">");
 fill_pos($pcr1R[4]-1,$pcr1R[3]-$from,"<");
 
 my ($pr1name, $pr2name, $pr2name, $pr4name);
 
 if ($basename) {
   $pr1name = fill_name($pcr1F[3],"1",$from);
   $pr4name = fill_name($pcr1R[3],"4",$from-3);
   $pr2name = fill_name($pcr2F[3],"2",$from-1) if $f2 ne "\n";
   $pr3name = fill_name($pcr2R[3],"3",$from-2) if $r2 ne "\n";
 }
 else {
   $pr1name = fill_name($pcr1F[3],"F",$from);
   $pr4name = fill_name($pcr1R[3],"R",$from-3);
   $pr2name = fill_name($pcr2F[3],"F",$from-1) if $f2 ne "\n";
   $pr3name = fill_name($pcr2R[3],"R",$from-2) if $r2 ne "\n";
 }
 
 $posstr=join("",@pos);
 $namestr=join("",@name);
# print $name,"\n";
 print  $amp;
 print "\n$pcr1size";
 print join("\t","",$pr1name,$pcr1F[3],$pcr1F[4],$pcr1F[5],$pcr1F[6],$pcr1F[7]);
 print join("\t","",$pr4name,$pcr1R[3],$pcr1R[4],$pcr1R[5],$pcr1R[6],$pcr1R[7]);
 push @ORDER, "$pr1name\t$pcr1F[7]";
 push @ORDER, "$pr4name\t$pcr1R[7]";
 print "$pcr1par";
 if ($f2 ne "\n") {
  print "\n$pcr2size"; 
  print join("\t","",$pr2name,$pcr2F[3],$pcr2F[4],$pcr2F[5],$pcr2F[6],$pcr2F[7]);
  print join("\t","",$pr3name,$pcr2R[3],$pcr2R[4],$pcr2R[5],$pcr2R[6],$pcr2R[7]);
  $pcr2F[7] = "$tail2 $pcr2F[7]" if $tail2;
  $pcr2R[7] = "$tail3 $pcr2R[7]" if $tail3;
  push @ORDER, "$pr2name\t$pcr2F[7]";
  push @ORDER, "$pr3name\t$pcr2R[7]";
  print "$pcr2par";
 }
 print "\n\n";
 printIt($ampseq,$posstr,$namestr,100,$from-1);
 print "\n--------------\n\n";

}

print "\n\nPRIMERS:\n\n";
foreach $pr(@ORDER){ print $pr;}
print "\n\n";

foreach my $cols (@ORDER) {
    my $col = substr $cols, 0, -1 if $cols =~ /\n$/;
    $col = "FAILED" unless $col;
    print OUT $col . "\t";
#    print OUT $col;
}
print OUT "FAILED" unless @ORDER;
print OUT "\n";
close OUT;
##################################

sub fill_pos{
  $len = $_[0];
  $shift = $_[1];
  $fill = $_[2];
  for (0..$len){ $pos[$_+$shift]=$fill;}
}

sub fill_name{
 my $coord = $_[0];
 my $direction = $_[1];
 my $shift = $_[2];
 my $fill;
 if ($direction=~/[FR]/) { $fill = "-$coord$direction" } else { $fill = "-$direction"}
 my @prname=split(//,$prefix.$fill);
 my $F=$coord-$shift;
 my $T=$F + (scalar @prname);
 @name[$F..$T]=@prname;
 return join("",@prname);
}


sub printIt{
 my $s=$_[0];
 my $p=$_[1];
 my $n=$_[2];
 my $shift=$_[4];
 my $wrap=$_[3];
 for ($i=0; $i < length($s); $i=$i+$wrap){
  $subn = substr($n,$i,$wrap);
  $subp = substr($p,$i,$wrap);
  $subs = substr($s,$i,$wrap)." ".($i+$wrap+$shift);
  if ($subp=~/[><]/){ print "\n$subn\n$subp\n$subs\n"; }
  else { print "$subs\n"; }
 }
}
 




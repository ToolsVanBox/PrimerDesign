#!/usr/bin/perl -w

use strict;
use sigtrap;
use diagnostics;

my $mispriming=$ARGV[0];
my $pcr_type=$ARGV[1];
my $psr=$ARGV[2];
my $amp = "";
$amp = $ARGV[5] if $ARGV[5];

#warn "amplicon3.pl params: ",join("\t",@ARGV),"\n";


my @params;
 $params[0] = " -mingc 40 -maxpolyx 3 -gcclamp 1";
 $params[1] = " -mingc 30 -maxpolyx 3 -gcclamp 1";
 $params[2] = " -mingc 25 -maxpolyx 3 -gcclamp 1";
 $params[3] = " -mingc 20 -maxpolyx 3 -gcclamp 1";
 $params[4] = " -mingc 20 -maxpolyx 4 -gcclamp 1";
 $params[5] = " -mingc 20 -maxpolyx 3 -gcclamp 0";
 $params[6] = " -mingc 20 -maxpolyx 5 -gcclamp 1";
 $params[7] = " -mingc 20 -maxpolyx 5 -gcclamp 0";

{
     my $seq;
     my $header = <STDIN>;
     while (my $line=<STDIN>){ $seq .= $line; }
     #warn $header;
     #warn $seq;
     print $header;
     $header =~s/[>\n\r]//g;
     $header =~s/[^a-z,A-Z,0-9].*//;
     warn "\n#########\nWorking on $header\n#########\n";
#     $header=~s/ .*//g;

     if($seq){
        print "$seq\n\/\/\n";
        $seq=~s/[ 0-9\n\/\r]//g;
        my @amplicons;
	my $i=0;
	my $realseq=$seq;
	   $realseq=~s/[\[\]]//g;
	my $seqlength=length($realseq);
	my @seqarr=split(//,$realseq);
        while ($seq=~s/\[(.*?)\]/$1/){
	  my $from = length($`);
	  my $length = length($1);
	  $i++;
	  push @amplicons,([$from,$length,"$header\_$i.tempseq"]);
	  open (FA,">$header\_$i.tempseq") or die "cannot open outfile\n";
	  print FA ">$header amplicon $from-",$from+$length,"\n";
          my $to=$from+$length+500;
	  $to=$seqlength if $to >$seqlength;
	  $from=$from-500;
	  $from=0 if $from < 0;
	  for (my $j=0; $j < $from; $j++){ print FA "n";}
#	  print FA "\n\n";
	  for (my $j=$from; $j < $to; $j++){ print FA $seqarr[$j];}
#	  print FA "\n\n";
	  for (my $j=$to; $j < $seqlength; $j++){ print FA "n";}
	  print FA "\n";
	  close FA;
	}

        foreach my $part(@amplicons){
	        print "\nAMPLICON $amp ",$part->[0],"-",($part->[0]+$part->[1]),":\n";
		warn  "\nAMPLICON $amp ",$part->[0],"-",($part->[0]+$part->[1]),":\n";
		my $msg;
		$msg = "  picking primers for single PCR." if ($pcr_type eq "single");
		$msg = "  picking internal primers for nested PCR." if ($pcr_type eq "nested");
		$msg = "  picking internal primers for tilling PCR." if ($pcr_type eq "tilling");
	        my @innerprimers = find_primers($part->[0],$part->[1],$part->[2],$msg);
		if ($pcr_type eq "single") {
		    print @innerprimers[3,5,6,8,12],"\n\n\n\n\n\n\n";
		    next;
#		    $amp++ unless $amp eq "";
		}
		my @outerprimers;
		undef @outerprimers;
                if($innerprimers[5]=~/^NO/){
		  for (0..11){ push @outerprimers,"NO SUTABLE PRIMERS FOUND FOR PCR1\n";}
                   push @outerprimers, "NOT even searched because inner primers failed\n";
		}
		else{
		  my $pcr2length=$innerprimers[5];
		     $pcr2length=~s/.* //; $pcr2length=~s/\n//;
		  my @pcr2start=split(/ +/,$innerprimers[6]);
		     $innerprimers[5]=~s/^   1/PCR2/;
		     $innerprimers[3]=~s/^#/ /;
		  my $msg;
		     $msg = "  picking external primers for nested PCR." if ($pcr_type eq "nested");
		     $msg = "  picking external primers for tilling PCR." if ($pcr_type eq "tilling");
		     @outerprimers=find_primers($pcr2start[3],$pcr2length,$part->[2],$msg);
		}
   	        $outerprimers[5]=~s/^   1/PCR1/;
		$outerprimers[3]=~s/^#/ /;
	        print @outerprimers[3,5,6,8,12],"\n",
                      @innerprimers[3,5,6,8,12],"\n";
		$amp++ unless $amp eq "";

        }
     }
}



sub find_primers{
  my $from=$_[0] - 15;
  my $to=$from + $_[1] + 30;
  my $file=$_[2];
    warn "$file\n";
  my $msg = $_[3];
  my @primers;
  my $rangemin=$to-$from; my $rangemax=$rangemin+200;
#  my $psr = $rangemin . "-" . $rangemax;
 # warn "\n\nPSR 1 $psr\n\n";
#  $psr = "160-200";
  my $param;
 for (0..7){
  $param = $params[$_];
  warn "$msg Trying with $param\n";
  my $primer3_parameters="eprimer3 -sequence $file -outfile primer3.out -pickanyway yes ".
                         "-target $from,$to -prange $psr ".
                         "-numreturn 1 -opttm 57 -mintm 53 -maxtm 64 -psizeopt 0 -maxdifftm 10 -maxgc 80 ".
			 "-minsize 17 -maxsize 27 -mispriminglibraryfile $mispriming ".$param." 2>p3.err";
  open OUT, ">./output.txt";
 print OUT $primer3_parameters . "\n";
  close OUT;




  system("$primer3_parameters");
  open (PRIMERS,"<primer3.out") or die "cannot open outfile primer3.out:$!\n";
  local $/="\n";
  undef @primers;
  while (my $line=<PRIMERS>){ push @primers, $line;}
  close PRIMERS;
  push @primers, "Parameters: $param\n";
  if ($primers[5] ne "\n"){
   warn "\tOK.\n";
   return @primers ;
  }
  warn "\tFailed. Retrying with less stringent conditions.\n";
 }
 undef @primers;
 for (0..11){ push @primers,"NO SUTABLE PRIMERS FOUND FOR PCR2\n";}
 push @primers, "Parameters: $param\n";
 return @primers;
}

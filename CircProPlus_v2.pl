#!/usr/bin/perl

use warnings;
use strict;
use File::Path;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Perl;
my ($circ_in, $out, $ref_file, $gtf, $help, $species, $max_thread, $ribo_adaptor, $ribo_in, $rrna_file, $overlap_len, $replace, $replace_list, $dont_remove_temp);
Getopt::Long::GetOptions(
        'circ_in|C=s'        =>        \$circ_in,
        'ribo_in|R=s'        =>        \$ribo_in,
        'out|O=s'        =>        \$out,
        'ref_file|ref=s'        =>        \$ref_file,
        'rRNA_file|rRNA=s'        =>        \$rrna_file,
        'gene_anno|G=s'        =>        \$gtf,
        'adaptor|A=s'        =>        \$ribo_adaptor,
        'help|H!'        =>        \$help,
        'model|M=s'        =>        \$species,
        'thread_num|T=i'        =>        \$max_thread,
        'overlap_len|L=i'        =>        \$overlap_len,
        'rep|E=s'        =>        \$replace,
        'rep_list|L=s'        =>        \$replace_list,
        'dont_remove_temp|D!'        =>        \$dont_remove_temp,
);
if(!defined($circ_in) and !defined($help) and !defined($ref_file) and !defined($ribo_in) and !defined($species) and !defined($gtf)){
        print "Please use the --help or -H option to get usage information.\n";
}elsif(defined($help)){
        print '
Program:  circpro-plus
Modified:  liuyunchang@tmmu.edu.cn
Contact:  xianwen@zju.edu.cn

Usage:    perl CircPro.pl -c TotalRNASeq.fastq -o OutputDir -ref Genome.fa -g Gene.gtf -m ve -r RiboSeq.fastq -rRNA rRNA.fa -a Adaptor -t N

Arguments:
    -C, --circ_in
          FASTQ file from total/poly(A)- RNA-Seq. Paired-end FASTQ files should be separated by ",", e.g. "-C file_1.fastq,file_2.fastq". Multiple FASTQ files should be separated by ":", e.g. "-C file1.fastq:file2.fastq:file3.fastq".
    -R, --ribo_in
          FASTQ file from Ribo-Seq. Paired-end FASTQ files should be separated by ",", e.g. "-R file_1.fastq,file_2.fastq". Multiple FASTQ files should be separated by ":", e.g. "-R file1.fastq:file2.fastq:file3.fastq".
    -O, --out
          output dir
    -ref, --ref_file
          reference genome sequence in FASTA format
    -rRNA, --rRNA_file
          rRNA sequence in FASTA format (optional)
    -G, --gene_anno
          input GTF/GFF3 formatted annotation file name
    -A, --adaptor
          adaptor string (optional)
    -H, --help
          show help information
    -M, --model
          CNCI classification model ("ve" for vertebrate species, "pl" for plant species)
    -T, --thread_num
          number of threads for parallel running (default: 1)
    -L, --overlap_len
          minimal overlap length between Ribo-Seq reads and junction region (in each direction). (default: 5)
    -D, --dont_remove_temp
          Donnot remove temp files
    -E, --rep
          Replace the circRNA sequence with the CirCompara2 sequence (optional)
';
}else{
        print "CircProPlus V2 started!\nSequence mapping ...\n";
        if(!defined $dont_remove_temp){
          print "Will Removing temp files ...\n";
        }else{
          print "Will NOT Removing temp files ...\n";
        }
        if(defined $replace){
          print "Will Also calculate with REPLACE fa...\n";
          if(!defined $replace_list){
            print "Please provide the replace list file! Run convert.pl first!\n";
            exit;
          }
        }
        if(!defined $out){
          $out="./";
        }
        mkpath("$out/temp");
        if(!defined $max_thread){
          $max_thread=30;
        }
        if(!defined $overlap_len){
          $overlap_len=5;
        }
        my $ref_input_dir;
        if(rindex($ref_file, "/") >= 0){
          $ref_input_dir = substr($ref_file, 0, rindex($ref_file, "/")+1);
        }else{
          $ref_input_dir = "./";
        }
	# system "bwa index -a bwtsw $ref_file"; ## you're totally crazy when you try to create index EVERYTIME you're running exactly SAME reference!!!! 
        my $input_dir;
        if(rindex($circ_in, "/") >= 0){
          $input_dir = substr($circ_in, 0, rindex($circ_in, "/")+1);
        }else{
          $input_dir = "./";
        }
        if($circ_in =~ /,/){
          $circ_in =~ s/,/ /;
          system "bwa mem -v 1 -t $max_thread $ref_file $circ_in > $out/temp/bwa_result.sam";
        }elsif($circ_in =~ /:/){
          $circ_in =~ s/:/ /g;
          system "cat $circ_in > $out/temp/circ_in.fastq";
          system "bwa mem -v 1 -t $max_thread $ref_file $out/temp/circ_in.fastq > $out/temp/bwa_result.sam";
        }else{
          system "bwa mem -v 1 -t $max_thread $ref_file $circ_in > $out/temp/bwa_result.sam";
        }
        print "Done!\n";
        ## Done

        ## Detecting circRNAs with CIRI
        print "Detecting circRNAs ...";
        system "perl /user/public/circpro/CircPro/CIRI_v2.0.5.pl -I $out/temp/bwa_result.sam -O $out/temp/ciri_result.out -F $ref_file -low -T $max_thread -A $gtf -Q";
        system "echo $out/temp/ciri_result.out.log";       
        open(IN1,"<$out/temp/ciri_result.out") || die $! ;
        my $line_n=0;
        my $check=0;
        while(<IN1>){
          $line_n++;
          if($line_n>1){
            $check=1;
            last;
          }
        }
        if($check == 0){
          print "\nWarning: CIRI2 detected 0 circRNA!\nStoped!\n";
          system "cat $out/temp/CIRIerror.log";
          exit;
        }
        close(IN1);
        print "Removing bwa_result sam\n";
        system "rm -f $out/temp/bwa_result.sam";
        print "Done\n";
        
        ## Done
        ## circRNA Sequence
        print "Extracting circRNA sequence ...\n";
        open(IN2,"<$ref_file") || die $!;
        my $pos;
        my @chrs;
        while(<IN2>){
          if(/>(.*?)[ |\s]/){
            push(@chrs,$1);
            open(OUT1,">$out/temp/$1.fa") || die $!;
            print OUT1 ">$1\n";
            while(<IN2>){
              if(/^>/){
                seek(IN2,$pos,0);
                last;
              }
              else{
                print OUT1;
                $pos=tell(IN2);
              }
            }
            close(OUT1);
          }
        }
        close(IN2);
        my (%hash,$seqIO,$seq);
        foreach my $chr (@chrs){
          $seqIO = Bio::SeqIO->new( -file => "$out/temp/$chr.fa", -format => 'fasta');
          $seq = $seqIO->next_seq()->seq();
          $hash{"$chr"}=$seq;
        }
        open(IN3,"<$gtf") || die $!;
        open(CIRI,"<$out/temp/ciri_result.out") || die $! ;
        open(SEQ,">$out/temp/circRNA_temp.fa") || die $! ;
        readline(CIRI);
        my %circs;
        my %circ_list;
        my %circpos;
        my %circ_seq_final;
        my %othercirc;
        my %introniccirc;
        my %intergeniccirc;
        my $pos_t;
        while(<CIRI>){
          my @row=split/\t/;
          $circs{"$row[1]\t$row[2]\t$row[3]\t$row[10]"}=1;
          $circ_list{"$row[1]:$row[2]|$row[3]"}="$row[1]:$row[2]|$row[3]\t$row[1]\t$row[2]\t$row[3]\t$row[10]\t$row[4]";
        }
        my %circ_host;
        my %allintrons;
        while(<IN3>){
          if(/^(.*?)\t.*\ttranscript\t(\d+)\t(\d+)\t.\t(.).*gene_id \"(.*?)\"/ || /^(.*?)\t.*\tmRNA\t(\d+)\t(\d+)\t.\t(.).*Parent=(.*?);/){
            my ($chrG,$startG,$endG,$strandG,$G)=($1,$2,$3,$4,$5);
            my $pos_gtf=tell(IN3);
            foreach (keys %circs){
            	seek(IN3,$pos_gtf,0);
              my ($chrC,$startC,$endC,$strandC)=split/\t/;
              if($chrC eq $chrG && $strandC eq $strandG){
                next if ($startC > $endG || $endC < $startG);
                my %introns=();
                my $circ_id="$G\t$chrC\t$startC\t$endC\t$strandC";
                my $circseq="";
                if($strandG eq "+"){
                  my $intron_s=0;
                  while(<IN3>){
                    if(/\texon\t(\d+)\t(\d+)/){
                     if($intron_s == 0){
                      $intron_s=$2+1;
                     }
                     else{
                      my $intron_e=$1-1;
                      $introns{"$intron_s\t$intron_e"}=1;
                      $allintrons{"$chrG\t$intron_s\t$intron_e\t$strandG"}=$G;
                      $intron_s=$2+1;
                     }
                    }
                    if(/\tgene\t/ || /\ttranscript\t/ || /\tmRNA\t/){
                      last;
                    }
                  }
                  my %temp=();
                  $temp{"$startC\t$endC"}=1;
                  my $intron_check=1;
                  foreach (sort{$a cmp $b} keys %introns){
                    my @intron_info=split/\t/;
                    if($intron_info[0]<$startC && $intron_info[1]>$endC){
                      $intron_check=0;
                      last;
                    }
                    foreach (sort{$a cmp $b} keys %temp){
                      my @temp_k=split/\t/;
                      if($intron_info[0]>$temp_k[0] && $intron_info[1]<$temp_k[1]){
                        $temp{$temp_k[0]."\t".($intron_info[0]-1)}=1;
                        $temp{($intron_info[1]+1)."\t".$temp_k[1]}=1;
                        delete $temp{$temp_k[0]."\t".$temp_k[1]};
                        last;
                      }
                    }
                  }
                  next if ($intron_check==0);
                  my %temp1=();
                  foreach (keys %temp){
                    my @ks=split/\t/;
                    $temp1{$ks[0]}=$_;
                  }
                  foreach (sort{$a<=>$b} keys %temp1){
                    my @position=split/\t/,$temp1{$_};
                    $circ_id .= "\t$position[0]\t$position[1]";
                    my $temp_Seq = substr($hash{$chrC},$position[0]-1,$position[1]-$position[0]+1);
                    $circseq .= $temp_Seq;
                  }
                  if($circ_host{"$chrC\t$startC\t$endC\t$strandC"}){
                    if($circ_host{"$chrC\t$startC\t$endC\t$strandC"} ne $G){
                      $othercirc{"$chrC\t$startC\t$endC\t$strandC"}=1;
                      next;
                    }
                  }else{
                    $circ_host{"$chrC\t$startC\t$endC\t$strandC"}=$G;
                  }
                  $circ_seq_final{$circ_id}=$circseq;
                }else{
                  my $intron_e=0;
                  while(<IN3>){
                    if(/\texon\t(\d+)\t(\d+)/){
                      if($intron_e == 0){
                        $intron_e=$1-1;
                      }else{
                        my $intron_s=$2+1;
                        $introns{"$intron_s\t$intron_e"}=1;
                        $allintrons{"$chrG\t$intron_s\t$intron_e\t$strandG"}=$G;
                        $intron_e=$1-1;
                      }
                    }
                    if(/\tgene\t/ || /\ttranscript\t/ || /\tmRNA\t/){
                      last;
                    }
                  }
                  my %temp=();
                  $temp{"$startC\t$endC"}=1;
                  my $intron_check=1;
                  foreach (sort{$a cmp $b} keys %introns){
                    my @intron_info=split/\t/;
                    if($intron_info[0]<$startC && $intron_info[1]>$endC){
                      $intron_check=0;
                      last;
                    }
                    foreach (sort{$a cmp $b} keys %temp){
                      my @temp_k=split/\t/;
                      if($intron_info[0]>$temp_k[0] && $intron_info[1]<$temp_k[1]){
                        $temp{$temp_k[0]."\t".($intron_info[0]-1)}=1;
                        $temp{($intron_info[1]+1)."\t".$temp_k[1]}=1;
                        delete $temp{$temp_k[0]."\t".$temp_k[1]};
                        last;
                      }
                    }
                  }
                  next if ($intron_check==0);
                  my %temp1=();
                  foreach (keys %temp){
                    my @ks=split/\t/;
                    $temp1{$ks[0]}=$_;
                  }
                  foreach (sort{$b<=>$a} keys %temp1){
                    my @position=split/\t/,$temp1{$_};
                    $circ_id .= "\t$position[0]\t$position[1]";
                    my $temp_Seq = substr($hash{$chrC},$position[0]-1,$position[1]-$position[0]+1);
                    my $rev = Bio::Seq->new(-id => 'testseq', -seq => $temp_Seq);
                    $rev = revcom($rev);
                    $temp_Seq = $rev->seq();
                    $circseq .= $temp_Seq;
                  }
                  if($circ_host{"$chrC\t$startC\t$endC\t$strandC"}){
                    if($circ_host{"$chrC\t$startC\t$endC\t$strandC"} ne $G){
                      $othercirc{"$chrC\t$startC\t$endC\t$strandC"}=1;
                      next;
                    }
                  }else{
                    $circ_host{"$chrC\t$startC\t$endC\t$strandC"}=$G;
                  }
                  $circ_seq_final{$circ_id}=$circseq;
                }
              }
            }
            seek(IN3,$pos_gtf,0);
          }
        }
        foreach my $k_circ_seq_final (keys %circ_seq_final){
          my @row=split/\t/,$k_circ_seq_final;
          delete $circs{"$row[1]\t$row[2]\t$row[3]\t$row[4]"};
          if($othercirc{"$row[1]\t$row[2]\t$row[3]\t$row[4]"}){
            delete $circ_seq_final{$k_circ_seq_final};
          }
        }
        foreach my $k_circ_seq_final (keys %circ_seq_final){
          my @row=split/\t/,$k_circ_seq_final;
          $circ_list{"$row[1]:$row[2]|$row[3]"} .= "\texonic\t$row[0]";
          print SEQ ">$row[1]:$row[2]|$row[3]\n$circ_seq_final{$k_circ_seq_final}\n";##exonic
        }
        foreach (keys %othercirc){
          my @row=split/\t/;
          $circ_list{"$row[0]:$row[1]|$row[2]"} .= "\tother\tn/a";
          if($row[3] eq "+"){
            my $temp_Seq = substr($hash{$row[0]},$row[1]-1,$row[2]-$row[1]+1);
            print SEQ ">$row[0]:$row[1]|$row[2]\n$temp_Seq\n";##other
          }else{
            my $temp_Seq = substr($hash{$row[0]},$row[1]-1,$row[2]-$row[1]+1);
            my $rev = Bio::Seq->new(-id => 'testseq', -seq => $temp_Seq);
            $rev = revcom($rev);
            $temp_Seq = $rev->seq();
            print SEQ ">$row[0]:$row[1]|$row[2]\n$temp_Seq\n";##other
          }
        }
        foreach (keys %circs){
          my $check_intron=0;
          my @row_circs=split/\t/;
          foreach (keys %allintrons){
            my @row_intron=split/\t/;
            if($row_circs[0] eq $row_intron[0] && $row_circs[3] eq $row_intron[3] && $row_circs[1] >= $row_intron[1] && $row_circs[2] <= $row_intron[2]){
              $check_intron=1;
              $circ_list{"$row_circs[0]:$row_circs[1]|$row_circs[2]"} .= "\tintronic\t".$allintrons{"$row_intron[0]\t$row_intron[1]\t$row_intron[2]\t$row_intron[3]"};
              if($row_circs[3] eq "+"){
                my $temp_Seq = substr($hash{$row_circs[0]},$row_circs[1]-1,$row_circs[2]-$row_circs[1]+1);
                print SEQ ">$row_circs[0]:$row_circs[1]|$row_circs[2]\n$temp_Seq\n";##intronic
              }else{
                my $temp_Seq = substr($hash{$row_circs[0]},$row_circs[1]-1,$row_circs[2]-$row_circs[1]+1);
                my $rev = Bio::Seq->new(-id => 'testseq', -seq => $temp_Seq);
                $rev = revcom($rev);
                $temp_Seq = $rev->seq();
                print SEQ ">$row_circs[0]:$row_circs[1]|$row_circs[2]\n$temp_Seq\n";##intronic
              }
            }
          }
          if($check_intron==0){
            $circ_list{"$row_circs[0]:$row_circs[1]|$row_circs[2]"} .= "\tintergenic\tn/a";
            if($row_circs[3] eq "+"){
              my $temp_Seq = substr($hash{$row_circs[0]},$row_circs[1]-1,$row_circs[2]-$row_circs[1]+1);
              print SEQ ">$row_circs[0]:$row_circs[1]|$row_circs[2]\n$temp_Seq\n";##intergenic
            }else{
              my $temp_Seq = substr($hash{$row_circs[0]},$row_circs[1]-1,$row_circs[2]-$row_circs[1]+1);
              my $rev = Bio::Seq->new(-id => 'testseq', -seq => $temp_Seq);
              $rev = revcom($rev);
              $temp_Seq = $rev->seq();
              print SEQ ">$row_circs[0]:$row_circs[1]|$row_circs[2]\n$temp_Seq\n";##intergenic
            }
          }
        }
        open(CIRCLIST,">$out/temp/circlist.out") || die $!;
        foreach (keys %circ_list){
          print CIRCLIST $circ_list{$_}."\n";
        }
        close(CIRCLIST);
        close(CIRI);
        close(SEQ);
        close(IN3);
        open(SEQIN,"<$out/temp/circRNA_temp.fa") || die $!;
        open(SEQOUT,">$out/circRNA.fa") || die $!;
        open(CPCFILE,">$out/temp/cpcIN.fa") || die $!;
        my %circids;
        my %circseq_final;
        while(<SEQIN>){
          if(/^(>.*)/){
            chomp(my $circseq_temp=<SEQIN>);
            if($circids{$1}){
              $circids{$1}++;
              $circseq_final{$1."_".$circids{$1}}=$circseq_temp;
            }
            else{
                    $circids{$1}=1;
              $circseq_final{$1."_1"}=$circseq_temp;
            }
          }
        }
        foreach (sort{$a cmp $b} keys %circseq_final){
          print SEQOUT $_."\n".$circseq_final{$_}."\n";
          print CPCFILE $_."\n".$circseq_final{$_}.$circseq_final{$_}."\n";
        }
        close(SEQIN);
        close(SEQOUT);
        close(CPCFILE);



        print "Done!\n";
        ## Done

        ## lyc:using CPC2 instead of CPC
        print "Predicting circRNA protein coding potential with CPC2 ... \n";
        system "mkdir $out/temp/CPC.out/";
        system "CPC2.py -i $out/temp/cpcIN.fa -o $out/temp/CPC.out/CPC --ORF ";
        print "Done!\n";

#        ## lyc:using CPAT instead of CPC
#        print "Predicting circRNA protein coding potential with CPAT ...";
#        system "mkdir $out/temp/CPAT.out/";
#        system "CPC2.py -i $out/temp/cpcIN.fa -o $out/temp/CPC.out/CPC --ORF ";
#        print "Done!\n";

#        ## CPC&Ribo-seq
#        print "Predicting circRNA protein coding potential ...";
#        system "mkdir $out/temp/CPC.out/";
#        system "sh ./cpc-0.9-r2/bin/run_predict.sh $out/temp/cpcIN.fa $out/temp/CPC.out/CPC $out/temp/CPC.out/ $out/temp/CPC.out/CPC_evidence";
#        print "Done!\n";

        print "Extracting junction reads from Ribo-Seq data ...\n";
        #system "bowtie2-build -q $ref_file $ref_file --threads $max_thread";
        my ($ribo_in1,$ribo_in2,$seq_len);
        if($ribo_in =~ /(.*),(.*)/){
          $ribo_in1=$1;
          $ribo_in2=$2;
          open (RIBO,"<$ribo_in1") || die $! ;
          readline(RIBO);
          chomp(my $seq_tem=<RIBO>);
          $seq_len=length($seq_tem);
          close(RIBO);
          if(defined $ribo_adaptor){
            system "fastx_clipper -a $ribo_adaptor -c -i $ribo_in1 -o $out/temp/ribo_in1.filter";
            system "fastx_clipper -a $ribo_adaptor -c -i $ribo_in2 -o $out/temp/ribo_in2.filter";
            if(defined $rrna_file){
              system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
              system "bowtie2 -p $max_thread -x $rrna_file -1 $out/temp/ribo_in1.filter -2 $out/temp/ribo_in2.filter -S $out/temp/ribo_rrna.out.sam";
              system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
              system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }else{
              system "bowtie2 -p $max_thread -x $ref_file -1 $out/temp/ribo_in1.filter -2 $out/temp/ribo_in2.filter | samtools view -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }else{
                  if(defined $rrna_file){
              system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
              system "bowtie2 -p $max_thread -x $rrna_file -1 $ribo_in1 -2 $ribo_in2 -S $out/temp/ribo_rrna.out.sam";
              system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
              system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }else{
              system "bowtie2 -p $max_thread -x $ref_file -1 $ribo_in1 -2 $ribo_in2 | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }
        }elsif($ribo_in =~ /:/){
          my @ribo_files=split/:/,$ribo_in;
          open(RIBO,"<$ribo_files[0]") || die $! ;
          readline(RIBO);
          chomp(my $seq_tem=<RIBO>);
          $seq_len=length($seq_tem);
          close(RIBO);
          my $ribo_input_file="";
          if(defined $ribo_adaptor){
            my $i=0;
            foreach (@ribo_files){
                    $i++;
              system "fastx_clipper -a $ribo_adaptor -c -i $_ -o $out/temp/ribo_in$i.filter";
              $ribo_input_file .= "$out/temp/ribo_in$i.filter,";
            }
            if(defined $rrna_file){
              system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
              system "bowtie2 -p $max_thread -x $rrna_file -U $ribo_input_file -S $out/temp/ribo_rrna.out.sam";
              system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
              system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }else{
              system "bowtie2 -p $max_thread -x $ref_file -U $ribo_input_file | samtools view -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }else{
            if(defined $rrna_file){
              $ribo_in =~ s/:/,/g;
              system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
              system "bowtie2 -p $max_thread -x $rrna_file -U $ribo_in -S $out/temp/ribo_rrna.out.sam";
              system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
              system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }else{
              $ribo_in =~ s/;/ /g;
              system "bowtie2 -p $max_thread -x $ref_file $ribo_in | samtools view -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }
        }else{
          open(RIBO,"<$ribo_in") || die $! ;
          readline(RIBO);
          chomp(my $seq_tem=<RIBO>);
          $seq_len=length($seq_tem);
          close(RIBO);
          if(defined $ribo_adaptor){
            system "fastx_clipper -a $ribo_adaptor -c -i $ribo_in -o $out/temp/ribo_in.filter";
            if(defined $rrna_file){
              system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
              system "bowtie2 -p $max_thread -x $rrna_file $out/temp/ribo_in.filter -S $out/temp/ribo_rrna.out.sam";
              system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
              system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }else{
              system "bowtie2 -p $max_thread -x $ref_file $ribo_in | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }else{
                  if(defined $rrna_file){
                    system "bowtie2-build -q $rrna_file $rrna_file --threads $max_thread";
                    system "bowtie2 -p $max_thread -x $rrna_file $ribo_in -S $out/temp/ribo_rrna.out.sam";
                    system "samtools view -@ $max_thread -f 4 $out/temp/ribo_rrna.out.sam | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/riboseq.fa";
                    system "bowtie2 -f -p $max_thread -x $ref_file $out/temp/riboseq.fa | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
                  }else{
              system "bowtie2 -p $max_thread -x $ref_file $ribo_in | samtools view -@ $max_thread -f 4 | awk -f /user/public/circpro/CircPro_old/awk-f > $out/temp/ribo_unaln_temp.fa";
            }
          }
        }
        open(IN4,"$out/temp/ribo_unaln_temp.fa") || die $!;
        open(RIBOAD,">$out/temp/ribo_unaln.fa") || die $!;
        while(<IN4>){
          if(/^>/){
            chomp;
            print RIBOAD;
            my $ribo_unaln_seq=<IN4>;
            my $length_riboseq=length($ribo_unaln_seq)-1;
            print RIBOAD "_$length_riboseq\n$ribo_unaln_seq";
          }
        }
        close(IN4);
        close(RIBOAD);
        open(IN5,"$out/circRNA.fa") || die $!;
        open(RIBOREF,">$out/temp/ribo_ref.fa") || die $!;
        while(<IN5>){
          chomp;
          print RIBOREF;
          if(/^>(.*)/){
            my $circname=$1;
            chomp(my $circ_seq=<IN5>);
            my $circ_len=length($circ_seq);
            print RIBOREF "_".$circ_len."\n";
            my $riboref_seq;
            if($circ_len > (2*$seq_len)){
              $riboref_seq=substr($circ_seq,$circ_len-$seq_len,$seq_len).substr($circ_seq,0,$seq_len);
            }else{
              $riboref_seq=substr($circ_seq,$circ_len-int($circ_len/2),int($circ_len/2)).substr($circ_seq,0,int($circ_len/2));
            }
            print RIBOREF $riboref_seq."\n";
          }
        }
        close(IN5);
        close(RIBOREF);
        system "bowtie2-build -q $out/temp/ribo_ref.fa $out/temp/ribo_ref.fa --threads $max_thread";
        system "bowtie2 -p $max_thread -f -x $out/temp/ribo_ref.fa $out/temp/ribo_unaln.fa | samtools view -@ $max_thread -Sb | bamToBed -i > $out/temp/ribo_aln.out.bed";

        my %circ_result;
        my %circ_ribo_counts;
        open(CIRCRE,"<$out/temp/circlist.out") || die $!;
        while(<CIRCRE>){
          chomp;
          my @row=split/\t/;
          $circ_ribo_counts{$row[0]}=0;
          $circ_result{$row[0]}="$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]";
        }
        close(CIRCRE);

        open(IN6,"$out/temp/ribo_aln.out.bed") || die $!;
        while(<IN6>){
          next if(/\t-$/);
          my @bed_rows=split/\t/;
          my ($circ_name,$circ_iso,$circ_length)=split/_/,$bed_rows[0];
          my ($ribo_name,$ribo_length)=split/_/,$bed_rows[3];
          next if($circ_length < $ribo_length);
          if($circ_length > 2*$seq_len){
            if($bed_rows[1] <= $seq_len-$overlap_len && $bed_rows[2] >= $seq_len+$overlap_len){
              $circ_ribo_counts{$circ_name}=$circ_ribo_counts{$circ_name}+1;
            }
          }else{
            if($bed_rows[1] <= 1/2*$circ_length-$overlap_len && $bed_rows[2] >= 1/2*$circ_length+$overlap_len){
              $circ_ribo_counts{$circ_name}=$circ_ribo_counts{$circ_name}+1;
            }
          }
        }
        close(IN6);
        open(FINAL,">$out/CircPro.out") || die $!;
        print FINAL "circRNA\tchr\tstart\tend\tstrand\tjunction reads (RNA-Seq)\ttype\tparent gene\tjunction reads (Ribo-Seq)\n";
        foreach (sort{$a cmp $b}keys %circ_result){
          print FINAL $circ_result{$_}."\t".$circ_ribo_counts{$_}."\n";
        }
        close(FINAL);
        print "Done!\n";


        print "Extracting circRNA isoforms and potential protein sequences ...\n";
        my %circ_cpc_result;
        my %classification;
        open(IN7,"$out/temp/CPC.out/CPC.txt") || die $!;
        while(<IN7>){
          chomp();
          my @cpc_out=split/\t/;
          $circ_cpc_result{$cpc_out[0]}=$cpc_out[1]."\t".$cpc_out[2]."\t".$cpc_out[6]."\t".$cpc_out[5]."\t".$cpc_out[7]."\t".$cpc_out[8];
          $classification{$cpc_out[0]}=$cpc_out[8]; # in CPC2
          #circRNA 1transcript_length	2peptide_length	6orf_start 5ORF_integrity	coding_probability	label
          #circRNA classification CPC_score ORF_length
        }
        close(IN7);

        open(CIRCCPC,">$out/circIsoform.out") || die $!;
        #print CIRCCPC "circRNA\tclassification\tCPC score\tORF length\n";
        print CIRCCPC "circRNA\ttranscript_length\tpeptide_length\torf_start\torf_integrity\tcoding_prob\tclassification\n";
        foreach (sort{$a cmp $b} keys %circ_cpc_result){
          print CIRCCPC "$_\t$circ_cpc_result{$_}\n";
        }
        close(CIRCCPC);
        if (defined $replace){
          # copy from source
          system "echo 'Now Processing Replaced fa'";
          system "cp $replace $out/temp/circRNA_temp_rep1.fa";
          system "seqkit rmdup -s $out/temp/circRNA_temp_rep1.fa -w 0 -o $out/temp/circRNA_temp_rep.fa";
          open(SEQINR,"<$out/temp/circRNA_temp_rep.fa") || die $!;
          open(SEQOUTR,">$out/circRNA_rep.fa") || die $!;
          open(CPCFILER,">$out/temp/cpcIN_rep.fa") || die $!;
          my %circidsr;
          my %circseq_finalr;
          while(<SEQINR>){
            if(/^(>.*)/){
              chomp(my $circseq_tempr=<SEQINR>);
              if($circidsr{$1}){
                $circidsr{$1}++;
                $circseq_finalr{$1."_".$circidsr{$1}}=$circseq_tempr;
              }
              else{
                $circidsr{$1}=1;
                $circseq_finalr{$1."_1"}=$circseq_tempr;
              }
            }
          }
          foreach (sort{$a cmp $b} keys %circseq_finalr){
            print SEQOUTR $_."\n".$circseq_finalr{$_}."\n";
            print CPCFILER $_."\n".$circseq_finalr{$_}.$circseq_finalr{$_}."\n";
          }
          close(SEQINR);
          close(SEQOUTR);
          close(CPCFILER);

          # CPC identify circRNA protein coding potential
          print "Predicting REPLACED circRNA protein coding potential with CPC2 ... \n";
          system "mkdir $out/temp/CPC_rep.out/";
          system "CPC2.py -i $out/temp/cpcIN_rep.fa -o $out/temp/CPC_rep.out/CPC --ORF ";
          # build reference fa for ribo-seq
          open(IN5R,"$out/circRNA_rep.fa") || die $!;
          open(RIBOREFR,">$out/temp/ribo_ref_rep.fa") || die $!;
          while(<IN5R>){
            chomp;
            print RIBOREFR;
            if(/^>(.*)/){
              my $circname=$1;
              chomp(my $circ_seq=<IN5R>);
              my $circ_len=length($circ_seq);
              print RIBOREFR "_".$circ_len."\n";
              my $riboref_seq;
              if($circ_len > (2*$seq_len)){
                $riboref_seq=substr($circ_seq,$circ_len-$seq_len,$seq_len).substr($circ_seq,0,$seq_len);
              }else{
                $riboref_seq=substr($circ_seq,$circ_len-int($circ_len/2),int($circ_len/2)).substr($circ_seq,0,int($circ_len/2));
              }
              print RIBOREFR $riboref_seq."\n";
            }
          }
          close(IN5R);
          close(RIBOREFR);
          # build index for ribo-seq
          system "echo 'Running bash: bowtie2-build  $out/temp/ribo_ref_rep.fa $out/temp/ribo_ref_rep.fa --threads $max_thread'";
          system "bowtie2-build -q $out/temp/ribo_ref_rep.fa $out/temp/ribo_ref_rep.fa --threads $max_thread";
          system "echo 'Running bash: bowtie2 -p $max_thread -f -x $out/temp/ribo_ref_rep.fa $out/temp/ribo_unaln.fa | samtools view -@ $max_thread -Sb | bamToBed -i > $out/temp/ribo_aln_rep.out.bed'";
          system "bowtie2 -p $max_thread -f -x $out/temp/ribo_ref_rep.fa $out/temp/ribo_unaln.fa | samtools view -@ $max_thread -Sb | bamToBed -i > $out/temp/ribo_aln_rep.out.bed";
          system "echo 'Replace Process Done!'";

          my %circ_resultr;
          my %circ_ribo_countsr;
          open(CIRCRER,"<$replace_list") || die $!;
          while(<CIRCRER>){
            chomp;
            my @row=split/\t/;
            $circ_ribo_countsr{$row[0]}=0;
            $circ_resultr{$row[0]}="$row[0]\t$row[1]\t$row[2]\t$row[3]\t$row[4]\t$row[5]\t$row[6]\t$row[7]";
          }
          close(CIRCRER);

          open(IN6R,"$out/temp/ribo_aln_rep.out.bed") || die $!;
          while(<IN6R>){
            next if(/\t-$/);
            my @bed_rowsr=split/\t/;
            my ($circ_namer,$circ_isor,$circ_lengthr)=split/_/,$bed_rowsr[0];
            my ($ribo_namer,$ribo_lengthr)=split/_/,$bed_rowsr[3];
            next if($circ_lengthr < $ribo_lengthr);
            if($circ_lengthr > 2*$seq_len){
              if($bed_rowsr[1] <= $seq_len-$overlap_len && $bed_rowsr[2] >= $seq_len+$overlap_len){
                $circ_ribo_countsr{$ribo_namer}=$circ_ribo_countsr{$ribo_namer}+1;
              }
            }else{
              if($bed_rowsr[1] <= 1/2*$circ_lengthr-$overlap_len && $bed_rowsr[2] >= 1/2*$circ_lengthr+$overlap_len){
                $circ_ribo_countsr{$ribo_namer}=$circ_ribo_countsr{$ribo_namer}+1;
              }
            }
          }
          close(IN6R);
          open(FINALR,">$out/CircPro_rep.out") || die $!;
          print FINALR "circRNA\tchr\tstart\tend\tstrand\tjunction reads (RNA-Seq)\ttype\tparent gene\tjunction reads (Ribo-Seq)\n";
          foreach (sort{$a cmp $b}keys %circ_resultr){
            print FINALR $circ_resultr{$_}."\t".$circ_ribo_countsr{$_}."\n";
          }
          close(FINALR);

          print "Extracting REPLACE circRNA isoforms and potential protein sequences ...\n";
          my %circ_cpc_resultr;
          my %classificationr;
          open(IN7R,"$out/temp/CPC_rep.out/CPC.txt") || die $!;
          while(<IN7R>){
            chomp();
            my @cpc_outr=split/\t/;
            $circ_cpc_resultr{$cpc_outr[0]}=$cpc_outr[1]."\t".$cpc_outr[2]."\t".$cpc_outr[6]."\t".$cpc_outr[5]."\t".$cpc_outr[7]."\t".$cpc_outr[8];
            $classificationr{$cpc_outr[0]}=$cpc_outr[8]; # in CPC2
          }
          close(IN7R);

          open(CIRCCPCR,">$out/circIsoform_rep.out") || die $!;
          #print CIRCCPCR "circRNA\tclassification\tCPC score\tORF length\n";
          print CIRCCPCR "circRNA\ttranscript_length\tpeptide_length\torf_start\torf_integrity\tcoding_prob\tclassification\n";
          foreach (sort{$a cmp $b} keys %circ_cpc_resultr){
            print CIRCCPCR "$_\t$circ_cpc_resultr{$_}\n";
          }
          close(CIRCCPCR);
          print "Process with REPLACE output ALL Done!\n";
        }

        print "Additional CPAT Running ...\n";
        system "mkdir $out/CPAT_data";
        system "cpat -g $out/temp/cpcIN.fa -x /user/public/circpro/CPAT_dat/Mouse_Hexamer.tsv -d /user/public/circpro/CPAT_dat/Mouse_logitModel.RData --top-orf=5 -o $out/CPAT_data/plus";
        system "cpat -g $out/temp/cpcIN_rep.fa -x /user/public/circpro/CPAT_dat/Mouse_Hexamer.tsv -d /user/public/circpro/CPAT_dat/Mouse_logitModel.RData --top-orf=5 -o $out/CPAT_data/plus_rep";
        print "Done!\n";

        if(!defined $dont_remove_temp){
          system "rm -rf $out/temp/";
          print "Removing temp files ...\n";
        }
        ## Done
        print "Done!\n";
        print "CircProplus completed!\n";
}

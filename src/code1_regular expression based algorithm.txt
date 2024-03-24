unlink ("/root/Desktop/Generating_Pipeline/CODE/Final_Conserved.txt");
unlink ("/root/Desktop/Generating_Pipeline/CODE/Conserved.txt");
$path1="/root/Desktop/Generating_Pipeline/Forward/Human_CNE_hs1344.txt";
$path2="/root/Desktop/Generating_Pipeline/Forward/Mouse_CNE_hs1344.txt";
$enhancer_name1="";
$enhancer_name2="";
$input1="";
$input2="";
$aligned_seq1="";
$aligned_seq2="";
$A_aligned;
$A_Rev_aligned;
$result="";
@conserved="";
@final_conserv;
@forward;
@reverse;
$count=0;
mkdir("Forward");
mkdir ("Reverse");

#####################  Retreiving Human Sequence  ################################
$file1=open(data1,$path1) or die "Could not open the specified file"; 
while(<data1>){if ($_=~/>/){
chomp;
$enhancer_name1=$_;}
else 
{chomp;
$input1=$input1.$_;}
}
$enhancer_name1=substr($enhancer_name1,1);
$input1=~s/\s//g;
#print "$enhancer_name1 \nInput_1 : (".length($input1).")\n$input1\n\n";
close(data1);

####################  Retreiving Mouse Sequence  ###############################
$file2=open(data2,$path2) or die "Could not open the specified file"; 
while(<data2>){if ($_=~/>/){
chomp;
$enhancer_name2=$_;}
else 
{chomp;
$input2=$input2.$_;}
}
$enhancer_name2=substr($enhancer_name2,1);
$input2=~s/\s//g;
#print "$enhancer_name2 \nInput_2 : (".length($input2).")\n$input2\n\n";
close(data2);



#if ($enhance_name1 =~/[]//g){print "contains []";}

###################  GENERATING ERROR   #########################################
if (($enhancer_name1) eq  ($enhancer_name2) || ($enhancer_name1=~ /$enhancer_name2/g) ||($enhancer_name2=~ /$enhancer_name1/g))
{print "Both Sequences should be different.\n";
exit;} 

#################################################################################
print "Sequence 1 : $enhancer_name1 ...... ". length($input1)." ntds\n";
print "Sequence 2 : $enhancer_name2 ...... ". length($input2)." ntds\n";


                
                
                
                sub Main_function{
                $seq1="";
                $seq2="";
                $aligned_seq1=""; 
                $aligned_seq2="";
$path1=shift;
$path2=shift;
$folder=shift;
$enhancer_name1=shift;
$enhancer_name2=shift;
$sign=shift;
#print "$path1  $path2  $folder  $enhancer_name1  $enhancer_name2  $sign\n";
###################  Pair-wise Alignment  ######################################
`needle -asequence $path1 -bsequence $path2 -gapopen 10 -gapextend 0.5 -outfile /root/Desktop/Generating_Pipeline/CODE/$folder/output.needle `;
print "Running alignment ..... \n";
print "gapopen : 10 \ngapextend : 0.5\n";
print "Finished\n";

#################  Reading output.needle file  ################################
$enhan_name1=substr($enhancer_name1,0,13);
$enhan_name2=substr($enhancer_name2,0,13);
open (data, "/root/Desktop/Generating_Pipeline/CODE/$folder/output.needle") or die("Could not open  file.");
while (<data>){
if ($_=~/$enhan_name1/){
$seq1=$seq1.$_;
 }
if ($_=~/$enhan_name2/){
$seq2=$seq2.$_;}
}
close (data);

################  Retreiving sequences in aligned form #######################
open (data,">", "/root/Desktop/Generating_Pipeline/CODE/$folder/align_seq1.txt") or die("Could not open  file.");
print data $seq1;
close(data);
open (data, "/root/Desktop/Generating_Pipeline/CODE/$folder/align_seq1.txt") or die("Could not open  file.");
while (<data>){if ($_=~/#/g){print "";}else {
$aligned_seq1=$aligned_seq1.$_;}}
#print "$aligned_seq1\n";
close (data);
open (data,">", "/root/Desktop/Generating_Pipeline/CODE/$folder/align_seq2.txt") or die("Could not open  file.");
print data $seq2;
close(data);
open (data, "/root/Desktop/Generating_Pipeline/CODE/$folder/align_seq2.txt") or die("Could not open  file.");
while (<data>){if ($_=~/#/g){print "";}else {
$aligned_seq2=$aligned_seq2.$_;}}
#print "$aligned_seq2\n";
close (data);
$aligned_seq1=~s/$enhan_name1//g;
$aligned_seq2=~s/$enhan_name2//g;
$aligned_seq1=~s/[0-9]+//g;
$aligned_seq2=~s/[0-9]+//g;
$aligned_seq1=~s/\s+//g;
$aligned_seq2=~s/\s+//g;
#print "\nAligned_seq1 : ".length($aligned_seq1)."\n$aligned_seq1\n";
#print "\nAligned_seq2 : ".length($aligned_seq2)."\n$aligned_seq2\n";
#unlink("/root/Desktop/Generating_Pipeline/CODE/align_seq1.txt");
#unlink("/root/Desktop/Generating_Pipeline/CODE/align_seq2.txt");
#foreach $item (keys %transcode){ $patlen=length $transcode{$item};
  #  print "$item : *$transcode{$item}  :  *$mod_trans{$item}*  : $patlen\n ";}
   searching($aligned_seq1);
 writing_result("$enhancer_name1",$folder);
 searching($aligned_seq2);
 writing_result("$enhancer_name2",$folder); 
 return $aligned_seq1;}

####################  Finding TFs patterns in aligned sequences ##################
#foreach $item (keys %transcode){ $patlen=length $transcode{$item};
  #  print "$item : *$transcode{$item}  :  *$mod_trans{$item}*  : $patlen\n ";}
  
   sub searching{
%ambcode=('A','A','C','C','G','G','T','T','R','[AG]','Y','[CT]','M','[AC]','K','[GT]','S','[CG]','W','[AT]','H','[ACT]','B','[CGT]','V','[ACG]','D','[AGT]','N','[ACGT]');
$factor_list='/root/Desktop/Forebrain/Forebrain_32Enhancers_sequences/TRANSFAC_Codes.txt';
open (handler,$factor_list) or die("Could not open  file.");
while (<handler>){if ($_=~TF){print "\n";}
else
	{@sp=split(":",$_);
	$sp[0]=~ s/\s//g;
	@key[$i]=$sp[0];
	$sp[1]=~ s/\s//g;
	@val[$i]=$sp[1];
	$transcode{$key[$i]}=$val[$i];	
	#print "$key[$i]  :   $val[$i]\n"; 
        $i=$i+1;
        }
        }
$size=keys %transcode;
#print "Number of Transcription factors to be scanned : $size\n";
%mod_trans=%transcode;
foreach $item(values %mod_trans){
foreach $num(keys %ambcode){
$item=~ s/$num/@ambcode{$num}/g;}
}


                          
my $inp= shift;
my $a=0;
my $i=0;
my $seq="";
$result="";
$len=length($inp);
#print "$len\n";
#print "\nNo.\tTF\t\t\tPosotion\tSequence\t\t\tPattern:  \n\n";
$result="\nNo.\tTF\t\t\tPosotion\tSequence\t\t\tPattern\t\t\tStarnd:  \n\n";
	
		foreach $item (keys %transcode){
 #   print "Item : $item *$transcode{$item}  :  *$mod_trans{$item}*\n";
    $temp=$item; @tem=split("_",$temp);
    $temp=substr($tem[0],2);
 #   print "$temp\n";
    $patlen=length $transcode{$item};
#print "^^$patlen^^\n";
for ($i=0;$i<=$len-$patlen;$i=$i+1){ #print "^^$patlen^^\n";
	$in=substr($inp,$i,$patlen);
	if ($in=~ /$mod_trans{$item}/ig)
	{  
	    $orgnl= ($len-($i+1));
	    $orgnl_pos=$len-$orgnl;
	    $pos[$a]=$orgnl_pos;
		$seq=substr($inp,($pos[$a]-1),$patlen);
	$ser_num=$a+1;
#	print $a+1,".\t",$item,"\t\t",$pos[$a],"\t\t$seq\t\t$mod_trans{$item}\t\n";	
	$result=$result."$ser_num.\t$item\t\t$pos[$a]\t\t$seq\t\t$mod_trans{$item}\t\t\t$sign\n";
$a=$a+1;	}}
}}  

                             sub writing_result{
$file_name=shift;
#unlink("C:/Users/Nus/Desktop/TFBStool_Results/".$enhancer_name."."."txt");
unlink("/root/Desktop/Generating_Pipeline/CODE/$folder/".$file_name."."."txt");
#$output="C:/Users/Nus/Desktop/TFBStool_Results/".$enhancer_name."."."txt";
$output="/root/Desktop/Generating_Pipeline/CODE/$folder/".$file_name."."."txt";
open(data,">",$output);
print data "Table View :\n";
print data $result,"\n\n";
close (data);
}

$A_aligned=Main_function($path1,$path2,"Forward",$enhancer_name1,$enhancer_name2,"+");

################# COMPLEMENT ############################################################
$path3="/root/Desktop/Generating_Pipeline/CODE/A_Rev.txt";
$path4="/root/Desktop/Generating_Pipeline/CODE/B_Rev.txt";

sub complement{
$input=shift;
$comp="";
for($i=0;$i<length $input;$i=$i+1){
	if (substr($input,$i,1) =~/[tT]/){
		$comp=$comp."A";}
		if (substr($input,$i,1)=~/[aA]/){
		$comp=$comp."T";
		}
		if (substr($input,$i,1) =~/[gG]/){
		$comp=$comp."C";
		}
		if (substr($input,$i,1) =~/[cC]/){
		$comp=$comp."G";
		}
	}
	$comp=reverse($comp); return $comp;}
	$A_Rev= complement("$input1");
	open (data,'>',$path3);
        print data ">$enhancer_name1\n";
	print "\n\nComplement (",length($A_Rev),") :\n"; 
	print data $A_Rev;
	close (data);
	$B_Rev= complement("$input2");
	open (data,'>',$path4);
        print data ">$enhancer_name2\n";
	print "\n\nComplement (",length($B_Rev),") :\n"; 
	print data $B_Rev;
	close (data);
	
$A_Rev_aligned=Main_function($path3,$path4,"Reverse",$enhancer_name1,$enhancer_name2,"-");
#print "\nA_aligned :\n$A_aligned\n\nA_Rev_aligned :\n$A_Rev_aligned\n";

####################### Finding Conserved patterns in Human and Mouse ####################
#########################################################################################
                    sub Conservation{
$folder=shift;
$sign=shift;
@TF_name_hum="";
@TF_pos_hum="";
@TF_seq_hum="";
@TF_pat_hum="";
$counter=1;

## EXTRACTING HUMAN TF NAME, POS, SEQ AND PATTERN FROM HUMAN FILE #######
open(data,"/root/Desktop/Generating_Pipeline/CODE/$folder/".$enhancer_name1."."."txt");
while (<data>){
#print "$_\n";
chomp;
if ($_=~/$counter/){
@sp=(split /\t/,$_)[0..7];
#print"0*$sp[0]   1*$sp[1]   2*$sp[2]   3*$sp[3]   4*$sp[4]   5*$sp[5]    6*$sp[6]    7*$sp[7]\n";
$TF_name_hum[$counter-1]=$sp[1];
$TF_pos_hum[$counter-1]=$sp[3];
$TF_seq_hum[$counter-1]=$sp[5];
$TF_pat_hum[$counter-1]=$sp[7];

$counter=$counter+1;

}
}
close (data);
$hum_length=@TF_name_hum;
#for ($i=0;$i<$hum_length;$i=$i+1) {print  "$TF_name_hum[$i]    $TF_pos_hum[$i]    $TF_seq_hum[$i]    $TF_pat_hum[$i]\n";    }

### EXTRACTING MOUSE TF NAME, POS, SEQ AND PATTERN FROM MOUSE FILE #####
open(data,"/root/Desktop/Generating_Pipeline/CODE/$folder/".$enhancer_name2."."."txt");
@TF_name_mou="";
@TF_pos_mou="";
$counter=1;
while (<data>){
chomp;
#print "$_\n";
if ($_=~/$counter/){
@sp=(split /\t/,$_)[0..7];
#print"0*$sp[0]   1*$sp[1]   2*$sp[2]   3*$sp[3]   4*$sp[4]   5*$sp[5]    6*$sp[6]    7*$sp[7]\n";
$TF_name_mou[$counter-1]=$sp[1];
$TF_pos_mou[$counter-1]=$sp[3];
$counter=$counter+1;

}
}
close (data);

$mou_length=@TF_name_mou;
#for ($i=0;$i<$mou_length;$i=$i+1) {print  "$TF_name_mou[$i]    $TF_pos_mou[$i] \n";    }


## FINDING CONSERVED POSTION BETWEEN HUMAN AND MOUSE ####
@conserv_name="";
@conserv_pos="";
@conserv_seq="";
@conserv_pat="";
$counter=0;


for ($i=0;$i<$hum_length;$i=$i+1){
for ($j=0;$j<$mou_length;$j=$j+1){
if ($TF_name_hum[$i] eq $TF_name_mou[$j] && $TF_pos_hum[$i] == $TF_pos_mou[$j] ){
#print "\n*$TF_name_hum[$i]  $TF_pos_hum[$i]  $TF_name_mou[$j]   $TF_pos_mou[$j]  \n";
$conserv_name[$counter]=$TF_name_hum[$i];
$conserv_pos[$counter]=$TF_pos_hum[$i];
$conserv_seq[$counter]=$TF_seq_hum[$i];
$conserv_pat[$counter]=$TF_pat_hum[$i];
$counter=$counter+1;}
else {}
}
}
$counter=1;
$conserv_length=@conserv_name;
if ($conserv_length <=1){print "No Conserved pattern is found.\n";}
else {
#print "\nNo.\tTF\t\tPosotion\tSequence(Strand)  \n\n";
for ($i=0;$i<$conserv_length;$i=$i+1){ #print "$counter. \t$conserv_pos[$i]\t$conserv_name[$i]\t$conserv_seq[$i]($sign)\n";
#print data "$counter. \t$conserv_name[$i]\t$conserv_pos[$i]\t$conserv_seq[$i]\t$sign\n";
$conserved[$count]="$conserv_pos[$i]\t$conserv_name[$i]\t$conserv_seq[$i]\t$sign\n";
$count=$count+1;
$counter=$counter+1;
}
}
}
   Conservation("Forward","+");
   Conservation("Reverse","-");
   @conserved = sort {$a <=> $b} @conserved;
   $size_lines= @conserved; print "############## $size_lines ###########\n";
   $num=1;
   open (data,">>", "/root/Desktop/Generating_Pipeline/CODE/Conserved.txt");
 #  print data "No.\tPosition\tTF\tSequence\tStrand  \n\n";
 for ($i=0;$i<$size_lines;$i=$i+1){if ($size_lines<=1){print data "No Conserved pattern is found \n"; exit;}else {print data "$num.\t$conserved[$i]\n";$num=$num+1;}}
 close (data);

 counting_for_rev_patterns("Forward");
 counting_for_rev_patterns("Reverse");
 #################################### Counting forward & reverse patterns #########################
 sub counting_for_rev_patterns{
 $strand=shift;
 $f=0;$r=0;
  open (data, "/root/Desktop/Generating_Pipeline/CODE/Conserved.txt") or die "Could not open Conserved.txt";
  while (<data>){
  chomp;
  if ($strand eq "Forward"){
 if ($_=~/\+/){ $forward[$f]=$_; #print "$f   $forward[$f]\n";
 $f=$f+1;}
 }
 } close (data);
   open (data, "/root/Desktop/Generating_Pipeline/CODE/Conserved.txt") or die "Could not open Conserved.txt";
  while (<data>){
  chomp;
  if ($strand eq "Reverse"){
 if ($_=~/\-/){ $reverse[$r]=$_; #print "$r   $reverse[$r]\n";
 $r=$r+1;}
 }
 } close (data);
 }
 #$for_len=@forward;print "for_len $for_len\n";
 #$rev_len=@reverse;print "rev_len $rev_len\n";
$count=0;
normalized($A_aligned,@forward);
normalized($A_Rev_aligned,@reverse);

###################################### Normalized Positions ##############################
sub normalized{
$input=shift;
@pos_count=@_;
$len=length($input);
$l=@pos_count;# print "l : $l\nlength : $len\n";
$gap=0;
for ($i=0;$i<$len;$i=$i+1){
$inp=substr($input,$i,1);
if ($inp eq '-'){$gap=$gap+1;}
else { for ($j=0;$j<$l;$j=$j+1){
@sp=(split/\t/,$pos_count[$j])[0..4]; #print "@sp\n";
$patlen=length($sp[3]);
$inp=substr($input,$i,$patlen);
if ($inp eq $sp[3] && $sp[1] == ($i+1)  )
{ #print "^^^^^^^^^^^^found \n";
$final_pos=($i-$gap)+1;
$final_conserv[$count]="$final_pos\t$sp[2]\t$sp[3]\t$sp[4]\n";
$count=$count+1;}

}
}
}
}

#@final_conserv=sort {$a<=>$b}@final_conserv;
$final_len=@final_conserv;
#print "final len : $final_len\n@final_conserv\n";
$size_lines=@final_conserv;
#print $size_lines;
open (data,">","/root/Desktop/Generating_Pipeline/CODE/Final_Conserved.txt");
for ($i=0;$i<$size_lines;$i=$i+1){
$num=$i+1;print data "$num\t$final_conserv[$i]"; }



################################ deleting ################################################
=begin
system('rm','-r',"/root/Desktop/Generating_Pipeline/CODE/Forward");
system('rm','-r',"/root/Desktop/Generating_Pipeline/CODE/Reverse");
system('rm','-r',$path3);
system('rm','-r',$path4);
=cut

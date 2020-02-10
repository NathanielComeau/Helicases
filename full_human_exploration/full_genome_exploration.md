## Full Human Genome Exploration

- Started with ERR013138_1.fastq obtained from https://www.internationalgenome.org/data-portal/sample/HG00114
- Tried to extract all identifiers and qa data from it with awk, still wasn't finished three minutes later...
- Killed it after ~two minutes, file was 1.2 GB. Implies total file will take ~5 minutes to extract.
```bash
(bio) w134-87-145-178:full_human_exploration nat$ ls -lh
total 18112232
-rw-r--r--@ 1 nat  staff   7.5G  9 Feb 12:03 ERR013138_1.fastq
-rw-r--r--@ 1 nat  staff   254B  9 Feb 15:22 full_genome_exploration.md
-rw-r--r--  1 nat  staff   1.2G  9 Feb 15:22 quality_scores.txt
```
- Extracted quality scores and identifier fields to text files with qa_manipulation.py and id_manipulation.py

- Looked at quality score data lined up with flow cell lane, flow cell number, x coord, and y coord.

```bash
paste quality_scores.txt identifierComponents/flow_cell_lane.txt identifierComponents/flow_cell_number.txt identifierComponents/x_coord.txt identifierComponents/y_coord.txt > side_by_side.txt
```

```txt
Quality Score                                                                                                 Lane Num  X-Cor.  Y-Cor.

GGAGGE;AEDE=7ADE;>CE;97;!3!!:!3354!4351(67855:C:9755468/81-3!5642433/63;((&20!(12!2-!!!!2!!!-/4&-(%$+!!!!!(!	5	1	1035	1408
A@5@>A;4:<<>;4;2><.93134!4!!-!110(!25..3:480642145/'1.31.)14!1%'('//1%'%5,(&&!()0!-)!!!!/!!!+)*$'&')%!!!!!&!	5	1	1035	8133
B9>A@@4@@>@<<@66@@@>-31,!4!!2!3424!12./3=>;><>D@@@>@3>178<51!:-3-44+/208'4985!453!2&!!!!'!!!56-/.,+1(!!!!!3!	5	1	1035	6544
@5;43==768=3:;=<899:0326!3!!1!/1,1!01*5)24434-747931422,0./(!1,0%132/16421&51!20-!.'!!!!'!!!'-,).%/)+!!!!!&!	5	1	1035	2520
GG@GGGGGGGAEC@A9FFBD>?>?!1!!4!2753!7;6<=DGCG9EFGAGGFEG7?>7<>!8628151356B;6722!449!32!!!!(!!!115-03323!!!!!2!	5	1	1036	15336
AA>A;7<==686@:=BBB>B57;5!3!!2!/134!7576798998;7638(322111093!-10+5&1-124.2501!20&!(1!!!!*!!!*+0-(+'6(!!!!!1!	5	1	1036	18571
A<GGGGGGGGGG@GG>DGEA9<:<!4!!2!3214!999;9BEAAG>9B@>>9GEG<<999!964><55631;?>9/5!444!/3!!!!2!!!+.5-2123+!!!!!)!	5	1	1036	7118
GGGGGGGGGGBG9=EFF>DF<A?=!7!!3!95-:!<:>?<EA?@>GGAGGGGGGA<;<>>!<9>7557470<;:960!974!25!!!!1!!!03.2.*3.3!!!!!3!	5	1	1036	12742
4AA<C96;DAA=;:AA9@;@;977!9!!0!'554!1447)774;)18570<<95044591!431*4443/7756)14!.+3!1(!!!!)!!!.'2-1.(24!!!!!5!	5	1	1036	11760
CE>65EC<9AC>@=@GGGAG9=>9!1!!5!'432!45472@@@;8>6ABADDD=D77778!83747(044;>;8642!305!+4!!!!$!!/343125034!!!1!+!	5	1	1036	5816

```

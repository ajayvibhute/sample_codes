#!/bin/csh

set path1 = "../log"
set path2 = "../shadows/preprocess"

set infile = "jobQueue"

cd $path2
#echo $path2

set cnt = `ls -d */ | wc -l` 
#echo "count =  $cnt"

set i = 1
set j = 1



while($i <= $cnt)
	set dir  = `ls -d * | head -$i | tail -1`
#	echo "dir = $dir"
	set rcnt = `wc -l $path1/$infile | awk '{print $1}'`

if($rcnt == 0) then 
echo $dir >> $path1/$infile
mkdir 	$dir/"output"	
#echo $infile
else
goto loop
endif

loop:
while($j <= $rcnt)
	set record = `awk -v var=$j '{if(NR==var) {print $1}}' $path1/$infile`
	if ("$dir" == "$record") then
	break
	else
	if ($j == $rcnt) then
	echo $dir >> $path1/$infile
	
	mkdir 	$dir/"output"
#	echo "$infile"
	else
	endif
	endif
	
     @ j++
   end

@ i++
end

cd -
#exit -1

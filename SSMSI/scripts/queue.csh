#!/bin/csh


set param = $argv[1]
set file = "inputfile"


if ($param == "add")  then
	set queue = $argv[2]
	echo $queue >>$file
endif


if($param == "remove") then
		set line = 1
		echo $line | sed -i "$line d" $file
endif
endif

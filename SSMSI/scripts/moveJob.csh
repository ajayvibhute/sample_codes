#!/bin/csh 
	set jobPath = "../shadows/log"
	set currentJobPath = "../shadows/preprocess"
	set jobQueue = "jobQueue"
	set currentJobs = "currentJobs"
	set processedJobs = "processedJobs"
	set dest = $argv[1]
	set record = $argv[2]
	
	if($dest == "process") then
		sed -i "/$record/d" $jobPath/$jobQueue
		echo $record >>$jobPath/$currentJobs
		#chmod 400 $currentJobPath/$record
		mv "../shadows/preprocess/"$record "../shadows/underprocess/"	
			
	endif
	
	if($dest == "done") then
		sed -i "/$record/d" $jobPath/$currentJobs
		echo $record >>$jobPath/$processedJobs
		#chmod 755 $currentJobPath/$record
	endif
	if($dest == "repository") then
		mv "../shadows/underprocess/"$record "../shadows/postprocess/"
	endif


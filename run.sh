#! /bin/bash
# Run on a quick test image

FILE="img.c"
BASE={FILE%.*}

if [ -e $BASE ]
then
	rm $BASE
fi

# Compile my program
gcc img.c -o img -O3 -g -Wall -lpng12 -lm

if [ -e img ]
# Successful compile
then
    DIR=$(pwd)

    #valgrind --leak-check=full /home/brandon/Documents/meta/img vid_frames/out_0001.png processed_vid_frames/2.png
    
    if [ -e processed_vid_frames/2.png ]
    then
        eog processed_vid_frames/2.png &
    fi

    ./plot.py $DIR/processed_vid_frames/$fName

    #avconv -r 30 -i "to_combine_2/out_out_%04d.png" -y $OUT
	
	# Send me an email to notify of the processes' completion
    mail -s "Process Finished" "bjd2385@aperiodicity.com" < /dev/null
fi

exit 0

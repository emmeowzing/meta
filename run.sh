#! /bin/bash
# Run on a quick test image

FILE="img.c"
BASE={FILE%.*}
#IMAGE="example.png"
IMAGE="out_out_1.png"


if [ -e $BASE ]
then
    rm $BASE
fi

# Compile my program
gcc $FILE -o $BASE -O3 -g -Wall -lpng12 -lm

if [ -e img ]
# Successful compile
then
    DIR=$(pwd)

    # Code has successfully passed valgrind with no memory leaks so far
    #valgrind --leak-check=full /home/brandon/Documents/meta/img vid_frames/out_0001.png processed_vid_frames/2.png
    
    # Show the frame we're working with
    if [ -e $IMAGE ]
    then
        eog $IMAGE &
    fi
    
    # Run on the example frame
    ./$BASE $IMAGE "out_$IMAGE"

    # Run plot.py on the frame (wherein matplotlib automatically scales the )
    ./plot.py "out_$IMAGE"
    
    eog "fig_out_$IMAGE"
    
    # So there are two images created here:
    # "out_$IMAGE" is where the processed data are written
    # "fig_out_$IMAGE" is created by plot.py

    # For looping these steps, etc.
    #avconv -r 30 -i "to_combine_2/out_out_%04d.png" -y $OUT
    
    # Send me an email to notify of the processes' completion
    mail -s "Process Finished" "bjd2385@aperiodicity.com" < /dev/null
fi

exit 0

#!/bin/sh

echo
echo "SHELL SCRIPT FOR MOVIE USING FFMPEG"
echo



# WRITE YOUR FILE NAME IN BETWEEN ""
imageDIR="../data/animated_images"
movieDIR="../data"
fps=20

################## SEARCH AND STORE ########################

ffmpeg -r 50 -i "$imageDIR/void_%04d.png" -c:v libx264 -vf fps=$fps -pix_fmt yuv420p "$movieDIR/animate.mp4"

echo "DONE!"

Program requirements:
=====================
*python3

*numpy

*astropy

*scipy

*matplotlib

Usage
=====
python analyze.py \[path to .fits image\] \[standard deviations above mean to count as star\]

Start with 1 standard deviation, and if the result is sloppy, gradually increase it until the result is acceptable.

Notes
=====
The star detection algorithm is very primitive: it simply looks for pixels a given number of standard deviations brighter than the mean. As a result, the program can't handle multiple stars in one image, so you'll have to crop the image down to one star before you run it through the program. If I wanted to make it better I could use some sort of Gaussian fit, but this works for my purposes.

Also, the program is really slow. It takes a good ten seconds or so to run on a 640x480 image. I didn't bother to optimize it all; as much as I love making code run fast, I'll only run this program a few times, so the amount of time I'll save is much less than the time it would take to optimize and debug it.

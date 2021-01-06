 !#usr/bin/env bash
 
 # setting the path of the parent directory 
 path="$1"
   
 # finding all the .csv files under path  
 for i in `find $path -type f -name '*.csv'` 
 do
    # sorting based on 4th columns(k4) by using comma(,) as the field separator
    sort -t"," -k1 $i 
 done

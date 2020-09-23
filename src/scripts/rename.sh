#rename the files that start with first arg (e.g LUAD) with second arg (e.g. BRCA)
for n in $1*.sh 
do 
    cp $n $(echo $n | sed -e "s/$1/$2/")
done

#iterate over the newly renamed files 
# 

for n in $2*.sh
do
   sed -i "s/$1/$2/I" $n
done


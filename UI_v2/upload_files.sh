

list=$(awk '{print $3}' ../../name_species)
for query in ${list}; do
echo "${query}"

count=1
find=true

while ${find} = true:
do
find=false
imagelink=$(wget --user-agent 'Mozilla/5.0' -qO - "www.google.be/search?q=${query}\&tbm=isch&tbs=il:cl" | sed 's/</\n</g' | grep '<img' | head -n"$count" | tail -n1 | sed 's/.*src="\([^"]*\)".*/\1/')

if wget $imagelink -O ${query}.png ; then
    echo "Command succeeded"
    wget $imagelink -O ${query}.png
else
    echo "Command failed"
    count=2
    find=true
fi

done
done


### Alternative method

list=$(cat ../name_sp.txt)
for query in ${list}; do
echo "${query}"

query="drosophila_melanogaster"
googliser -p ${query} --format "png" -n 1 -m 2mp
cp ./${query}/* ./${query}.jpeg
convert ${query}.jpeg ${query}.png
rm ${query}.jpeg
rm -r ${query}
done

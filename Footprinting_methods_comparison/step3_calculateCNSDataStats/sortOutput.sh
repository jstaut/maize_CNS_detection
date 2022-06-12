
# USAGE: ./sortOutput.sh <output folder>

cd $1
mkdir sizeInfo
find * -maxdepth 0 -type f -name '*_CNSSizes.txt' -exec sh -c 'mv "$0" sizeInfo' {} \;
find * -maxdepth 0 -type f -name '*.txt' -exec sh -c 'sort -b -o "$0" "$0"' {} \;

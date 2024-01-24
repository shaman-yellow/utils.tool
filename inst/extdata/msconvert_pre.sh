path=`readlink -f $1`
docker run -it --rm \
  -e WINEDEBUG=-all \
  -v ${path}:/data \
  chambm/pwiz-skyline-i-agree-to-the-vendor-licenses \
  wine msconvert $2 --mzML --verbose --64 \
  --filter "peakPicking vendor" \
  --filter "titleMaker <RunId>.<ScanNumber>.<ScanStartTimeInMinutes>.<ChargeState>" \
  --filter "polarity $4" \
  -o $3
## bash msconvert_pre.sh /home/echo/operation file.raw targetfile
# sudo chown -R $3

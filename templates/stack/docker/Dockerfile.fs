FROM alpine

RUN apk update && apk add wget tar
RUN wget https://github.com/seaweedfs/seaweedfs/releases/download/3.64/linux_amd64_large_disk.tar.gz
RUN tar -xf linux_amd64_large_disk.tar.gz
RUN chmod +x weed
RUN mv weed /usr/bin/
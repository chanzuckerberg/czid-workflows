FROM golang:1.19.2-bullseye

WORKDIR /go/src/s3quilt
COPY . .

RUN apt-get update && apt-get install -y python3-pip

RUN go get golang.org/x/tools/cmd/goimports
RUN go install golang.org/x/tools/cmd/goimports
RUN go get github.com/go-python/gopy
RUN go install github.com/go-python/gopy

RUN pip3 install pybindgen setuptools wheel

RUN gopy build -vm=python3 --output s3quilt .
RUN cp -r s3quilt/ /usr/lib/python3/dist-packages/ && \
    cp __init__.py /usr/lib/python3/dist-packages/s3quilt/
